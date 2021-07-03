# Copyright 2021 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import os
import sys
import subprocess

from emeraldbgc import _params
log = logging.getLogger(f"EMERALD.{__name__}")

from distutils.spawn import find_executable


class Preprocess:
    """External tools needed for emerald bgc detection"""

    def __init__(self, seq_file, ip_file, meta, cpus, outdir):

        self.seq_file = seq_file
        self.ip_file = ip_file
        self.meta = meta
        self.cpus = int(cpus)
        self.outdir = outdir if outdir else "temp"

    def runProdigal(self):
        """Predict genes using prodigal"""
        log.info("Progial gene prediction...")

        if not find_executable("prodigal"):
            log.exception("Parodigal is not installed or in PATH")

        if not os.path.isfile(self.seq_file):
            log.exception(f"{self.seq_file} file not found")

        outFaa = os.path.join(
            self.outdir, "{}.prodigal.faa".format(os.path.basename(self.seq_file))
        )
        cmd = ["prodigal", "-i", self.seq_file, "-a", outFaa] + (
            ["-p", "meta"] if self.meta == "True" else []
        )

        outs, errs = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()

        if "Error:  Sequence must be 20000 characters" in errs.decode("utf8"):

            log.info(
                "Sequence must be 20000 characters when running Prodigal in normal mode. Trying -p meta"
            )

            cmd = ["prodigal", "-i", self.seq_file, "-a", outFaa, "-m"] + (["-p", "meta"])

            outs, errs = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()

        log.info(errs.decode("utf8"))

        """ Remove asterixs from faa
            """
        log.info("Removing asterix from prodigal faa")
        with open(outFaa, "r") as h:
            noAstFaa = [x.replace("*", "") for x in h]
        with open(outFaa, "w") as h:
            for l in noAstFaa:
                h.write(f"{l}")

        return os.path.abspath(outFaa)

    def gbkToProdigal(self):
        """Transform gbk to faa. This enables preprocessing with sequences"""
        log.info("write gbk as faa")
        from Bio import SeqIO

        recs = list(SeqIO.parse(open(self.seq_file, "r"), "gb"))

        base = os.path.basename(self.seq_file)

        outFaa = os.path.join(
            self.outdir, "{}.prodigal.faa".format(os.path.basename(self.seq_file))
        )

        with open(outFaa, "w") as h:

            for rec in recs:
                ct = 0
                for f in rec.features:

                    
                    if f.type == "CDS" and "translation" in f.qualifiers:
                        ct += 1
                        pid = (
                            f.qualifiers["protein_id"][0]
                            if "protein_id" in f.qualifiers
                            else f.qualifiers["locus_tag"][0]
                        )
                        seq = f.qualifiers["translation"][0]
                        h.write(f">{pid}\n{seq}\n")
                        
                if ct == 0:
                    log.info("{} CDS found with translation in {}".format(ct, rec.name) )

        return os.path.abspath(outFaa)

    def check_fmt(self):
        """ Evaluate if input format is FNA or GBK"""
        with open(self.seq_file) as h:
            if h.readline()[0] == ">":
                self.fmt = "fna"
                log.info("FASTA sequence file detected")
            else:
                self.fmt = "gbk"
                log.info(f"FASTA sequence file NOT detected; trying GBK")

    def process_sequence(self):
        """ CDS prediction on sequence file"""
        if 'linux' not in sys.platform and self.ip_file == None:
            log.info("internal run of InterProScan only available for Linux OS. Make sure to use --ip-file option")
            raise ValueError('Non Linux OS must be run with --ip-file option')

        self.check_fmt()
        if self.fmt == "fna":
            self.outFaa = self.runProdigal()
        elif self.fmt == "gbk":
            self.outFaa = self.gbkToProdigal()
        else:
            log.info("missing sequence file format")

        ih_f = self.runHmmScan()
        ip_f = self.runInterproscan() if self.ip_file == None else self.ip_file
        return self.outFaa, ip_f, ih_f

    def runHmmScan(self):
        """annotate functionally with emerald hmm library and hmmScan"""
        log.info("emerald functional annotation...")
        hmmLib = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "models",
            "hmm_lib",
            "emerald.hmm",
        )

        if not find_executable("hmmscan"):
            log.exception("hmmscan is not installed or in PATH")

        if not os.path.isfile(self.outFaa):
            log.exception(f"{self.outFaa} file not found")

        if not os.path.isfile(hmmLib):
            log.exception(f"{hmmLib} file not found")

        outTsv = os.path.join(
            self.outdir, "{}.emerald.tsv".format(os.path.basename(self.outFaa))
        )

        cmd = (
            ["hmmscan", "--domtblout", outTsv, "--cut_ga"]
            + (["--cpu", str(self.cpus)] if self.cpus else [])
            + ([hmmLib, self.outFaa])
        )
        log.info(" ".join(cmd))
        try:
            outs, errs = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            log.info(errs.decode("utf8"))
        except subprocess.CalledProcessError as err:
            log.exception(err.output)

        return os.path.abspath(outTsv)

    def runInterproscan(self):
        """annotate functionally with InterproScan"""
        log.info("InterProScan")

        if not find_executable('interproscan.sh'):
            log.exception(f"interproscan.sh not found, only available for Linux OS")

        if not os.path.isfile(self.outFaa):
            log.exception(f"self.outFaa file not found")

        outGff = os.path.join(
            self.outdir, "{}.ip.tsv".format(os.path.basename(self.outFaa))
        )
        cmd = (
            ["interproscan.sh", "-i", self.outFaa, "-o", outGff, "-f", "TSV"]
            + (["-appl", ",".join(_params["ip_an"])] if _params["ip_an"] else [])
            + (["-cpu", str(self.cpus)] if self.cpus else [])
        )
        log.info(" ".join(cmd))

        outs, errs = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        log.info(outs.decode("utf8"))
        log.info(errs.decode("utf8"))

        return os.path.abspath(outGff)
