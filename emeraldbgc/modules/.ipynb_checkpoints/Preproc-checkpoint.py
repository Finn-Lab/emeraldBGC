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
import subprocess

from emeraldbgc import _params

log = logging.getLogger(f"EMERALD.{__name__}")

from distutils.spawn import find_executable


class Preprocess:
    """External tools needed for emerald bgc detection"""

    def __init__(self, seq_file, meta, cpus, outdir):

        self.seq_file = seq_file
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
            ["-p", "meta"] if self.meta else []
        )

        outs, errs = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()

        if "Error:  Sequence must be 20000 characters" in errs.decode("utf8"):

            log.info(
                "Sequence must be 20000 characters when running Prodigal in normal mode. Trying -p meta"
            )

            cmd = ["prodigal", "-i", self.seq_file, "-a", outFaa] + (["-p", "meta"])

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

        with open(f"{base}.faa", "w") as h:

            for rec in recs:

                for f in rec.features:

                    if f.type == "CDS" and "translation" in f.qualifiers:

                        pid = (
                            f.qualifiers["protein_id"]
                            if "protein_id" in f.qualifiers
                            else f.qualifiers["locus_tag"]
                        )
                        seq = f.qualifiers["translation"][0]
                        h.write(f">{pid}\n{seq}\n")

        return os.path.abspath(f"{base}.emerald.faa")

    def check_fmt(self):
        """ Evaluate if input format is FNA or GBK"""
        with open(self.seq_file) as h:
            if h.readline()[0] == ">":
                self.fmt = "fna"
                log.info("FASTA sequence file detected")
            else:
                self.fmt = "gbk"
                log.info(f"FASTA sequence file NOT detected; trying {self.fmt}")

    def process_sequence(self):
        """ CDS prediction on sequence file"""
        self.check_fmt()
        if self.fmt == "fna":
            self.outFaa = self.runProdigal()
        elif self.fmt == "gbk":
            self.outFaa = self.gbkToProdigal()
        else:
            log.info("missing sequence file format")

        ih_f = self.runHmmScan()
        ip_f = self.runInterproscan()
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
        log.info(cmd)
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

        ip_exec = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "interproscan",
            "interproscan.sh",
        )

        if not find_executable(ip_exec):
            log.exception(f"{ip_exec} not found, run emrald_download_data")

        if not os.path.isfile(self.outFaa):
            log.exception(f"{self.outFaa} file not found")

        outGff = os.path.join(
            self.outdir, "{}.ip.tsv".format(os.path.basename(self.outFaa))
        )
        cmd = (
            [ip_exec, "-i", self.outFaa, "-o", outGff, "-f", "TSV"]
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