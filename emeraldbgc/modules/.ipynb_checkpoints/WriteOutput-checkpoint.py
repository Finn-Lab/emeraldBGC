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

import json
import logging
import os
import pickle
from itertools import groupby

import numpy as np
import tensorflow as tf
from joblib import load

from emeraldbgc import __version__

log = logging.getLogger(f"EMERALD.{__name__}")


class Outputs:

    """ write results in file """

    def __init__(self, obj, minimal_output, refined_borders, outfile):

        self.annotation = obj
        self.minimal = minimal_output
        self.ref_b = refined_borders
        self.outfile = outfile
        self.gff3 = []
        self.clusts = []
        self.writeGff3()

    def writeGff3(self):

        log.info("build GFF3")
        log.info(self.outfile)
        type_code = {
            0: "Alkaloid",
            1: "NRP",
            2: "Polyketide",
            3: "RiPP",
            4: "Saccharide",
            5: "Terpene",
            6: "Other",
        }

        for contig in self.annotation.looseClst:

            if self.minimal != "True":
                for ix, f in enumerate(self.annotation.contigsDct[contig]):

                    ID, (start, end) = f
                    emrldProb = "{:.3f}".format(self.annotation.annResults[contig][ix])
                    self.gff3.append(
                        f"{contig}\tEMERALDv{__version__}\tCDS\t{start}\t{end}\t.\t.\t.\tID={ID};emerald_probability={emrldProb}"
                    )

            ct = 1
            for k, g in groupby(
                enumerate(
                    self.annotation.looseClst[contig][
                        : len(self.annotation.contigsDct[contig])
                    ]
                ),
                key=lambda x: x[1],
            ):

                if k == 0:
                    continue

                ID = f"{contig}_emrld_{ct}"
                ct += 1

                gg = list(list(zip(*g))[0])

                start, end = (
                    self.annotation.contigsDct[contig][gg[0]][1][0],
                    self.annotation.contigsDct[contig][gg[-1]][1][1],
                )

                edge = "{}{}".format(
                    1 if gg[0] == 0 else 0,
                    1 if gg[-1] == len(self.annotation.contigsDct[contig]) - 1 else 0,
                )

                typs = self.annotation.typesClst[contig][gg[0]]
                # typs = "nearest_MiBIG={};nearest_MiBIG_class={};nearest_MiBIG_jacardDistance={}".format(typs[0],typs[1],typs[2])
                self.gff3.append(
                    f"{contig}\tEMERALDv{__version__}\tCLUSTER\t{start}\t{end}\t.\t.\t.\tID={ID};{typs};partial={edge}"
                )

                if self.ref_b == "True":

                    ct2 = 1
                    for k2, g2 in groupby(
                        zip(gg, self.annotation.borderClst[contig][gg]),
                        key=lambda x: x[1],
                    ):
                        if k2 == 0:
                            continue
                        gg2 = list(list(zip(*g2))[0])
                        start, end = (
                            self.annotation.contigsDct[contig][gg2[0]][1][0],
                            self.annotation.contigsDct[contig][gg2[-1]][1][1],
                        )
                        edge = "{}{}".format(
                            1 if gg[0] == 0 else 0,
                            1
                            if gg[-1] == len(self.annotation.contigsDct[contig]) - 1
                            else 0,
                        )

                        self.gff3.append(
                            f"{contig}\tEMERALDv{__version__}\tCLUSTER_border\t{start}\t{end}\t.\t.\t.\tID={ID}_{ct2};{typs};partial={edge}"
                        )
                        ct2 += 1

        log.info(f"Writing output to file {self.outfile}")
        with open(self.outfile, "w") as h:

            h.write("##gff-version 3\n")

            key = lambda x: (
                x.split("\t")[0],
                int(x.split("\t")[3]),
                "Z" if x.split("\t")[2] == "CDS" else x.split("\t")[2],
            )
            for l in sorted(self.gff3, key=key):

                h.write(f"{l}\n")
