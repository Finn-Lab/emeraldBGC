import os
import sys
import filecmp
import pytest
import emeraldbgc

emrld_dir = os.path.dirname(os.path.abspath(emeraldbgc.__file__))
test_files_dir = os.path.join(emrld_dir, "test", "files")
sys.path.append(emrld_dir)

from modules.Preproc import BGCdection

def test_bgc_prediction():
    
    prodigal_file = os.path.join(test_files_dir, "BGC0001472.fna.prodigal.faa")
    ips_file = os.path.join(test_file_dir, "BGC0001472.fna.prodigal.faa.emerald.tsv")
    hmm_file = os.path.join(test_file_dir, "BGC0001472.fna.prodigal.faa.ip.tsv")

    ann = AnnotationFilesToEmerald()
    ann.transformIPS(ips_file)
    ann.transformEmeraldHmm(hmm_file)

    ann.transformCDSpredToCDScontigs(
            prodigal_file,
            preprocess.fmt)

    ann.buildMatrices()

    ann.predictAnn()

    ann.defineLooseClusters(score=args.score, g=args.greed)

    ann.predictType()

    outp = Outputs(
        ann,
        args.minimal_out,
        args.ref_b,
        args.outfile if args.outfile else f"{outdir}/{base}.emerald.full.gff",
    )
