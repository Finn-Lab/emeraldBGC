import os
import sys
import filecmp
import pytest
import emeraldbgc

emrld_dir = os.path.dirname(os.path.abspath(emeraldbgc.__file__))
test_files_dir = os.path.join(emrld_dir, "test", "files")
sys.path.append(emrld_dir)

from modules.Preproc import Preprocess

def test_prodigal():
    pp = Preprocess(
        os.path.join(test_files_dir, "BGC0001472.fna"),
        None,
        False,
        1,
        os.path.join(test_files_dir, "Test") 
    )
    test_prodigal_file =  pp.runProdigal()
    assert filecmp.cmp(test_prodigal_file, os.path.join(test_files_dir, "BGC0001472.fna.prodigal.faa")) 
    
def test_gbk_transform():
    pp = Preprocess(
        os.path.join(test_files_dir, "BGC0001472.fna.prodigal.faa.gb"),
        None,
        False,
        1,
        os.path.join(test_files_dir, "Test") 
    )
    test_prodigal_file =  pp.gbkToProdigal()
    assert filecmp.cmp(test_prodigal_file, os.path.join(test_files_dir, "BGC0001472.fna.prodigal.faa.gb.prodigal.faa"))
     
def test_hmmscan_ih():
    pp = Preprocess(
        os.path.join(test_files_dir, None),
        None,
        False,
        1,
        os.path.join(test_files_dir, "Test") 
    )
    pp.outFaa = os.path.join(test_files_dir, "BGC0001472.fna.prodigal.faa")
    test_hmm_file =  pp.runHmmScan()
    assert filecmp.cmp(test_hmm_file, os.path.join(test_file_dir, "BGC0001472.fna.prodigal.faa.emerald.tsv"))
     
def test_ipscan_ih():
    
    if 'linux' not in sys.platform:
        assert True
    else: 
        pp = Preprocess(
            os.path.join(test_files_dir, None),
            None,
            False,
            1,
            os.path.join(test_files_dir, "Test") 
        )
        pp.outFaa = os.path.join(test_files_dir, "BGC0001472.fna.prodigal.faa")
        test_ip_file =  pp.runInterproscan()
        assert filecmp.cmp(test_ip_file, os.path.join(test_file_dir, "BGC0001472.fna.prodigal.faa.ip.tsv"))
