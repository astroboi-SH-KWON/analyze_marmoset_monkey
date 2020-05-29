from Bio import SeqIO
import glob
from time import clock
import re
import numpy as np

import Util
import Logic
import LogicPrep
import Valid

############### start to set env ################
WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200527/WORK_DIR/marmoset/"

REF_PATH = "D:/000_WORK/000_reference_path/monkey/marmoset/"
CDS_FILE = "cds/Callithrix_jacchus.ASM275486v1.cds.all.fa"
DNA_FILE = "dna/Callithrix_jacchus.ASM275486v1.dna.nonchromosomal.fa"
ANNO_FILE = "genbank_anno/Callithrix_jacchus.ASM275486v1.100.nonchromosomal.dat"

TEST_CHR = "D:/000_WORK/000_reference_path/monkey/chlorocebus_sabaeus/"
chr_file_name = "chromosome/Chlorocebus_sabaeus.ChlSab1.1.dna.chromosome.X.fa"
genbank_file_name = "genbank_anno/Chlorocebus_sabaeus.ChlSab1.1.99.chromosome.Y.dat"

FRONT_WIN_LEN = 4
gRNA_LEN = 20
PAM_SEQ = "NGG"
BACK_WIN_LEN = 20

INIT_DEEP_PE = [PAM_SEQ, FRONT_WIN_LEN, gRNA_LEN, BACK_WIN_LEN]
A_or_C_IDX = [4, 10]
ACTG_RULE = ['A', 'C']
############## make_deep_pe_input ##############
BE_BACK_WIN_LEN = 3
CLEAVAGE_SITE = 3
MAX_MISMATCH = 3
REF_SRV_PATH = "FASTA/marmoset"
INIT_BE = [PAM_SEQ, FRONT_WIN_LEN, gRNA_LEN, BE_BACK_WIN_LEN, CLEAVAGE_SITE]
INITIAL_CAS_OFF = ['NGG', gRNA_LEN, MAX_MISMATCH, 66, WORK_DIR + "CAS_OFF_FINDER/marmoset_monkey_off_", REF_SRV_PATH, INIT_BE]

############### end setting env ################

def merge():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    trgt_seq_dict = logic_prep.get_target_seq_with_clvg_site(REF_PATH + CDS_FILE, INIT_BE)
    chr_dict = logic_prep.target_seq_with_clvg_site_group_by_chromosome(trgt_seq_dict,
                                                                                      "primary_assembly:ASM275486v1:")
    a_c_dict = logic.filter_out_by_ACGTU_rule(chr_dict, A_or_C_IDX, ACTG_RULE)

    abe_score_dict = logic_prep.get_deep_base_ed_score(WORK_DIR + "deep_ABE/ABE_Efficiency.txt")
    print(len(abe_score_dict))
    cbe_score_dict = logic_prep.get_deep_base_ed_score(WORK_DIR + "deep_CBE/CBE_Efficiency.txt")
    print(len(cbe_score_dict))
    cs9_score_dict = logic_prep.get_deep_cas9_tupl(WORK_DIR + "deep_cas_9/", "RANK_final_DeepCas9_Final.txt", "sample.txt")
    print(len(cs9_score_dict))  # 2421639



def test():
    logic_prep = LogicPrep.LogicPreps()










start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
merge()
# test()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))







