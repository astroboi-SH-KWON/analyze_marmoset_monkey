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
WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200527/WORK_DIR/"

REF_PATH = "D:/000_WORK/000_reference_path/monkey/marmoset/"
CDS_FILE = "cds/Callithrix_jacchus.ASM275486v1.cds.all.fa"
DNA_FILE = "dna/Callithrix_jacchus.ASM275486v1.dna.nonchromosomal.fa"
ANNO_FILE = "genbank_anno/Callithrix_jacchus.ASM275486v1.100.nonchromosomal.dat"


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
#################### top N #####################
TOP_N = 10
TOP_N_ALL = 100
INIT_MERGE_BY_CHAR = [REF_PATH, CDS_FILE, A_or_C_IDX, ACTG_RULE, WORK_DIR, TOP_N]
INIT_MERGE_BY_ALL = [REF_PATH, CDS_FILE, A_or_C_IDX, ACTG_RULE, WORK_DIR, TOP_N_ALL]
############### end setting env ################

def sort_n_merge_by_all():
    logic = Logic.Logics()
    logic.sort_n_merge_by_all(INIT_MERGE_BY_ALL, INIT_BE)

def sort_n_merge_by_chr():
    logic = Logic.Logics()
    logic.sort_n_merge_by_chr(INIT_MERGE_BY_CHAR, INIT_BE)

def merge_cas9_abe_cbe():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    trgt_seq_dict = logic_prep.get_target_seq_with_clvg_site(REF_PATH + CDS_FILE, INIT_BE)
    chr_dict = logic_prep.target_seq_with_clvg_site_group_by_chromosome(trgt_seq_dict,
                                                                                      "primary_assembly:ASM275486v1:")
    a_c_dict = logic.filter_out_by_ACGTU_rule(chr_dict, A_or_C_IDX, ACTG_RULE)

    abe_score_dict = logic_prep.get_deep_base_ed_score(WORK_DIR + "deep_ABE/ABE_Efficiency.txt")
    cbe_score_dict = logic_prep.get_deep_base_ed_score(WORK_DIR + "deep_CBE/CBE_Efficiency.txt")
    cs9_score_dict = logic_prep.get_deep_cas9_tupl(WORK_DIR + "deep_cas_9/", "RANK_final_DeepCas9_Final.txt", "sample.txt")

    # util.make_merge_tab_txt(WORK_DIR + "marmoset_merge_abe_cbe_cas9", [a_c_dict, abe_score_dict, cbe_score_dict, cs9_score_dict], INIT_BE)
    util.make_merge_excel_by_chr(WORK_DIR + "merge_cas9_abe_cbe/marmoset_merge_abe_cbe_cas9", [a_c_dict, abe_score_dict, cbe_score_dict, cs9_score_dict], INIT_BE)

def make_deep_cas9_base_editor_input():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    trgt_seq_dict = logic_prep.get_target_seq_with_clvg_site(REF_PATH + CDS_FILE, INIT_BE)
    chr_dict = logic_prep.target_seq_with_clvg_site_group_by_chromosome(trgt_seq_dict,
                                                                        "primary_assembly:ASM275486v1:")
    a_c_dict = logic.filter_out_by_ACGTU_rule(chr_dict, A_or_C_IDX, ACTG_RULE)

    util.make_cas_off_finder_input(a_c_dict, INITIAL_CAS_OFF)
    util.make_deep_cas9_input(WORK_DIR + "deep_cas_9/sample", [a_c_dict], INIT_BE)

def make_deep_pe_input():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    trgt_seq_dict = logic_prep.get_target_seq(REF_PATH + CDS_FILE, INIT_DEEP_PE)
    result_dict = logic_prep.group_by_chromosome(trgt_seq_dict, "primary_assembly:ASM275486v1:")

    # util.make_Deep_PE_input_excel(WORK_DIR + "crab_eating/", result_dict, INIT_DEEP_PE)
    util.make_Deep_PE_input_tb_txt(WORK_DIR + "crab_eating/deep_pe_input_crab_eating", result_dict)

start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# make_deep_pe_input()
# make_deep_cas9_base_editor_input()
## merge_cas9_abe_cbe()
# sort_n_merge_by_chr()
sort_n_merge_by_all()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))