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
FILE_NUM_LIST = ['X']
############### end setting env ################

def make_deep_pe_input():
    logic = Logic.Logics()
    util = Util.Utils()

    tmp_dict = logic.get_Deep_PE_input(REF_PATH + CDS_FILE, INIT_DEEP_PE)
    result_dict = logic.group_by_chromosome(tmp_dict, "primary_assembly:ASM275486v1:")

    util.make_Deep_PE_input_excel(WORK_DIR + "crab_eating/", result_dict, INIT_DEEP_PE)