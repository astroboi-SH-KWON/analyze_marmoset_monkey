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

TEST_CHR = "D:/000_WORK/000_reference_path/monkey/chlorocebus_sabaeus/"
chr_file_name = "chromosome/Chlorocebus_sabaeus.ChlSab1.1.dna.chromosome.X.fa"
genbank_file_name = "genbank_anno/Chlorocebus_sabaeus.ChlSab1.1.99.chromosome.Y.dat"

FRONT_WIN_LEN = 4
gRNA_LEN = 20
PAM_SEQ = "NGG"
BACK_WIN_LEN = 20

INIT_DEEP_PE = [PAM_SEQ, FRONT_WIN_LEN, gRNA_LEN, BACK_WIN_LEN]
FILE_NUM_LIST = ['X']
############### end setting env ################

def test():
    logic = Logic.Logics()
    util = Util.Utils()



def test2():
    tmp_dict = {}

    for idx in range(1, 20):
        FILE_NUM_LIST.append(str(idx))

    for seq_record in SeqIO.parse(REF_PATH + CDS_FILE, "fasta"):

        dscript = seq_record.description
        chrsm = dscript[dscript.index("primary_assembly:ASM275486v1:") + len("primary_assembly:ASM275486v1:"):].split(":")[0]

        if chrsm not in tmp_dict:
            tmp_dict.update({chrsm: seq_record.seq})
    with open(WORK_DIR + "marmoset_list.txt", "a") as f:

        for chrsm_key, val in tmp_dict.items():
            # if chrsm_key in FILE_NUM_LIST:
            #     continue
            f.writelines(str(chrsm_key) + "\n")
            f.writelines(str(val) + "\n")
            f.writelines(" \n")

start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
test()
# test2()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))







