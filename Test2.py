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

MUT_FILE = "Mutation_summary"
WINDOW_SIZE = 1
MAX_SEQ_LEN = 9

# INITIAL_MAIN = [CHR_PATH, MAX_SEQ_LEN, WINDOW_SIZE]
############### end setting env ################

def test():
    idx = 1
    flag = False
    for seq_record in SeqIO.parse(REF_PATH + ANNO_FILE, "genbank"):
        if flag:
            break
        print("id : " + seq_record.id)
        # ['data_file_division', 'date', 'accessions', 'sequence_version', 'keywords', 'source', 'organism', 'taxonomy', 'comment']
        print(seq_record.features)
        seq_feature = seq_record.features
        seq_f = seq_feature.__iter__()
        while True:
            try:
                gene = next(seq_f)
                print("")
                print(gene)
                idx += 1
                print("idx  ::::::::::: " + str(idx))
                if idx > 20:
                    flag = True
                    break
                # if gene.type == 'gene':
                #     print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                #     if 'note' in gene.qualifiers:
                #         print("Description : " + gene.qualifiers['note'][0])
                #     if 'locus_tag' in gene.qualifiers:
                #         print("Target gene name : " + gene.qualifiers['locus_tag'][0])
                #
                # if gene.type == 'CDS':
                #     print(gene.location)
                #     if 'gene' in gene.qualifiers:
                #         print("Ensembl Gene ID : " + gene.qualifiers['gene'][0])
                #     if 'note' in gene.qualifiers:
                #         print("Ensembl transcript ID : " + gene.qualifiers['note'][0])
            except StopIteration:
                print("end file")
                break

def read_genbank():
    for seq_record in SeqIO.parse(TEST_CHR + genbank_file_name, "genbank"):
        print("id : " + seq_record.id)
        # ['data_file_division', 'date', 'accessions', 'sequence_version', 'keywords', 'source', 'organism', 'taxonomy', 'comment']
        print(seq_record.features)
        seq_feature = seq_record.features
        seq_f = seq_feature.__iter__()
        while True:
            try:
                gene = next(seq_f)
                # print(gene)
                if gene.type == 'gene':
                    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                    if 'note' in gene.qualifiers:
                        print("Description : " + gene.qualifiers['note'][0])
                    if 'locus_tag' in gene.qualifiers:
                        print("Target gene name : " + gene.qualifiers['locus_tag'][0])

                if gene.type == 'CDS':
                    print(gene.location)
                    if 'gene' in gene.qualifiers:
                        print("Ensembl Gene ID : " + gene.qualifiers['gene'][0])
                    if 'note' in gene.qualifiers:
                        print("Ensembl transcript ID : " + gene.qualifiers['note'][0])
            except StopIteration:
                print("end file")
                break

def read_FASTA_all_at_once():
    idx = 1
    for seq_record in SeqIO.parse(REF_PATH + CDS_FILE, "fasta"):
        print("id : " + seq_record.id)
        print("name : " + seq_record.name)
        print("description : " + seq_record.description)
        print(repr(seq_record.seq))
        print(len(seq_record))
        idx +=1
        if idx == 10:
            break

def read_FASTA_head():
    for seq_record in SeqIO.parse(REF_PATH + DNA_FILE, "fasta"):
        print("id : " + seq_record.id)
        print("description : " + seq_record.description)

def just_read():
    idx = 1
    with open(REF_PATH + DNA_FILE, "r") as f:
        while True:
            idx += 1
            if idx == 20:
                break
            str_line = f.readline()
            print(str_line)
            if str_line == "":
                break

start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# test()
# just_read()
read_FASTA_all_at_once()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))









