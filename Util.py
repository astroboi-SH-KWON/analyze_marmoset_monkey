from os import listdir
from os.path import isfile, join
import pandas as pd
import openpyxl
from time import clock
import random
import math

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    def read_excel_2_dataframe(self, path):
        return pd.read_excel(path + self.ext_xlsx)

    def make_Deep_PE_input_excel(self, path, data_dict, init_arr):
        pam_seq = init_arr[0]
        add_seq1_len = init_arr[1]
        spacer_len = init_arr[2]
        add_seq2_len = init_arr[3]
        pam_len = len(pam_seq)

        for chr_key, val_dict in data_dict.items():

            workbook = openpyxl.Workbook()
            sheet = workbook.active

            row = 1
            sheet.cell(row=row, column=1, value="INDEX")
            sheet.cell(row=row, column=2, value='Target gene name')
            sheet.cell(row=row, column=3, value='Ensembl transcript ID')
            sheet.cell(row=row, column=4, value='Target sequence')
            # sheet.cell(row=row, column=5, value='gRNA binding')
            # sheet.cell(row=row, column=6, value='PAM')

            for trnscrpt_id, vals_arr in val_dict.items():
                gene_name = ""
                dscript = vals_arr[0]

                if "gene_symbol:" in dscript:
                    gene_name = dscript[dscript.index("gene_symbol:") + len("gene_symbol:"):].split(" ")[0]
                for trgt_idx in range(1, len(vals_arr)):
                    row += 1
                    sheet.cell(row=row, column=1, value=str(row - 1))
                    sheet.cell(row=row, column=2, value=gene_name)
                    sheet.cell(row=row, column=3, value=trnscrpt_id)
                    sheet.cell(row=row, column=4, value=vals_arr[trgt_idx])
                    # sheet.cell(row=row, column=5, value=vals_arr[trgt_idx][add_seq1_len: - pam_len - add_seq2_len])
                    # sheet.cell(row=row, column=6, value=vals_arr[trgt_idx][add_seq1_len + spacer_len: - add_seq2_len])

            # workbook.save(filename=path + "deep_pe_input_chr_" + str(chr_key) + "_" + str(clock()) + self.ext_xlsx)
            workbook.save(filename=path + "deep_pe_input_chr_" + str(chr_key) + self.ext_xlsx)

    def make_cas_off_finder_input(self, result_dict, init_arr):
        pam_seq = init_arr[0]
        spacr_len = init_arr[1]
        mis_mtch_num = init_arr[2]
        sub_set_num = init_arr[3]
        path = init_arr[4]
        ref_srv_path = init_arr[5]
        init_be_list = init_arr[6]

        seq_set_pool = set()
        for chrm_key, val_dict in result_dict.items():
            for trnscrpt_id, val_arr in val_dict.items():
                val_arr_len = len(val_arr)
                for idx in range(1, val_arr_len):
                    seq_set_pool.add(val_arr[idx][0][init_be_list[1]:init_be_list[1] + init_be_list[2]])

        raw_each_sub_len = len(seq_set_pool) / sub_set_num
        each_sub_len = int(raw_each_sub_len)
        # * 1000) / 1000 : raw_each_sub_len 값이 반올림되는 경우의 error 방지
        raw_last_sub_len_to_add = int(((raw_each_sub_len - each_sub_len) * sub_set_num) * 1000) / 1000
        last_sub_len_to_add = math.ceil(raw_last_sub_len_to_add)
        print("total len : " + str(len(seq_set_pool)) + ", raw_each_sub_len : " + str(raw_each_sub_len))
        print("raw_last_sub_len_to_add : " + str(raw_last_sub_len_to_add) + " , last_sub_len_to_add : " + str(
            last_sub_len_to_add))

        for fname_idx in range(sub_set_num):
            if fname_idx == sub_set_num - 1:
                each_sub_len += last_sub_len_to_add
            tmp_sub_set = set(random.sample(seq_set_pool, each_sub_len))
            seq_set_pool -= tmp_sub_set
            print("tmp_sub_set_" + str(fname_idx) + " len : " + str(len(tmp_sub_set)))
            print("rest total len : " + str(len(seq_set_pool)))
            print(" ")
            with open(path + str(fname_idx) + self.ext_txt, 'a') as f:
                f.write(ref_srv_path + "\n")
                f.write("N" * spacr_len + pam_seq + "\n")
                for data_str in tmp_sub_set:
                    f.write(data_str + "NNN " + str(mis_mtch_num) + "\n")

    def make_Deep_PE_input_tb_txt(self, path, data_dict):
        row = 0
        with open(path + self.ext_txt, 'a') as f:
            f.write("INDEX\tTarget gene name\tEnsembl transcript ID\tTarget sequence\n")
            for chr_key, val_dict in data_dict.items():
                for trnscrpt_id, vals_arr in val_dict.items():
                    gene_name = ""
                    dscript = vals_arr[0]

                    # get gene name
                    if "gene_symbol:" in dscript:
                        gene_name = dscript[dscript.index("gene_symbol:") + len("gene_symbol:"):].split(" ")[0]
                    for trgt_idx in range(1, len(vals_arr)):
                        row += 1
                        f.write(str(row) + "\t" + gene_name + "\t" + trnscrpt_id + "\t" + vals_arr[trgt_idx] + "\n")

    def make_deep_cas9_input(self, path, dict_arr, init_arr):
        pam_seq = init_arr[0]
        add_seq1_len = init_arr[1]
        spacer_len = init_arr[2]
        add_seq2_len = init_arr[3]
        pam_len = len(pam_seq)

        # del duplicates
        tmp_set = set()
        for data_dict in dict_arr:
            for chr_key, val_dict in data_dict.items():
                for trnscrpt_id, vals_arr in val_dict.items():
                    for trgt_idx in range(1, len(vals_arr)):
                        tmp_set.add(vals_arr[trgt_idx][0])

        row = 0
        with open(path + self.ext_txt, 'a') as f:
            f.write("Target number\t30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp)\n")
            for seq_str in tmp_set:
                row += 1
                f.write(str(row) + "\t" + seq_str + "\n")

