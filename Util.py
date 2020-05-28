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

    def make_Deep_PE_input_excel(self, path, data_dict, init_arr):
        pam_seq = init_arr[0]
        add_seq1_len = init_arr[1]
        spacer_len = init_arr[2]
        add_seq2_len = init_arr[3]
        pam_len = len(pam_seq)

        workbook = openpyxl.Workbook()
        sheet = workbook.active

        row = 1
        sheet.cell(row=row, column=1, value="INDEX")
        sheet.cell(row=row, column=2, value='Target gene name')
        sheet.cell(row=row, column=3, value='Ensembl transcript ID')
        sheet.cell(row=row, column=4, value='Target sequence')
        # sheet.cell(row=row, column=5, value='gRNA binding')
        # sheet.cell(row=row, column=6, value='PAM')

        for trnscrpt_id, vals_arr in data_dict.items():
            gene_name = ""
            if "gene_symbol:" in vals_arr[0]:
                gene_name = vals_arr[0][vals_arr[0].index("gene_symbol:") + len("gene_symbol:"):].split(" ")[0]
            for trgt_idx in range(1, len(vals_arr)):
                row += 1
                sheet.cell(row=row, column=1, value=str(row - 1))
                sheet.cell(row=row, column=2, value=gene_name)
                sheet.cell(row=row, column=3, value=trnscrpt_id)
                sheet.cell(row=row, column=4, value=vals_arr[trgt_idx])
                # sheet.cell(row=row, column=5, value=vals_arr[trgt_idx][add_seq1_len: - pam_len - add_seq2_len])
                # sheet.cell(row=row, column=6, value=vals_arr[trgt_idx][add_seq1_len + spacer_len: - add_seq2_len])


        workbook.save(filename=path + "_deep_pe_input_" + str(clock()) + self.ext_xlsx)
