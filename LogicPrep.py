from Bio import SeqIO
import re

import Logic
import Util

class LogicPreps:

    def __init__(self):
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"

    def check_strnd(self, loc_str):
        if "-" in loc_str:
            return "-"
        return "+"

    """
    :param
        path : "~.cds.all.fa" file that has seq_record.description and CDS seq
        init_arr : [PAM_SEQ, FRONT_WIN_LEN, gRNA_LEN, BACK_WIN_LEN]
    """
    def get_target_seq(self, path, init_arr):
        logic = Logic.Logics()
        tmp_dict = {}

        pam_seq = init_arr[0]
        add_seq1_len = init_arr[1]
        spacer_len = init_arr[2]
        add_seq2_len = init_arr[3]
        pam_len = len(pam_seq)

        std_tot_len = add_seq1_len + spacer_len + pam_len + add_seq2_len
        for seq_record in SeqIO.parse(path, "fasta"):

            trncrpt_id = seq_record.id
            if trncrpt_id not in tmp_dict:
                tmp_dict.update({trncrpt_id: [seq_record.description]})

            tmp_p_str = ""
            for c in seq_record.seq:
                tmp_p_str = tmp_p_str + c.upper()

                if len(tmp_p_str) > std_tot_len:
                    tmp_p_str = tmp_p_str[-std_tot_len:]

                if len(tmp_p_str) == std_tot_len:
                    if 'N' not in tmp_p_str:
                        if logic.match(0, tmp_p_str[-(add_seq2_len + pam_len):-add_seq2_len], pam_seq):
                            tmp_dict[trncrpt_id].append(tmp_p_str)

        return tmp_dict

    def group_by_chromosome(self, data_dict, deli_str):
        result_dict = {}
        for trnscrpt_id, vals_arr in data_dict.items():
            dscript = vals_arr[0]
            chrsm = dscript[dscript.index(deli_str) + len(deli_str):].split(":")[0]

            if chrsm in result_dict:
                if trnscrpt_id not in result_dict[chrsm]:
                    result_dict[chrsm].update({trnscrpt_id: vals_arr})
            else:
                result_dict.update({chrsm: {trnscrpt_id: vals_arr}})

        return result_dict

    def get_target_seq_with_clvg_site(self, path, init_arr):
        logic = Logic.Logics()
        tmp_dict = {}

        pam_seq = init_arr[0]
        add_seq1_len = init_arr[1]
        spacer_len = init_arr[2]
        add_seq2_len = init_arr[3]
        clvg_site = init_arr[4]
        pam_len = len(pam_seq)

        std_tot_len = add_seq1_len + spacer_len + pam_len + add_seq2_len

        for seq_record in SeqIO.parse(path, "fasta"):
            tot_cds_len = len(seq_record.seq)

            trncrpt_id = seq_record.id
            if trncrpt_id not in tmp_dict:
                tmp_dict.update({trncrpt_id: [seq_record.description]})

            tmp_p_str = ""
            idx = 0
            for c in seq_record.seq:
                idx += 1
                tmp_p_str = tmp_p_str + c.upper()

                if len(tmp_p_str) > std_tot_len:
                    tmp_p_str = tmp_p_str[-std_tot_len:]

                if len(tmp_p_str) == std_tot_len:
                    if 'N' not in tmp_p_str:
                        if logic.match(0, tmp_p_str[-(add_seq2_len + pam_len):-add_seq2_len], pam_seq):
                            tmp_dict[trncrpt_id].append([tmp_p_str, ((idx - clvg_site - pam_len - add_seq2_len) / tot_cds_len) * 100])

        return tmp_dict

    def target_seq_with_clvg_site_group_by_chromosome(self, trgt_seq_dict, deli_str):
        result_dict = {}
        for trnscrpt_id, vals_arr in trgt_seq_dict.items():
            dscript = vals_arr[0]
            chrsm = dscript[dscript.index(deli_str) + len(deli_str):].split(":")[0]

            if chrsm in result_dict:
                if trnscrpt_id not in result_dict[chrsm]:
                    result_dict[chrsm].update({trnscrpt_id: vals_arr})
            else:
                result_dict.update({chrsm: {trnscrpt_id: vals_arr}})

        return result_dict

    def get_deep_base_ed_score(self, path):
        tmp_dict = {}
        with open(path, "r") as f:
            while True:
                tmp_line = f.readline().replace("\n","")
                if tmp_line == "":
                    break

                tmp_arr = tmp_line.split("\t")
                tmp_dict.update({tmp_arr[0]: tmp_arr[-1]})

        return tmp_dict

    def get_cs9_scre(self, scre_txt_path):
        tmp_tpl = ()
        with open(scre_txt_path, "r") as f:
            f.readline().replace("\n", "")
            f.readline().replace("\n", "")
            f.readline().replace("\n", "")
            f.readline().replace("\n", "")
            # make result as tuple
            tmp_tpl += eval(f.readline().replace("\n", ""))
        return tmp_tpl

    def get_deep_cas9_tupl(self, path, scre_txt_path, seq_txt_path):
        tmp_dict = {}
        tmp_tpl = self.get_cs9_scre(path + scre_txt_path)

        with open(path + seq_txt_path) as f:
            f.readline().replace("\n", "")
            idx = 0
            while True:
                tmp_line = f.readline().replace("\n", "")
                if tmp_line == "":
                    break

                tmp_arr = tmp_line.split("\t")
                tmp_dict.update({tmp_arr[1]: tmp_tpl[idx]})
                idx += 1

        return tmp_dict





