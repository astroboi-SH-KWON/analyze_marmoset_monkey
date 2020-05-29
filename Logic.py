from Bio import SeqIO

import Util
import LogicPrep

class Logics:
    def __init__(self):
        self.tmp = ""

    """
    checkSeqByChar : match sequences by char with rules
    :param
        dna_char :
        rule_char : rules with "A", "C", "G", "T", "U", "N", "R",...
    :return
        boolean
    """
    def checkSeqByChar(self,dna_char, rule_char):
        flag = False
        if rule_char == 'N':
            return True
        elif rule_char in 'ACGTU':
            if dna_char == rule_char:
                return True
        elif rule_char == 'R':
            if dna_char in 'AG':
                return True
        # elif rule_char == 'r':
        #     if dna_char in 'CT':
        #         return True
        """
        add more rules of "ACTGU"
        """
        return flag

    """
    match : match sequence with same length strings
    :param
        i : index of seq
        dna_seq : targeted DNA/RNA sequence 
        rule_str : rules with "ACGTU", "N", "R",...
    :return
        boolean
    """
    def match(self, i, dna_seq, rule_str):
        if len(dna_seq) == i:
            return True
        if self.checkSeqByChar(dna_seq[i], rule_str[i]):
            return self.match(i + 1, dna_seq, rule_str)
        else:
            return False

    def filter_out_by_ACGTU_rule(self, input_dict, window_idx_arr, rule_acgt_arr):
        result_dict = {}
        for chr_key, val_dict in input_dict.items():
            result_dict.update({chr_key: {}})
            for trnscrpt_id, vals_arr in val_dict.items():
                result_dict[chr_key].update({trnscrpt_id: [vals_arr[0]]})
                for trgt_idx in range(1, len(vals_arr)):
                    trgt_seq = vals_arr[trgt_idx][0][window_idx_arr[0] -1: window_idx_arr[1]]
                    for rule_acgt in rule_acgt_arr:
                        if rule_acgt in trgt_seq:
                            result_dict[chr_key][trnscrpt_id].append(vals_arr[trgt_idx])
                            break

        return result_dict


