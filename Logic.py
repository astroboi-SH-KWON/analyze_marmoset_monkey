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

    def get_Deep_PE_input(self, path, init_arr):
        tmp_dict = {}

        pam_seq = init_arr[0]
        add_seq1_len = init_arr[1]
        spacer_len = init_arr[2]
        add_seq2_len = init_arr[3]
        pam_len = len(pam_seq)

        std_tot_len = add_seq1_len + spacer_len + pam_len + add_seq2_len
        for seq_record in SeqIO.parse(path, "fasta"):

            # print(seq_record)
            trncrpt_id = seq_record.id
            if trncrpt_id not in tmp_dict:
                tmp_dict.update({trncrpt_id: [seq_record.description]})

            tmp_p_str = ""
            # cds_seq_len = len(seq_record)
            for c in seq_record.seq:
                tmp_p_str = tmp_p_str + c.upper()

                if len(tmp_p_str) > std_tot_len:
                    tmp_p_str = tmp_p_str[-std_tot_len:]

                if len(tmp_p_str) == std_tot_len:
                    if 'N' not in tmp_p_str:
                        if self.match(0, tmp_p_str[-(add_seq2_len + pam_len):-add_seq2_len], pam_seq):
                            tmp_dict[trncrpt_id].append(tmp_p_str)

        return tmp_dict
