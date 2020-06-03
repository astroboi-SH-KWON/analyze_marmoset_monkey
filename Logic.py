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

    """
    :return
        result_dict = {chr_key: {trnscrpt_id: [
                                    'description'
                                    , ['seq', clvg_site]
                                    , ['seq', clvg_site] ...
                                    ]
                                }
                        'NTIC01035965.1': {'ENSCJAT00000001565.3': [
                                                'ENSCJAT00000001565.3 cds primary_assembly:ASM275486v1:NTIC01035965.1:203:857:1 gene:ENSCJAG00000000864.4 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:TMEM190 description:transmembrane protein 190 [Source:HGNC Symbol;Acc:HGNC:29632]'
                                                , ['GGGGCGCTGTGAAGAAGCCAGCACTGGTGG', 8.384458077709612]
                                                , ['GCGCTGTGAAGAAGCCAGCACTGGTGGTTG', 8.997955010224949]
                                                ...
                                                ]
                        }
    """
    def filter_out_by_ACGTU_rule(self, input_dict, window_idx_arr, rule_acgt_arr):
        result_dict = {}
        for chr_key, val_dict in input_dict.items():
            result_dict.update({chr_key: {}})
            for trnscrpt_id, vals_arr in val_dict.items():
                result_dict[chr_key].update({trnscrpt_id: [vals_arr[0]]})
                for trgt_idx in range(1, len(vals_arr)):
                    trgt_seq = vals_arr[trgt_idx][0][window_idx_arr[0] - 1: window_idx_arr[1]]
                    for rule_acgt in rule_acgt_arr:
                        if rule_acgt in trgt_seq:
                            result_dict[chr_key][trnscrpt_id].append(vals_arr[trgt_idx])
                            break

        return result_dict

    def sort_n_merge_by_chr(self, init_merge, init_be):
        ref_path = init_merge[0]
        cdf_file = init_merge[1]
        a_or_c_idx = init_merge[2]
        a_c_rule = init_merge[3]
        work_dir = init_merge[4]
        top_n = init_merge[5]

        logic_prep = LogicPrep.LogicPreps()
        util = Util.Utils()

        trgt_seq_dict = logic_prep.get_target_seq_with_clvg_site(ref_path + cdf_file, init_be)
        chr_dict = logic_prep.target_seq_with_clvg_site_group_by_chromosome(trgt_seq_dict,
                                                                            "primary_assembly:ASM275486v1:")
        a_c_dict = self.filter_out_by_ACGTU_rule(chr_dict, a_or_c_idx, a_c_rule)

        abe_score_dict = logic_prep.get_deep_base_ed_score(work_dir + "deep_ABE/ABE_Efficiency.txt")
        cbe_score_dict = logic_prep.get_deep_base_ed_score(work_dir + "deep_CBE/CBE_Efficiency.txt")
        cs9_score_dict = logic_prep.get_deep_cas9_tupl(work_dir + "deep_cas_9/", "RANK_final_DeepCas9_Final.txt",
                                                       "sample.txt")

        top_n_abe_list = []
        top_n_cbe_list = []
        for chr_key, trnscrpt_list in a_c_dict.items():
            result_list = []
            result_list = logic_prep.merge_cas9_abe_cbe_to_list(chr_key, [trnscrpt_list, abe_score_dict, cbe_score_dict,
                                                                          cs9_score_dict], result_list)

            sort_by_abe_list = logic_prep.sort_by_idx_element(result_list, -2, [])
            sort_by_cbe_list = logic_prep.sort_by_idx_element(result_list, -1, [])

            """
            # extend TOP N lists to (top_n_abe_list, top_n_cbe_list)
            it needs filter out same context seq in different trnscrpt
            """
            top_n_abe_list.extend(sort_by_abe_list[:top_n])
            top_n_cbe_list.extend(sort_by_cbe_list[:top_n])

        util.make_excel_after_sorting(work_dir + "merge_cas9_abe_cbe_top_N/merge_by_ABE_top_" + str(top_n),
                                      top_n_abe_list, init_be)
        util.make_excel_after_sorting(work_dir + "merge_cas9_abe_cbe_top_N/merge_by_CBE_top_" + str(top_n),
                                      top_n_cbe_list, init_be)

    def sort_n_merge_by_all(self, init_merge, init_be):
        ref_path = init_merge[0]
        cdf_file = init_merge[1]
        a_or_c_idx = init_merge[2]
        a_c_rule = init_merge[3]
        work_dir = init_merge[4]
        top_n = init_merge[5]

        logic_prep = LogicPrep.LogicPreps()
        util = Util.Utils()

        trgt_seq_dict = logic_prep.get_target_seq_with_clvg_site(ref_path + cdf_file, init_be)
        chr_dict = logic_prep.target_seq_with_clvg_site_group_by_chromosome(trgt_seq_dict,
                                                                            "primary_assembly:ASM275486v1:")
        a_c_dict = self.filter_out_by_ACGTU_rule(chr_dict, a_or_c_idx, a_c_rule)

        abe_score_dict = logic_prep.get_deep_base_ed_score(work_dir + "deep_ABE/ABE_Efficiency.txt")
        cbe_score_dict = logic_prep.get_deep_base_ed_score(work_dir + "deep_CBE/CBE_Efficiency.txt")
        cs9_score_dict = logic_prep.get_deep_cas9_tupl(work_dir + "deep_cas_9/", "RANK_final_DeepCas9_Final.txt",
                                                       "sample.txt")

        result_list = []
        for chr_key, trnscrpt_list in a_c_dict.items():
            result_list = logic_prep.merge_cas9_abe_cbe_to_list(chr_key, [trnscrpt_list, abe_score_dict, cbe_score_dict,
                                                                          cs9_score_dict], result_list)

        """
        # sorting then get TOP_N_ALL
        it needs filter out same context seq in different trnscrpt
        """
        sort_by_abe_list = logic_prep.sort_by_idx_element(result_list, -2, [])[:top_n]
        sort_by_cbe_list = logic_prep.sort_by_idx_element(result_list, -1, [])[:top_n]

        util.make_excel_after_sorting(work_dir + "merge_cas9_abe_cbe_top_N/merge_by_ABE_top_" + str(top_n) + "_all",
                                      sort_by_abe_list, init_be)
        util.make_excel_after_sorting(work_dir + "merge_cas9_abe_cbe_top_N/merge_by_CBE_top_" + str(top_n) + "_all",
                                      sort_by_cbe_list, init_be)

