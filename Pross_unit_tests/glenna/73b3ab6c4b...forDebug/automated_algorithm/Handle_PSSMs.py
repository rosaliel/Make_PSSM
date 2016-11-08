import urllib,os,re, subprocess, sys
import linecache

from make_all_numbers_equal_in_length import *
from Align2seq_class import *
from path_variables import *

#The class gets a bunch of sub_PSSM files, each relevant for different residues in the protein chain and creates a
#single PSSM file representing them all (copies the relevant lines from each sub_PSSM file into a new file.

class Handle_PSSMs:

    def __init__(self):
        pass

    def Unite_PSSMs(self, all_loops_vec, number_correction_map, all_seq_pssm_name, path_to_all_seq_PSSM, path_to_sub_PSSMs):
        self.number_correction_map = number_correction_map
        print 'SELF CORRECTION MAP IS'
        print self.number_correction_map
        print 'ALL_LOOPS_VEC IS'
        print all_loops_vec
        self.path_to_all_seq_PSSM = path_to_all_seq_PSSM
        self.all_seq_pssm_name = all_seq_pssm_name
        self.path_to_sub_PSSMs = path_to_sub_PSSMs

        self.unified_pssm = open(self.path_to_all_seq_PSSM + '/' 'unified_pssm', 'a')

        #The loop prints the first 3 lines of the final pssm (common to all pssms).
        for i in range (0,3):
            line = linecache.getline(self.path_to_all_seq_PSSM + '/' + all_seq_pssm_name, i + 1)
            self.unified_pssm .write(line)

        #Used at the end of the next for loop
        previous_loop_no_gaps = [-1,-1]

        #The following loop is meant to deal with a special case where there is hypen in the MSA at the 1st position.
        #See also the same idea in the function define_name_for_sub_MSA (Make_sub_MSA_files class)
        for loop in all_loops_vec:
            counter = 0
            for i in range(loop[0], loop[1]):
                if self.number_correction_map[loop[0] + counter] == '-':
                    counter += 1
                else:
                    break

            loop_no_gaps = [self.number_correction_map[loop[0] + counter], self.number_correction_map[loop[1] -1]]
            pos_diff = loop_no_gaps[0] - previous_loop_no_gaps[1]

            #pos_diff = 1 means that there are no residues between the current sub_PSSM_file and the previous one.
            if pos_diff == 1:
                self.define_lines_to_copy_sub_PSSM(loop_no_gaps)

            #pos_diff > 1 means that there is at least 1 residue that is not represented in any sub_PSSM,
            #therefore its line is copied from the full pssm.
            elif pos_diff > 1:
                self.define_line_to_copy_full_pssm(loop_no_gaps[0], pos_diff)
                self.define_lines_to_copy_sub_PSSM(loop_no_gaps)

            else:
                print 'ERROR: sub_PSSMs NUMBERING IS WRONG!!!'
                print 'ERROR MAY RESULT FROM sub_PSSMs WITH OVERLAPPING POSITIONS'

            previous_loop_no_gaps[0] = loop_no_gaps[0]
            previous_loop_no_gaps[1] = loop_no_gaps[1]

        self.unified_pssm .close()

    def define_line_to_copy_full_pssm(self, next_loop_beg_pos, pos_diff_from_previous_loop):
        pssm_to_copy_path = self.path_to_all_seq_PSSM
        pssm_to_copy_from = self.all_seq_pssm_name

        ########################################################################################################
        # Explanation for the next 2 variables:                                                                #
        # 1."first_line_to_copy": the term in the brackets gets the 1st residue that is not in the previous     #
        #   loop/sub_MSA and also not yet in the next loop/sub_MSA. The +3 accounts for the fact that residue  #
        #   number 1 is printed in line 4 in a standard PSSM file. The +1 accounts for the fact that a file    #
        #   begins with line 1 while our indexing begins with 0.                                               #
        # 2."last_line_to_copy" follows similar logic to "first_line_to_copy". The extra +1 there is for the   #
        #   use of this variable in a for loop in the function self.define_lines_to_copy_sub_PSSM()           #                                                        #
        ########################################################################################################
        first_line_to_copy =  (next_loop_beg_pos - pos_diff_from_previous_loop + 1) + 3 + 1
        last_line_to_copy = (next_loop_beg_pos - 1) + 3 + 1 + 1

        self.copy_lines_from_pssm(pssm_to_copy_from, pssm_to_copy_path, first_line_to_copy, last_line_to_copy)

    def define_lines_to_copy_sub_PSSM(self, sub_PSSM_range):

        #The following loop is necessary to determine the sub_PSSM identity (see variable pssm_to_copy_from)
        range_no_gaps_str = []
        for i in range(0, len(sub_PSSM_range)):
            range_no_gaps_str.append(make_all_numbers_equal_in_length(sub_PSSM_range[i] + 1,\
                                     len(self.number_correction_map)))
        pssm_to_copy_from = 'sub_PSSM_pos_' + range_no_gaps_str[0] + 'to' + range_no_gaps_str[1]
        pssm_to_copy_path = self.path_to_sub_PSSMs
        first_line_to_copy =  sub_PSSM_range[0] + 3 + 1
        last_line_to_copy = sub_PSSM_range[1] + 3 + 1 + 1

        self.copy_lines_from_pssm(pssm_to_copy_from, pssm_to_copy_path, first_line_to_copy, last_line_to_copy)

    def copy_lines_from_pssm(self, pssm_to_copy_from, pssm_to_copy_path, first_line_to_copy, last_line_to_copy):

        for line_num in range (first_line_to_copy, last_line_to_copy):
            line = linecache.getline(pssm_to_copy_path + '/' + pssm_to_copy_from, line_num)
            self.unified_pssm.write(line)

    def remove_non_crystallized_res(self, pdb_seq, full_seq):
        counter = 0
        res_num_to_remove = ''
        res_id_to_remove = ''
        for res_num in range (0, len(pdb_seq)):
            counter += 1
            if pdb_seq[res_num] == '-':
                res_num_to_remove += str(counter) + ','
                res_id_to_remove += full_seq[res_num] + ','
        
	print pdb_seq
        print "RES_NUM_TO_REMOVE"
        print res_num_to_remove
        res_num_to_remove = '\"' + res_num_to_remove[:-1] + '\"'
        res_id_to_remove = '\"' + res_id_to_remove[:-1] + '\"'
        path_to_pssm = '\"' + self.path_to_all_seq_PSSM + '/' + 'unified_pssm' + '\"'

        os.system('sh ' + gen_path + 'rm_non_crystallized_res_from_pssm.sh '
                   + res_num_to_remove + ' '
                   + res_id_to_remove + ' '
                   + path_to_pssm + ' '
                   + '>tmp_cmd')
        os.system('chmod +x tmp_cmd')
        os.system('sh tmp_cmd >' + self.path_to_all_seq_PSSM + '/' + 'final_pssm_for_Rosetta')
        os.system('rm -rf tmp_cmd')


