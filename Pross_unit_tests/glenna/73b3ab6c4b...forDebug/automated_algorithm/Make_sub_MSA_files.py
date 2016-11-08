from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import urllib,os,re, subprocess, sys

#from PSSM_class import *

class Make_sub_MSA_files:

    def __init__(self, input_seq_name, input_seq_in_msa, loops_vec, regexp_vec, msa_name2seq_map):
        self.input_seq_name = input_seq_name
        self.input_seq_in_msa = input_seq_in_msa
        self.all_regexp = regexp_vec
        self.msa_name2seq_map = msa_name2seq_map
        self.all_loops = loops_vec
        self.sub_msa_dir = ''
        self.res_num_with_hyphens_2_res_num_no_hyphens = {}
        self.all_files_in_sub_MSA_dir = []

    def create_path_for_output(self, path):
        self.sub_msa_dir = path + '/' + 'sub_MSAs'
        if not os.path.exists(self.sub_msa_dir): os.makedirs(self.sub_msa_dir)

    def iter_over_regexps_to_find_matching_seqs_in_msa(self):
        for loop in range(0,len(self.all_loops)):
            seq_per_file_counter = 0
            regexp = re.compile('^'+(self.all_regexp[loop])+'$')

            matching_seqs = {}
            for key, val in self.msa_name2seq_map.items():
                substr_to_compare = self.msa_name2seq_map[key][self.all_loops[loop][0]:self.all_loops[loop][1]]
                if regexp.match(substr_to_compare):
                    matching_seqs[key] = val
                    seq_per_file_counter += 1

            self.make_sub_msa_files_fa(loop, matching_seqs)   #Easier for automatic processing
            self.make_sub_msa_files_clw(loop, matching_seqs)  #Easier for manual inspection

    def res_renum_ignore_gaps(self):
        counter = 0
        for res_num in range (0, len(self.input_seq_in_msa)):
            key = res_num
            if self.input_seq_in_msa[res_num] == '-':
                counter += 1
                self.res_num_with_hyphens_2_res_num_no_hyphens[key] = '-'
            else:
                self.res_num_with_hyphens_2_res_num_no_hyphens[key] = res_num - counter

    def make_sub_msa_files_fa(self, loop, matching_seqs):
        path = self.get_sub_MSA_path() + '/'
        msa_file_name = self.define_name_for_sub_MSA(loop)
        msa_file = open(path + msa_file_name, 'a')

        #Taking care of printing the Input_seq first (all blocks besides the last)
        msa_file.write('>' + self.input_seq_name + '\n')
        msa_file.write(self.msa_name2seq_map[self.input_seq_name] + '\n')

        for key_, val_ in matching_seqs.items():
            if key_ != self.input_seq_name: #Input was already printed
                msa_file.write('>' + key_ + '\n')
                msa_file.write(val_ + '\n')

        msa_file.close()


    def make_sub_msa_files_clw(self, loop, matching_seqs):

        aa_per_line = 60
        if len(self.input_seq_in_msa) % aa_per_line == 0:
            num_of_blocks_in_msa = len(self.input_seq_in_msa)/aa_per_line
        else:
            num_of_blocks_in_msa = len(self.input_seq_in_msa)/aa_per_line + 1

        path = self.get_sub_MSA_path() + '/'
        msa_file_name = self.define_name_for_sub_MSA(loop) + '_clw'
        msa_file = open(path + msa_file_name, 'a')
        msa_file.write('CLUSTAL\n\n')

        for i in range(0, num_of_blocks_in_msa):
            #The 'if' takes care of printing all the blocks besides the last one (as it might have a diff length)
            if i < num_of_blocks_in_msa - 1:
                line_range = (i+1)*aa_per_line
            else:
                line_range = len(self.input_seq_in_msa)

            #Taking care of printing the Input_seq first (all blocks besides the last)
            msa_file.write('{:<63}'.format(self.input_seq_name) + '{:<63}'.format(self.msa_name2seq_map\
                          [self.input_seq_name][i*aa_per_line:line_range]))
            msa_file.write('\n')

            for key_, val_ in matching_seqs.items():
                if key_ != self.input_seq_name: #Input was already printed
                    #if key_ is too long the msa becomes unreadable so I trim it [0:57].
                    msa_file.write('{:<63}'.format(key_[0:57]) + '{:<63}'.format(val_[i*aa_per_line:line_range]))
                    msa_file.write('\n')

            msa_file.write('\n')

        msa_file.close()

        #When the MSA has only 1 sequence the PSSM module does not work on clustal so I convert to fasta format
        #The original file will get the ending "_clustal" and the fasta file will get the original name.
        # if len(matching_seqs.items()) == 1:
        #     os.system('mv ' + path + msa_file_name + ' ' + path + msa_file_name + '_clustal')
        #     os.system('convert_clustal2fasta.sh ' + path + msa_file_name + '_clustal' '>' + path + msa_file_name)

    def define_name_for_sub_MSA(self, loop):

	#The next loop deals with the special case of hypen in the MSA's 1st position
        #See also the same idea in the function Unite_PSSMs (Handle_PSSMs class)
        counter_b = 0
        for i in range (self.all_loops[loop][0], self.all_loops[loop][1]):
            if self.res_num_with_hyphens_2_res_num_no_hyphens[self.all_loops[loop][0] + counter_b] == '-':
                counter_b += 1
            else:
                break
        #The next loop deals with the special case of hypen in the MSA's last position.
        counter_e = 0
        for i in range (self.all_loops[loop][1]-1, self.all_loops[loop][0], -1):
            if self.res_num_with_hyphens_2_res_num_no_hyphens[self.all_loops[loop][1]-1 + counter_e] == '-':
                counter_e -= 1 #zero or negative number
            else:
                break

        begin_pos = str(self.res_num_with_hyphens_2_res_num_no_hyphens[self.all_loops[loop][0] + counter_b] + 1)
        end_pos = str(self.res_num_with_hyphens_2_res_num_no_hyphens[self.all_loops[loop][1]-1 + counter_e] + 1)
        begin_pos = self.make_all_numbers_equal_in_length(begin_pos, len(self.input_seq_in_msa))
        end_pos = self.make_all_numbers_equal_in_length(end_pos, len(self.input_seq_in_msa))

        msa_file_name = 'sub_MSA_for_pos_' + begin_pos + 'to' + end_pos

        self.all_files_in_sub_MSA_dir.append(msa_file_name)

        return msa_file_name

    def make_all_numbers_equal_in_length(self, num_to_handle, longest_num):
        num_to_str = str(num_to_handle)
        num_to_str_len = len(num_to_str)
        longest_num_len = len(str(longest_num))
        zeros_to_add_before = ''

        for i in range (num_to_str_len, longest_num_len):
            zeros_to_add_before += '0'

        num_to_str = zeros_to_add_before + num_to_str

        return num_to_str

    def get_sub_MSA_path(self):
        return self.sub_msa_dir

    def get_all_sub_MSA_file_names(self):
        return self.all_files_in_sub_MSA_dir

    def get_res_num_with_hyphens_2_res_num_no_hyphens_map(self):
        return self.res_num_with_hyphens_2_res_num_no_hyphens
