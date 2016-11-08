from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import urllib,os,re, subprocess, sys

from path_variables import *

class Handle_dssp:

    def __init__(self, pdb_file_name, path_to_dir_of_pdb):
        self.pdb_file_name = pdb_file_name
        self.path_to_dir = path_to_dir_of_pdb

        self.protocol = protocols_path + 'dssp.xml'
        self.flags = flags_path + 'flags_for_dssp'

        self.ori_dssp = ''
        self.prolonged_dssp = ''

        self.all_loops = []
        self.all_regexp = []

    def find_Rosetta_dssp(self):
        print 'About to submit ' + self.pdb_file_name + ' to dssp.xml'

        os.system('/home/labs/fleishman/adig/Rosetta/main/source/bin/rosetta_scripts.static.linuxgccrelease '
                    + ' -s ' + self.path_to_dir + '/' + self.pdb_file_name
                    + ' -parser:protocol ' + self.protocol
                    + ' @' + self.flags
                    + ' >' + self.path_to_dir + '/' + 'dssp.log'
                    + ' -overwrite')

        self.extract_dssp_from_log()

        os.system('rm -rf *_0001.pdb')
        os.system('rm -rf score.sc')

    def extract_dssp_from_log(self):
        dssp_log_file = open(self.path_to_dir + '/dssp.log')
        lines = dssp_log_file.readlines()
        dssp_log_file.close()

        for line in lines:
            if len(line) > 50:
                dssp_line = line.split()
                self.ori_dssp = dssp_line[1]

    def prolong_dssp_to_match_a_seq_with_gaps(self, current_dssp, seq_with_gaps, ignore_gaps_in_end='No'):
        self.prolonged_dssp = current_dssp

        #The following block takes care of the case where there are gaps at the end of the input sequence in the MSA.
	#This is probably a bad solution and I fixed that (see Make_sub...py file). I want to see that it is ok before I delete this block. 28/5/15
        counter = 0
        if ignore_gaps_in_end is 'Ignore':
            for char in reversed(seq_with_gaps):
                if char is '-':
                    counter += 1
                else:
                    break

        for char in range (0, len(seq_with_gaps) - counter):
            if seq_with_gaps[char] == '-':
                self.prolonged_dssp = self.prolonged_dssp[:char] + 'L' + self.prolonged_dssp[char:]

        print self.prolonged_dssp

    def find_loop_motifs_in_dssp(self, dssp_string): #Loop is a sequence of L chars flanked by H or E chars
        flanking_bases = 2
        loop_pos_with_flank = [0, 0]     # Just for initiation
        previous_loop_position = [0, 0]  # Just for initiation

        iter = re.finditer('L{1,}', dssp_string)

        for match in iter:
            previous_loop_position[0] = loop_pos_with_flank[0]
            previous_loop_position[1] = loop_pos_with_flank[1]
            loop_position = match.span()

            if loop_position[0] >= flanking_bases:
                loop_pos_with_flank[0] = loop_position[0] - flanking_bases
            else:
                loop_pos_with_flank[0] = loop_position[0]

            if loop_position[1] <= (len(dssp_string) - flanking_bases):
                loop_pos_with_flank[1] = loop_position[1] + flanking_bases
            else:
                loop_pos_with_flank[1] = loop_position[1]


            #This part takes care of cases where the flanking bases contain L themselves.
            if loop_pos_with_flank[0] - previous_loop_position[1] <= 1: #changed from 2 to 1 on 1.2.15
                if self.all_loops:
                    del self.all_loops[-1]
                    loop_pos_with_flank[0] = previous_loop_position[0]

            self.all_loops.append(loop_pos_with_flank[:]) # [:] allows to append a list object to a list by value

    def match_loops_to_their_regexp_in_the_input_msa(self, input_seq_in_msa):
        hyphen_char = '-'
        letter_char = '[A-Z]'

        for i in range(0, len(self.all_loops)):
            current_loop_regex = ''
            for j in range(self.all_loops[i][0], self.all_loops[i][1]):
                if input_seq_in_msa[j] == hyphen_char:
                    current_loop_regex += hyphen_char
                elif input_seq_in_msa[j].isalpha():
                    current_loop_regex += letter_char
            self.all_regexp.append(current_loop_regex)

    def get_original_dssp(self):
        #By original dssp I mean directly from Rosetta before it is being manipulated
        return self.ori_dssp

    def get_prolonged_dssp(self):
        return self.prolonged_dssp

    def get_all_loops_vec(self):
        return self.all_loops

    def get_all_regexp_vec(self):
        return self.all_regexp
