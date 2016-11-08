#pdbid = ''

#running_notes = open('running_notes', 'a')
#running_notes.write('Running the sub_pssm code for ' + pdbid + '\n\n')

from Thermostab_query import *

class Handle_msa:

    def __init__(self, msa_path, path_for_output):
        self.msa_path = msa_path
        self.name2seq_map = {}
        self.input_seq_name = ''
        self.path_for_output = path_for_output

    def read_fasta_msa(self):
        msa = open(self.msa_path, 'r')
        line_num = 1 #first line in text file is 1 in python
        line_num_vec = []

        homolog_name = ''
        homolog_seq = ''
        for line in msa:
            line = line.strip()

            if '>' in line:
                if 'THIS_IS_QUERY' in line:
                    self.input_seq_name = line
                if len(homolog_name) > 1:
                    self.name2seq_map[homolog_name] = homolog_seq

                homolog_name = line
                homolog_seq = ''
            else:
                homolog_seq += line

        self.name2seq_map[homolog_name] = homolog_seq

        self.make_head_seq_file(self.get_input_seq_name(), self.get_input_seq())

    def read_clustal_msa(self):
        msa = open(self.msa_path, 'r')
        line_number = 0

        for line in msa:
            line = line.strip()
            if (line != '') and (line != 'CLUSTAL') and not ('/apps/RH6U4/muscle/3.8.425/bin/muscle' in line):
                name_seq_pair = line.split()
                line_number += 1
                if 'THIS_IS_QUERY' in name_seq_pair[0]:
                    self.input_seq_name = name_seq_pair[0]

                #Later add code here with this if statement that throws an error if a name_seq_pair is not really a pairs.
                if len(name_seq_pair) == 2:
                    key = name_seq_pair[0]
                    if key in self.name2seq_map:
                        self.name2seq_map[name_seq_pair[0]] += name_seq_pair[1]
                    else:
                        self.name2seq_map[name_seq_pair[0]] = name_seq_pair[1]
                else:
                    print 'WARNING (Handle_msa_class): LINE NUMBER ' + str(line_number) + ' IS NOT IN CORRECT FORMAT!!!'
                    print 'THE CORRECT FORMAT IS: SEQ_NAME + SPACE(S) + SEQ'

        #running_notes.write('{:<40}'.format('Input sequence length in the MSA is: ') \
        #                   + '{:<40}'.format(str(len(self.name2seq_map[self.input_seq_name]))) + '\n\n')

        self.make_head_seq_file(self.get_input_seq_name(), self.get_input_seq())


    def get_name2seq_map(self):
        return self.name2seq_map


    def get_input_seq_name(self):
        if len(self.input_seq_name) > 0:
            return self.input_seq_name
        else:
            print 'ERROR: NO SEQUENCE IN THE MSA WAS FOUND TO BE THE INPUT SEQUENCE'
            print 'THE CORRECT FORMAT IS: Input_pdb_SEQRES_<CHAIN>/<1-FINAL_RES>'


    def get_input_seq(self):
        key = self.get_input_seq_name()
        if key in self.name2seq_map:
            return self.name2seq_map[key]
        else:
            print 'ERROR: ' + key + ' IS NOT IN MAP!'


    def make_head_seq_file(self, input_seq_name, msa_input_seq):
        head_seq_file = open(self.path_for_output + '/' + 'head_seq', 'a')
        head_seq_file.write('>' + input_seq_name + '\n')
        head_seq_file.write(msa_input_seq)
        head_seq_file.close()


    def get_head_seq_path(self):
        return self.path_for_output
