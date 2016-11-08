from Bio.SubsMat import MatrixInfo as matlist
from prody import *
from pylab import *
import urllib,os,re, subprocess, sys
from os.path import basename
from Align2seq_class import *
from Blast_class import *
from path_variables import *

ion()

class Thermostab_query:

    def __init__(self, args):

        self.args = args
        if len(args['pdbid']) == 4:
            self.pdbid = args['pdbid']
            self.uploaded_full_path = None
        else:
            tmp = basename(args['pdbid']).strip()
            self.pdbid = tmp[:-4]
            self.uploaded_full_path = args['pdbid']

        self.input_seq = args['input_seq'].replace(" ", "").replace("'\n'","")
        self.chain_id = args['chain_id']
        self.small_molecules = args['small_ligands']
        self.interacting_chains = args['large_ligands']
        self.res_to_fix_str = args['res_to_fix']
        self.symm_type = args['symm_interfaces'][0:2]
        self.oligomers = args['symm_interfaces'][3:]
        self.idealize = args['idealize']
        self.msa = args['input_msa']
        self.need_refine = args['refined_file']
        self.uploaded_fasta = args['uploaded_fasta']
        self.cst_value = args['cst_value']
        self.homology_id = 100 #default for queries with structure
        self.ex1_t = 700 #for proteins longer than ex1_t, the ex1 flags won't be used
        self.ex2_t = 500 #for proteins longer than ex2_t, the ex2 flags won't be used

    def __str__(self):
        report_str = '\nYour input structure is ' + self.pdbid
        if self.input_seq != None:
            report_str += ' and your input sequence is: \n' + str(self.input_seq)
        report_str += '\nChain to design is: ' + self.chain_id + '\n'
        if self.small_molecules != None:
            report_str += 'The following small ligands are interacting with the protein: '
            for idx in range (0, len(self.small_molecules)):
                report_str += str(self.small_molecules[idx]) + ' '
            report_str += '\n'
        if self.interacting_chains != None:
            report_str += 'The following chains are interacting with chain ' + self.chain_id + ': '
            for idx in range (0, len(self.interacting_chains)):
                report_str += str(self.interacting_chains[idx]) + ' '
            report_str += '\n'
        return report_str

    def fetch_pdb(self):
        # make directory and sub directories.
        self.create_dirs()

        # Structure file is uploaded by user
        if self.uploaded_full_path:
            os.system('cp ' + self.uploaded_full_path + ' ' + self.get_path_in() + '/')
            self.pdbid = subprocess.check_output(['ls', '-1', self.get_path_in()]).rstrip('\n')[:-4]
            self.pdb_file_name = str(self.pdbid) + '.pdb'
        # Structure is downloaded from the protein data bank
        else:
            url_pdb = 'http://www.rcsb.org/pdb/files/%s.pdb' % self.pdbid
            self.pdb_file_name = str(self.pdbid) + '.pdb'
            urllib.urlretrieve(url_pdb, self.path_in + '/' + self.pdb_file_name)

        self.get_chain_identifiers()

    def _download_fasta(self, pdbid, path):
        """Returns the fasta file name"""
        url_fasta = 'http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=' + pdbid
        fasta_file_name = pdbid + ".fasta"
        urllib.urlretrieve(url_fasta, path + '/' + fasta_file_name)
        return fasta_file_name

    def get_chain_identifiers(self):
        """The method stores the single letter code of all pdb chains in a list (method called by fetch_pdb())
        and downloads the fasta (or copy it)
        """
        if self.uploaded_full_path:
            if self.uploaded_fasta:
                self.fasta_file_name = 'uploaded_fasta.fa'
                os.system('cp ' + self.uploaded_fasta + ' ' + self.get_path_in() + '/' + self.fasta_file_name)
                os.system('dos2unix -c Mac ' + self.path_in + '/' + self.fasta_file_name)
            else: # uploaded pdb but not fasta
                self.fasta_file_name = self._download_fasta(self.pdbid[:4], self.path_in)
        else:
            print 'Downloading fasta file from web'
            self.fasta_file_name = self._download_fasta(self.pdbid, self.path_in)

        self.pdb_chain_list = []
        with open(self.path_in + '/' + self.fasta_file_name) as file:
            for line in file:
                if '>' in line:
                    for i in range (0, len(line)):
                        if (line[i] == ':'):
                            self.pdb_chain_list.append(line[i+1])
                            break
            file.close()

        self.chain_id_to_num_map = {}
        for chain_num in range (0, len(self.pdb_chain_list)):
            self.chain_id_to_num_map[self.pdb_chain_list[chain_num]] = chain_num

    def extract_seq_from_pdb_file(self):
        #The following command calls the pdb2seq.sh script (slightly modified) and stores the output in a variable.
        pdb2seq_output = subprocess.check_output([gen_path + 'pdb2seq.sh', self.path_in + '/' + self.pdb_single_chain])
        self.chain_seq = pdb2seq_output.rstrip()

    def extract_pdb_full_seq(self):
        #Extracting the full_sequence of the pdb structure from the rcsb display Fasta sequence file.
        #Make sure the pdb file does not contain BXXX residues (for example: BGLN)
        fasta=open(self.path_in + '/' + self.fasta_file_name,'r')
        list_of_lines = fasta.readlines()

        lines_counter = 0     #counts all the lines in list_of_lines
        arrows_counter = 0    #counts all lines with > sign in list_of_lines
        arrows_in_lines = []  #stores the number of lines with > sign

        for line in list_of_lines:
            if '>' in line:
                arrows_in_lines.insert(arrows_counter, lines_counter)
                arrows_counter += 1

            lines_counter += 1

        first_seq_line = (arrows_in_lines[self.chain_id_to_num_map[self.chain_id]]) + 1
        if self.chain_id !=  self.pdb_chain_list[-1]:
            last_seq_line = arrows_in_lines[self.chain_id_to_num_map[self.chain_id] + 1]
        else: #special case for last chain since there is no > after the last seq line.
            last_seq_line = lines_counter

        self.full_seq = ''
        for i in range (first_seq_line, last_seq_line):
            self.full_seq += list_of_lines[i].strip()
        
        self.full_seq = self.full_seq.upper()
        print 'NOW PRINTING FASTA SEQUENCE FILE UPPERCASE'
        print self.full_seq

        return self.full_seq

    def create_dirs(self):
        """The method creates a new path to store future files (method called by fetch_pdb())."""
        self.path = query_path + self.pdbid + '_'
        while os.path.exists(self.path):
            self.path += '_'

        self.path_in = self.path + '/in'
        self.path_msa = self.path + '/msa'
        sub_dirs = [self.path_in, self.path_msa]
        os.makedirs(self.path)
        for dir in sub_dirs:
            os.makedirs(os.path.join(self.path, dir))

        # Print the args from flag file used in the command line into a file.
        f = open(self.path + '/used_flags', 'w')
        for arg in self.args.keys():
            if isinstance(self.args[arg], list):
                line = ''
                for item in self.args[arg]:
                    if line == '':
                        line = item
                    else:
                        line += ',' + item
                    # f.write(arg_vec[counter] + line + '\n')
                f.write('-{}={}\n'.format(arg, line))
            else:
                f.write('-{}={}\n'.format(arg, str(self.args[arg])))
        f.close()

    def find_interacting_residues(self):
        #The method finds all the residues in the protein of interest that are interacting with ligands.
        #The functions here are taken fro ProDy (http://prody.csb.pitt.edu/index.html)
        self.parsed_pdb = parsePDB(self.path_in + '/' + self.pdbid + '.pdb')

        #The loop finds residues interacting with small ligand that their name was given as input by the user.
        self.interacting_res = []
        for molecule in self.small_molecules:
            if molecule in ('ZN', 'CA', 'MG', 'SE', 'CU', 'FE', 'CO'):
                cutoff = 5
            else:
                cutoff = 8
            print molecule + ': cutoff = ' + str(cutoff)

            try:
                select_res = self.parsed_pdb.select('protein chain ' + self.chain_id + ' within ' + str(cutoff) + ' of resname ' + molecule)
                self.interacting_res.append(sorted(set(select_res.getResnums())))

            except AttributeError:
                print "The small ligand " + molecule + " is not in contact with the chain of interest"

        #The loop finds residues interacting with other (dna/protein) chains in the structure if their id was given by the user.
        for chain in self.interacting_chains:
            try:
                select_res = self.parsed_pdb.select ('protein chain ' + self.chain_id + ' within 5 of chain ' + chain)
                self.interacting_res.append(sorted(set(select_res.getResnums())))
            except AttributeError:
                print "The large ligand " + chain + " is not in contact with the chain of interest"

        #The next block adds specific residues that the user wants to fix to the - self.interacting_res vector.
        if self.res_to_fix_str:
            fix_res_from_user_str = sorted(set(self.res_to_fix_str.split(',')))
            if len(fix_res_from_user_str) != 0:
                fix_res_from_user_int = []
                for res in fix_res_from_user_str:      #The loop converts the strings to integers.
                    fix_res_from_user_int.append(int(res))

                self.interacting_res.append(fix_res_from_user_int)

        print self.interacting_res
        print '' #can't put '\n' sign in previous line due to concatenation

    def make_symm_file(self, input_pdb, in_out_path):
        self.symm_filename = str(self.pdbid) + '_' + self.symm_type + '.symm'
        print 'SYMMETRIC FILE NAME IS:' + self.symm_filename
        self.symm_full_path = in_out_path + '/' + self.symm_filename

        cmd_line = 'perl ' + gen_path + 'make_symmdef_file.pl -m NCS' \
                   + ' -p ' + in_out_path + '/' + input_pdb \
                   + ' -a ' + self.chain_id \
                   + ' -i ' + self.oligomers \
                   + ' >' + in_out_path + '/' + self.symm_filename

        os.system(cmd_line)

        self.pdb_file_name = input_pdb[:-4] + '_INPUT.pdb'

    def make_pdb_single_chain(self, input_pdb = None, in_out_path = None):
        #This method prepares the downloaded pdb file for rosettascripts.
        #Gets rid of all irrelevant lines, chains, multiple rotamers and non amino acid HETATMs
        if input_pdb:
            print 'In make_pdb_single_chain with parameters'
            self.pdb_clean = input_pdb[:-4] + 'clean.pdb'
            self.pdb_single_chain = input_pdb[:-4] + '_sc.pdb'
            pdb_file_single_chain = open(in_out_path + '/' + self.pdb_single_chain, 'a')
            file_to_read = in_out_path + '/' + input_pdb
        else:
            print 'In make_pdb_single_chain DEFAULT VERSION '
            self.pdb_clean = self.pdbid + '_' + 'clean.pdb'
            self.pdb_single_chain = self.pdbid + '_' + self.chain_id + '.pdb'
            pdb_file_single_chain = open(self.path_in + '/' + self.pdb_single_chain, 'a')
            file_to_read = self.path_in + '/' + self.pdb_file_name

        os.system('sh ' + gen_path + 'clean_pdb.sh '
                  + file_to_read + ' '
                  + self.path_in + '/' + self.pdb_clean)

        for line in open (self.path_in + '/' + self.pdb_clean):
            if line[21:22] == self.chain_id and line [16:17] != 'B':
                pdb_file_single_chain.write(line)

        pdb_file_single_chain.close()

        return self.pdb_single_chain

    def renumber_pdb(self, pdb_name, pdb_path):
        #########################################Artificial!!!############################################
        # 17/8/15                                                                                        #
        # I wished to get rid of renumber_pdb with minimum changes in the code. So I call it but instead #
        # of renumbering I just copy the input file and give it the renumbered name.                     #
        # Note that if renumber_pdb is not used, avoid matching functions - no need to match anymore     #
        ##################################################################################################

        #This method works only with Adi's version of renumber.sh in which the output is printed out.
        self.renum_pdb_file = self.pdbid + '_' + self.chain_id + '_renum.pdb'

        #Relevant for homolougy cases where a renum file already exists from a previous step.
        os.system('rm -rf ' + self.path_in + '/' + self.renum_pdb_file)

        #Comment the next cp command if you wish to reuse renumber.sh (see comment from 17/8/15 above)
        os.system('cp ' + pdb_path + '/' + pdb_name + ' ' + self.path_in + '/' + self.renum_pdb_file)

        #Uncomment the next 2 lines if you wish to reuse renumber.sh (see comment from 17/8/15 above)
        #renumbered_file = open(self.path_in + '/' + self.renum_pdb_file, 'a')
        #subprocess.call([gen_path + 'renumber_auto.sh', pdb_path + '/' + pdb_name], stdout = renumbered_file)

    def make_coor_cst_file(self, pdb_name, pdb_path):
        self.cst_full_path = self.path_in + '/' + self.pdbid + '_coorcst'
        cst_file = open(self.cst_full_path, 'a')
        subprocess.call([gen_path + 'make_csts.sh', pdb_path + '/' + pdb_name], stdout = cst_file)

        return self.cst_full_path

    def match_interacting_res_to_renum_pdb(self):
        #This method gets the vector of interacting residues (that were defined on the original pdb)
        #and finds the matching residue number in the renumbered pdb file.

        all_interacting_res = [] #store all interacting residues as a simple list (rather than list of lists)
        ori_num_vec = []         #store all residue numbers from the original pdb chain - len(this_vec) = len(renum_vec)
        renum_vec = []           #store all residue numbers from the renumbered pdb file - len(this_vec) = len(ori_num_vec)

        for interacting_group in self.interacting_res:
            for res in interacting_group:
                all_interacting_res.append(res)

        all_interacting_res = sorted(set(all_interacting_res))


        for line in open(self.path_in + '/' + self.pdb_single_chain):
            columns = line.split()
            if columns[2] == 'CA' and columns[4] == self.chain_id:
                ori_num_vec.append(int(columns[5]))

        for line in open(self.path_in + '/' + self.renum_pdb_file):
            columns = line.split()
            if columns[2] == 'CA' and columns[4] == self.chain_id:
                renum_vec.append(int(columns[5]))

        if len(ori_num_vec) != len(renum_vec):
            print 'WARNING: original file and renumbered file have different number of CA atoms'
            #This can indicate an error or a homolougy case where residue were deleted from the pdb file.
        counter = 0

        self.interacting_res_renumbered = []
        for res in all_interacting_res:
            while counter != len(ori_num_vec):
                if res == ori_num_vec[counter]:
                    self.interacting_res_renumbered.append(renum_vec[counter])
                    break
                counter += 1

        #The next line is meant to make sure that interacting_res_renumbered is not empty -
        #the refinement protocol can't get an empty vector and this is the easiest solution.
        self.interacting_res_renumbered.append(1)

        print 'length of original residue numbers vector followed by that of the renumbered vec: ' + str(len(ori_num_vec)) + ' ' + str(len(renum_vec))
        print all_interacting_res
        print self.interacting_res_renumbered
        print len(all_interacting_res), len(self.interacting_res_renumbered)

    def match_interacting_res_to_threaded_renum_pdb(self, gaps_vec):
        #If the threaded sequence is shorter than the pdb in some places, residues will be deleted from the pdb file
        #using the function self.remove_res_from_pdb(). In such a case numbers are changing and we must re-match the
        #interacting residue files to the new pdb file.
        #This is done differently than in self.match_interacting_res_to_renum_pdb()

        inter_res_to_remove = []
        for res_idx in range (0, len(self.interacting_res_renumbered)):
            counter = 0
            for gap in gaps_vec:
                if self.interacting_res_renumbered[res_idx] < gap:
                    self.interacting_res_renumbered[res_idx] = self.interacting_res_renumbered[res_idx] - counter
                    break
                elif self.interacting_res_renumbered[res_idx] == gap:
                    inter_res_to_remove.append(res_idx)
                else:
                    counter += 1
                    if counter == len(gaps_vec):
                        self.interacting_res_renumbered[res_idx] = self.interacting_res_renumbered[res_idx] - counter

        reversed(sorted(inter_res_to_remove))
        for idx in inter_res_to_remove:
            self.interacting_res_renumbered.pop(idx)

        print self.interacting_res_renumbered

    def handle_homology_case(self, pdb_seq, input_seq):
        #For homolog case this function does the follow:
        # 1. Removes extra residues that may exist in the input_seq (threaded) but not in the structure
        # 2. Adds to the threaded seq residues that are present in the pdb but not in the threaded.
        # 3. Finds the identity between input && pdb to determine the coorcst value later

        #Align
        thread_align = Align2seq(pdb_seq, input_seq)
        thread_align.align2seq()

        aligned_pdb = thread_align.get_aligned_seq1()
        aligned_thread = thread_align.get_aligned_seq2()
        print "aligned_pdb is "
        print aligned_pdb
        print "alinged_thread is "
        print aligned_thread

        # 1. Finds all gaps in the aligned_pdb variable, store their indices in a vector
        # 2. Removes these gaps from the threaded sequence since they do not exist in the pdb
        # 3. Removes these gaps from the pdb sequence - important for the next 'if' section
        gaps_in_struct = thread_align.find_gaps(aligned_pdb)
        if len(gaps_in_struct) > 0:
            thread_align.remove_gaps(aligned_pdb, gaps_in_struct)
            aligned_thread = thread_align.remove_gaps(aligned_thread, gaps_in_struct)

        # 1. Finds all gaps in the threaded sequence and store their indices in a vector
        # 2. Replaces these gaps (-) with the matching amino acid in the aligned_pdb sequence.
        gaps_in_input_seq = thread_align.find_gaps(aligned_thread)
        if len(gaps_in_input_seq) > 0:
            for gap_idx in gaps_in_input_seq:
                aligned_thread = aligned_thread[:gap_idx] + aligned_pdb[gap_idx:gap_idx + 1] + aligned_thread[gap_idx + 1:]

        self.input_seq = aligned_thread
        os.system('echo ' + self.input_seq + ' >' + self.path_in + '/threaded_seq')

        cmp_f = open(self.path_in + '/bl2seq2compare',   'w')
        cmp_f.write('>pdb_seq' + '\n' + self.full_seq)
        cmp_f.close()

        in_f = open(self.path_in + '/bl2seq_input', 'w')
        in_f.write('>input_seq' + '\n' + self.input_seq)
        in_f.close()

        blast = Blast('bl2seq_input', self.path_in + '/', 'default', 'default', 'default', 'default', self.path_in + '/bl2seq2compare')
        parsed_bl2seq = blast.parse_xml_hit(blast.parse_blast_output(blast.blast2seq_p()))
        self.homology_id = parsed_bl2seq[0][5]

    def remove_res_from_pdb(self, res_list, chain_id, pdb_file_name = None):
        #Function is used in the case of homoloug structure in case there are residue that are present in the
        #structure but not in the input sequence. Removes the extra residues from the pdb file to match the input.
        if pdb_file_name == None:
            pdb_file_name = self.pdbid + '_' + chain_id + '_renum'

        str_for_bash = 'cat ' + self.path_in + '/' + pdb_file_name + '.pdb'

        #The next loop takes care of the different number of spaces between col5 and col6 of a pdb file
        for res in res_list:
            if res < 10:
                spaces = '   '
            elif res < 100:
                spaces = '  '
            elif res < 1000:
                spaces = ' '

            str_for_bash += ' | grep -v \" ' + chain_id + spaces + str(res) + '\"'

        threaded_pdb_file_name = pdb_file_name + '_thread.pdb'

        os.system(str_for_bash + '>' + self.path_in + '/' + threaded_pdb_file_name)

        return threaded_pdb_file_name

    def get_ex1_t(self):
        return self.ex1_t

    def get_ex2_t(self):
        return self.ex2_t

    def get_pdbid(self):
        return self.pdbid

    def get_pdb_file_name(self):
        return self.pdb_file_name

    def get_path_in(self):
        return self.path_in

    def get_pdb_single_chain(self):
        return self.pdb_single_chain

    def get_path_msa(self):
        return self.path_msa

    def get_pdb_renumbered(self):
        return self.renum_pdb_file

    def get_cst_file(self):
        return self.cst_full_path

    def get_coorcst_value(self):
        return self.cst_value

    def get_res_to_fix(self):

        all_interacting_res = [] #store all interacting residues as a simple list (rather than list of lists)

        for interacting_group in self.interacting_res:
            for res in interacting_group:
                all_interacting_res.append(res)

        res_to_fix_string = '-1'
        if all_interacting_res:
            res_to_fix_string = str(all_interacting_res[0]) + self.chain_id #special case: no comma before.
            for i in range (1, len(all_interacting_res)): #all other residues.
                res_to_fix_string += ',' + str(all_interacting_res[i]) + self.chain_id

        return res_to_fix_string

    def get_symm_full_path(self):
        return self.symm_full_path

    def get_chain_id(self):
        return self.chain_id

    def get_pdb_seq(self):
        return self.chain_seq

    def get_pdb_length(self):
        return len(self.chain_seq)

    def get_full_seq(self):
        return self.full_seq

    def get_input_seq(self):
        #Get the sequence to thread on the structure in case it's a homoloug case.
        return self.input_seq

    def idealize_YorN(self):
        return self.idealize
