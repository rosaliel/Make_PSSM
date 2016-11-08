from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import urllib,os,re, subprocess, sys
import linecache
import time
import sys

from path_variables import *
from Get_input_class import *
from Input_class import *
from Blast_class import *
from Handle_msa_class import *
from Thermostab_query import *
from Align2seq_class import *
from Handle_dssp_class import *
from Make_sub_MSA_files import *
from PSSM_class import *
from Handle_PSSMs import *
from Idealize_class import *
from Refinement_class import *
from Filterscan_class import *
from renum_aligned_func import *
from Design_class import *

input = Input()
args = input.get_input()
my_query = input.arrange_input(args)
blast_args = input.blast_args_f(args)

print my_query

if my_query[6]:
    print 'SYMMETRIC DESIGN IS DESIRED. SYMMETRY RULES!!'
    symmetry = 'SYMMETRIC'
else:
    symmetry = 'NON-SYMMETRIC'

pose = Thermostab_query(my_query)
pose.fetch_pdb()
pose.find_interacting_residues()

#Handle symmetry - make symmetric file
if my_query[6]:
    pose.make_symm_file(pose.get_pdb_file_name(), pose.get_path_in())

pose.make_pdb_single_chain()
pose.extract_pdb_full_seq()
pose.extract_seq_from_pdb_file()

pdb_seq = pose.get_pdb_seq()
#input_seq is empty unless this is homology
input_seq = pose.get_input_seq()
if input_seq:
    full_seq = input_seq
    pose.handle_homology_case(pdb_seq, input_seq)
else:
    full_seq = pose.get_full_seq()

if my_query[8] == 'No':
    f = open(pose.get_path_msa() + '/' + 'blast_infile', 'w')
    f.write('>' + my_query[0] + 'Input_pdb_seqres\n')
    f.write(full_seq)
    f.close()
    blast = Blast('blast_infile', pose.get_path_msa() + '/', blast_args[0], blast_args[1], blast_args[2], blast_args[3])
    blast.get_seq_data_as_MSA()
    msa_name = 'query_msa.fa'
    msa_full_path = pose.get_path_msa() + '/' + msa_name

else:
    os.system('cp ' + my_query[8] + ' ' + pose.get_path_msa() + '/')
    msa_name = subprocess.check_output(['ls', '-1', pose.get_path_msa()]).rstrip('\n')
    msa_full_path = pose.get_path_msa() + '/' + msa_name

#The following two lines appear in few locations in the code to enable Jaim to stop jobs
if os.path.isfile(stop_flag):
    sys.exit()

path_in = pose.get_path_in()
pdb_name = pose.get_pdb_single_chain()

pose.renumber_pdb(pdb_name, path_in)
renumbered_pdb_name = pose.get_pdb_renumbered()

#pose.match_interacting_res_to_renum_pdb()

#Uncomment the next call to match... if you wish to reuse renumber.sh (see comment from 17/8/15 on Thermostab_query)
pose.make_coor_cst_file(renumbered_pdb_name, path_in)

msa = Handle_msa(msa_full_path, path_in)
msa.read_fasta_msa()
msa_map = msa.get_name2seq_map()
msa_input_name = msa.get_input_seq_name()
msa_input_seq = msa.get_input_seq()
head_seq_file_path = msa.get_head_seq_path()
subject_filename = 'head_seq'

alignment = Align2seq(pdb_seq, full_seq)
alignment.align2seq()

aligned_pdb_vs_full = alignment.get_aligned_seq1()
aligned_full_seq = alignment.get_aligned_seq2()
print aligned_pdb_vs_full
print ""
print aligned_full_seq


dssp = Handle_dssp(pdb_name, path_in)
dssp.find_Rosetta_dssp()
struct_dssp = dssp.get_original_dssp()
dssp.prolong_dssp_to_match_a_seq_with_gaps(struct_dssp, aligned_pdb_vs_full)
dssp_to_match_fasta = dssp.get_prolonged_dssp()
dssp.prolong_dssp_to_match_a_seq_with_gaps(dssp_to_match_fasta, msa.get_input_seq(), 'Ignore')
dssp_to_match_MSA = dssp.get_prolonged_dssp()
dssp.find_loop_motifs_in_dssp(dssp_to_match_MSA)
dssp.match_loops_to_their_regexp_in_the_input_msa(msa.get_input_seq())

all_loops = dssp.get_all_loops_vec()
all_regexp = dssp.get_all_regexp_vec()

sub_MSA_object = Make_sub_MSA_files(msa_input_name, msa_input_seq, all_loops, all_regexp, msa_map)
sub_MSA_object.create_path_for_output(path_in)
sub_MSA_object.res_renum_ignore_gaps()
sub_MSA_object.iter_over_regexps_to_find_matching_seqs_in_msa()
all_sub_MSA_vec = sub_MSA_object.get_all_sub_MSA_file_names()
sub_MSAs_path = sub_MSA_object.get_sub_MSA_path()

resnum_map = sub_MSA_object.get_res_num_with_hyphens_2_res_num_no_hyphens_map()

pdb_gap_map = renum_aligned_seqs(aligned_pdb_vs_full)
full_seq_gap_map = renum_aligned_seqs(aligned_full_seq)

all_seq_pssm_object = PSSM()
all_seq_pssm_object.create_PSSM_from_MSA(pose.get_path_msa(), msa_name, head_seq_file_path, subject_filename, path_in)

sys.stdout.flush()

for file_name in all_sub_MSA_vec:
    sub_PSSMs_obj = PSSM()
    sub_PSSMs_obj.create_PSSM_from_MSA(sub_MSAs_path, file_name, head_seq_file_path, subject_filename)

pssm_handler = Handle_PSSMs()
pssm_handler.Unite_PSSMs(all_loops, resnum_map, 'query_full_pssm', path_in, sub_PSSMs_obj.get_pssm_output_path())
pssm_handler.remove_non_crystallized_res(aligned_pdb_vs_full, aligned_full_seq)

ex1_t = pose.get_ex1_t()
ex2_t = pose.get_ex2_t()

#############################################Preparations for Refinement##########################################
if os.path.isfile(stop_flag):
    sys.exit()

#define whether needed already done, wheter symmetric and whether need to thread sequence
if my_query[9] == 'No':
    refine = Refinement(pose, symmetry, ex1_t, ex2_t)  #If desired, idealization is also called here.
    refine.refine()
    input_for_fs = refine.best_struct()

else:
    os.system('mkdir ' + pose.get_path_refinement() + '/pdbs')
    os.system('cp ' + my_query[9] + ' ' + pose.get_path_refinement() + '/pdbs')
    input_for_fs = subprocess.check_output(['ls', '-1', pose.get_path_refinement() + '/pdbs']).rstrip('\n')

##############################################################################################################

#####################################Preparations for filterscan##############################################
if os.path.isfile(stop_flag):
    sys.exit()

#define if symmetric and if yes prepare file again after refinement
if my_query[6]:
    pose.make_symm_file(input_for_fs, pose.get_path_refinement() + '/pdbs')
    input_for_fs = pose.make_pdb_single_chain(input_for_fs[:-4] + '_INPUT.pdb', pose.get_path_refinement() + '/pdbs' )

#define keep_n
with open (msa_full_path, "r") as file:
    data=file.read()
    number_of_seqs = data.split('>')
    file.close

    if len(number_of_seqs) < 200:
        keep_native = 0
    else:
        keep_native = 1

#define residues to restrict to repacking
pdb_gaps = re.finditer('-{1,}', aligned_pdb_vs_full)
full_seq_gaps = re.finditer('-{1,}', aligned_full_seq)

restrict_around = 3
restrict_edge_case = 1
match_to_pdb_num = 1 #loop will count from zero but the pdb starts with residue 1
restrict_res = []

all_gaps = [pdb_gaps, full_seq_gaps]

for gaps in all_gaps:

    for match in gaps:
        gap_pos = match.span()

        #Restrict residues before each gap
        if gap_pos[0] > 0: #0 means it is not a gap but a missing N terminus
            gap_pos_copy = gap_pos[0] #can't make assignments for tuple (gap_pos[0])
	    if gap_pos[0] > 3:
                third_before_gap = pdb_gap_map[gap_pos[0]-restrict_around]
	    else:
		third_before_gap = pdb_gap_map[gap_pos[0]-restrict_edge_case]
            #Loop added on 9Feb16 following an edge-case of gap after 2 residues (3oe8)
            while (third_before_gap =="-"):
                gap_pos_copy+=1
                third_before_gap = pdb_gap_map[gap_pos_copy-restrict_around]

            for res_num in range (third_before_gap, third_before_gap + restrict_around):
                restrict_res.append(res_num + match_to_pdb_num)

        #Restrict residues after each gap
        if gap_pos[1]  <= len(aligned_pdb_vs_full) - restrict_around \
                and pdb_gap_map[gap_pos[1]] < pose.get_pdb_length() - restrict_around:
            first_after_gap = pdb_gap_map[gap_pos[1]]
            for res_num in range (first_after_gap, first_after_gap + restrict_around):
                restrict_res.append(res_num + match_to_pdb_num)

        #Restrict gap residues in case the gap is in the threaded input and therefore represents residues in the pdb.
        if gaps == full_seq_gaps:
            gap_len = gap_pos[1] - gap_pos[0]
            for res_num in range (third_before_gap + restrict_around, first_after_gap):
                restrict_res.append(res_num + match_to_pdb_num)

res_to_restrict_str = '-1'
if restrict_res:
    res_to_restrict_str = str(restrict_res[0])
    for i in range (1, len(restrict_res)):
        res_to_restrict_str += ',' + str(restrict_res[i])

    print '\nresidues to restrict to repack:'
    print restrict_res
    print res_to_restrict_str

#run filterscan
filterscan = Filterscan(pose, symmetry, input_for_fs, keep_native, res_to_restrict_str, ex1_t, ex2_t)
filterscan.filterscan()
filterscan.verify_output()

######################################################################################
sys.stdout.flush()

if os.path.isfile(stop_flag):
    sys.exit()

design = Design(pose, symmetry, input_for_fs, keep_native, filterscan.get_resfiles(), ex1_t, ex2_t)
design.design()

if os.path.isfile(stop_flag):
    sys.exit()

print '\nThe algorithm has reached the end of the process'

print '\nNow analyzing results and sequence diversity'

sys.stdout.flush()
time.sleep(4)

# Analysis of design results
path_analyze = pose.get_path_analysis()
path_fs = pose.get_path_fs()

if my_query[8] == 'No':
    os.system('sh ' + gen_path + 'analyze_seq_info.sh '
              + msa_full_path + ' '
              + path_in + ' '
              + '>' + path_analyze + '/' + 'MSA_analysis')

#time.sleep(2)

#threshold = "-0.45"
#os.system('rm -rf ' + path_fs + '/scores/score.sc')
#os.system('sh ' + gen_path + 'filtered_score_files.sh '
#          + threshold + ' '
#          + path_fs + '/scores/ '
#          + path_analyze + '/' )

#time.sleep(60)

#os.system('sh ' + gen_path + 'analyze_fs_results.sh '
#          + threshold + ' '
#          + path_analyze + '/ '
#          + path_fs + '/ '
#          + input_for_fs)

