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

# Parsing input
input = Input()
args = input.get_input()
blast_args = input.blast_args_f(args)

# Download PDB & fasta, clean and extract sequences
pose = Thermostab_query(args)
pose.fetch_pdb()
pose.make_pdb_single_chain()
pose.extract_pdb_full_seq()
pose.extract_seq_from_pdb_file()
pdb_seq = pose.get_pdb_seq() # seq of pdb
full_seq = pose.get_full_seq() # seq of fasta

# Creating an MSA
if args['input_msa'] == 'No':
    f = open(pose.get_path_msa() + '/' + 'blast_infile', 'w')
    f.write('>' + args['pdbid'] + 'Input_pdb_seqres\n')
    f.write(full_seq)
    f.close()
    blast = Blast('blast_infile', pose.get_path_msa() + '/', blast_args[0], blast_args[1], blast_args[2], blast_args[3])
    blast.get_seq_data_as_MSA()
    msa_name = 'query_msa.fa'
    msa_full_path = pose.get_path_msa() + '/' + msa_name
else:
    os.system('cp ' + args['input_msa'] + ' ' + pose.get_path_msa() + '/')
    msa_name = subprocess.check_output(['ls', '-1', pose.get_path_msa()]).rstrip('\n')
    msa_full_path = pose.get_path_msa() + '/' + msa_name

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
