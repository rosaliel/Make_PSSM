from Bio.SubsMat import MatrixInfo as matlist
from prody import *
from pylab import *
import urllib,os,re, subprocess, sys

ion()

from Thermostab_query import *
from sleep import *
from path_variables import *

class Filterscan:

    def __init__(self, pose, symm, pdb_input, keep_native, restrict_res, ex1_t, ex2_t):
        self.pose = pose
        self.symm = symm
        self.pdb_input = pdb_input
        self.pdb_full_path = self.pose.get_path_refinement() + '/pdbs/' + self.pdb_input	
        self.restrict_res = restrict_res
        self.ex1_t = ex1_t #t for threshold
        self.ex2_t = ex2_t
        self.t_vec = [-2, -1.8, -1.5, -1.25, -1, -0.75, -0.45]
	self.delay_t = 400

        #Define which xml and flags file to use
        if not self.symm == 'NONE':
            self.script_path = protocols_path + 'filterscan_symm_auto.xml'
        else:
            self.script_path = protocols_path + 'filterscan_auto.xml'

        self.flags_path = flags_path + 'flags_classic'

        self.path_fs = self.pose.get_path_fs()
        self.pdb_len = self.pose.get_pdb_length()
        self.pssm_path = self.pose.get_path_in() + '/' + 'final_pssm_for_Rosetta'
        self.fix_res = self.pose.get_res_to_fix()
        self.keep_n = str(keep_native)

        #I make a new cst_file for homology since I want to avoid bias to the original structure.
        if self.pose.get_input_seq():
            self.cst_path = pose.make_coor_cst_file(self.pdb_input, self.pose.get_path_refinement() + '/pdbs')
        else:
            self.cst_path = self.pose.get_cst_file()

    def __str__(self):
        pass

    def filterscan(self):

        self.scores_path = self.path_fs  + '/scores/'
        self.resfiles_path = self.path_fs + '/resfiles/'
        self.pdbs_path_and_name = self.path_fs + '/pdbs/' + self.pdb_input
        job_command = 'sh ' + gen_path + 'create_job_auto.sh ' \
                        + self.path_fs \
                        + ' -s ' + self.pose.get_path_refinement() + '/pdbs/' + self.pdb_input \
                        + ' -parser:protocol ' + self.script_path \
                        + ' -parser:script_vars keep_n=' + self.keep_n \
                        + ' -parser:script_vars res_to_fix=' + self.fix_res \
                        + ' -parser:script_vars pdb_reference=' + self.pdb_full_path \
                        + ' -parser:script_vars res_to_restrict=' + self.restrict_res \
                        + ' -parser:script_vars cst_full_path=' + self.cst_path \
                        + ' -parser:script_vars cst_value=' + str(self.pose.get_coorcst_value()) \
                        + ' -parser:script_vars pssm_full_path=' + self.pssm_path \
                        + ' -parser:script_vars scores_path=' + self.scores_path \
                        + ' -parser:script_vars pdb_path_and_name=' + self.pdbs_path_and_name \
                        + ' -parser:script_vars resfiles_path=' + self.resfiles_path \
                        + ' @' + self.flags_path \
                        + ' -overwrite '

        if self.pdb_len < self.ex1_t:
            print '\nAdding The ex1 Flag to Rosetta Simultations'
            add_to_command = '-ex1 '
            job_command = job_command + add_to_command

        if self.pdb_len < self.ex2_t:
            print '\nAdding The ex2 Flag to Rosetta Simultations'
            add_to_command = '-ex2 '
            job_command = job_command + add_to_command

        if self.pdb_len > self.delay_t:
            delay = int(self.pdb_len * 1.2)
            print '\nAdding Delay Flags: -random_delay ' + str(delay)
            add_to_command = '-random_delay ' + str(delay) + ' '
            job_command = job_command + add_to_command

        else:
            print '\nAdding Delay Flag: -random_delay 400'
            add_to_command = '-random_delay 400 '
            job_command = job_command + add_to_command 

        if not self.symm == 'NONE':
            self.t_vec = [x*(int(self.symm[1])) for x in self.t_vec]  #multiply the fs thresholds by the number of monomers (if dihedral represents one layer)

            if self.symm[0] == 'D':
                self.t_vec = [x*2 for x in self.t_vec] #multiply fs thresholds by two if symmetry is dihedral           
            
            print self.t_vec

            add_to_command = ' -parser:script_vars symmetry_file=' + self.pose.get_symm_full_path() \
                           + ' -parser:script_vars symm_thresholds=1,' #a positive threshold to make sure that wt PDBs are dumped. 

            for threshold in self.t_vec: 
                if threshold == self.t_vec[-1]:
                    add_to_command += str(threshold)
                else:
                    add_to_command += str(threshold) + ','

            job_command = job_command + add_to_command

        # Prepare the final command for each residue in the protein.
        for j in range(0, self.pdb_len):
            job_command_per_res = job_command + ' -parser:script_vars current_res=' + str(j + 1)
            os.system(job_command_per_res)

        print '\nAbout to submit ' + self.pdb_input + ' to filterscan'

        #run the command file
        os.system('chmod +x ' + self.path_fs + '/command')
        os.system(self.path_fs + '/command')

    def verify_output(self):

        self.out_files = []
        jobs = subprocess.Popen(['ls', '-1', self.path_fs + '/job/'], stdout=subprocess.PIPE).stdout.read().split('\n')
        for job in jobs:
            if 'job.' in job:
                self.out_files.append('out.' + job[4:])

        self.t_files = []
        for threshold in self.t_vec:
            if (threshold - int(threshold)) == 0:
                threshold = int(threshold)
                print 'Rounding threshold to ' + str(threshold)

            name = 'designable_aa_resfile.' + str(threshold)
            self.t_files.append(name)

        repeats = 12000
        sleep_and_check('Filterscan', repeats, self.out_files, self.path_fs + '/out/')
        for out_file in self.out_files:
            line = subprocess.Popen(['grep', 'reported success', self.path_fs + '/out/' + out_file], stdout=subprocess.PIPE).stdout.read()
            if not 'success' in line:
                print 'job.' + out_file[4:] + ' FAILED! check err file'

        repeats = 100
        
        sleep_and_check('Filterscan_resfiles', repeats, self.t_files, self.resfiles_path)


    def get_resfiles(self):
        return self.t_files


# What to add to this class:
# 1. Get filtered scored files
# 2. detect artifact
