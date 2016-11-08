from Bio.SubsMat import MatrixInfo as matlist
from prody import *
from pylab import *
import urllib,os,re, subprocess, sys

ion()

from Thermostab_query import *
from sleep import *
from path_variables import *

class Filterscan:

    def __init__(self, pose, symm, pdb_input, keep_native, restrict_res):
        self.pose = pose
        self.symm = symm
        self.pdb_input = pdb_input
        self.pdb_full_path = self.pose.get_path_refinement() + '/pdbs/' + self.pdb_input	
        self.restrict_res = restrict_res

        #Define which xml and flags file to use
        if self.symm == 'SYMMETRIC':
            self.script_path = protocols_path + 'filterscan_symm_auto.xml'
        else:
            self.script_path = protocols_path + 'filterscan_auto.xml'

        self.flags_path = flags_path + 'flags_classic_fs'

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

        if self.symm == 'SYMMETRIC':
            add_to_command = ' -parser:script_vars symmetry_file=' + self.pose.get_symm_full_path()
            job_command = job_command + add_to_command


        print 'About to submit ' + self.pdb_input + ' to filterscan'

        # Prepare the final command for each residue in the protein.
        for j in range(0, self.pdb_len):
            job_command_per_res = job_command + ' -parser:script_vars current_res=' + str(j + 1)
            os.system(job_command_per_res)

        #run the command file
        os.system('chmod +x ' + self.path_fs + '/command')
        os.system(self.path_fs + '/command')

    def verify_output(self):

        self.err_files = []
        self.out_files = []
        jobs = subprocess.Popen(['ls', '-1', self.path_fs + '/job/'], stdout=subprocess.PIPE).stdout.read().split('\n')
        for job in jobs:
            if 'job.' in job:
                self.err_files.append('err.' + job[4:])
                self.out_files.append('out.' + job[4:])

        if self.symm == 'SYMMETRIC':
            t_vec = [-4, -3, -2.5, -2, -1.5, -0.9, -0.1]
        else:
            t_vec = [-2, -1.8, -1.5, -1.25, -1, -0.75, -0.45]

        self.t_files = []
        for threshold in t_vec:
            name = 'designable_aa_resfile.' + str(threshold)
            self.t_files.append(name)

        #Err files are printed at the end of the process so if they exist the process ended.
        repeats = 30000
        sleep_and_check('Filterscan_err', repeats, self.err_files, self.path_fs + '/err/')
        for out_file in self.out_files:
            line = subprocess.Popen(['grep', 'reported success', self.path_fs + '/out/' + out_file], stdout=subprocess.PIPE).stdout.read()
            if not 'success' in line:
                print 'job.' + out_file[4:] + ' FAILED! check err file'

        repeats = 100
        sleep_and_check('Filterscan_resfiles', repeats, self.t_files, self.resfiles_path)


    def get_resfiles(self):
        return self.t_files


# What to add to this class:
# 1. Verification that all jobs did not fail (look for success in the out or use the err directory)
# 2. Get filtered scored files
# 3. detect artifact
