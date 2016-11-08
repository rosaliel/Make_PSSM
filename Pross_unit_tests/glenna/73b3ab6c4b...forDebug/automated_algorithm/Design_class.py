import urllib,os,re, subprocess, sys
from Thermostab_query import *
from sleep import *

from Thermostab_query import *
from path_variables import *

class Design:

    def __init__(self, pose, symm, pdb_input, keep_native, in_resfiles):
        self.pose = pose
        self.symm = symm
        self.pdb_input = pdb_input
        self.pdb_full_path = self.pose.get_path_refinement() + '/pdbs/' + self.pdb_input
        self.keep_n = str(keep_native)
        self.in_resfiles = in_resfiles
        self.path_design = self.pose.get_path_design()
        self.pssm_path = self.pose.get_path_in() + '/' + 'final_pssm_for_Rosetta'
        self.cst_path = self.pose.get_cst_file()
        self.fix_res = self.pose.get_res_to_fix()

        #Define which xml and flags file to use
        if self.symm == 'SYMMETRIC':
            self.script_path = protocols_path + 'design_symm_auto.xml'
        else:
            self.script_path = protocols_path + 'design_auto.xml'

        self.flags_path = flags_path + 'flags_classic'


    def __str__(self):
        pass

    def design(self):
        all_thresholds_path = []

        for file in self.in_resfiles:
            os.system('mkdir ' + self.path_design + '/' + file)
            all_thresholds_path.append(self.path_design + '/' + file + '/pdbs')
            self.resfile_full_path = self.pose.get_path_fs() + '/resfiles/' + file

            job_command = 'sh ' + gen_path + 'create_job_auto.sh ' \
                            + self.path_design + '/' + file \
                            + ' -s ' + self.pose.get_path_refinement() + '/pdbs/' + self.pdb_input \
                            + ' -parser:protocol ' + self.script_path \
                            + ' -parser:script_vars in_resfile=' + self.resfile_full_path \
                            + ' -parser:script_vars keep_n=' + self.keep_n \
                            + ' -parser:script_vars res_to_fix=' + self.fix_res \
                            + ' -parser:script_vars pdb_reference=' + self.pdb_full_path \
                            + ' -parser:script_vars cst_full_path=' + self.cst_path \
                            + ' -parser:script_vars cst_value=' + str(self.pose.get_coorcst_value()) \
                            + ' -parser:script_vars pssm_full_path=' + self.pssm_path \
                            + ' @' + self.flags_path \
                            + ' -overwrite '

            if self.symm == 'SYMMETRIC':
                add_to_command = ' -parser:script_vars symmetry_file=' + self.pose.get_symm_full_path()
                job_command = job_command + add_to_command

            print 'About to submit ' + self.pdb_input + ' to design with in_resfile ' + file

            self.output_files = []
            for job in range(1,3):
                os.system(job_command + ' -out:suffix _' + str(job))
                name = self.pdb_input[:-4] + '_' + str(job) + '_0001.pdb'
                self.output_files.append(name)

            os.system('sh ' + self.path_design + '/' + file + '/' + '/command')

        repeats = 30000
        for threshold in all_thresholds_path:
            sleep_and_check('Design', repeats, self.output_files, threshold + '/')

        #Change the output pdb names. numbers 0-8 are FS thresholds. 0 is 0.5, 1 is -0.45 and so forth.
        for i in range (0, len(all_thresholds_path)):
            for j in range (0, len(self.output_files)):
                abc_vec = ['a','b','c','d','e','f','g','h']
                old_name = all_thresholds_path[i] + '/' + self.output_files[j]
                new_name = self.pose.get_pdbid() + '_designed_' + str(i + 1) + '_' + abc_vec[j] + '.pdb'
                new_name = all_thresholds_path[i] + '/' + new_name
                os.system('mv ' + old_name + ' ' + new_name)

