from Bio.SubsMat import MatrixInfo as matlist
from prody import *
from pylab import *
import urllib,os,re, subprocess, sys

ion()

from Thermostab_query import *
from Idealize_class import *
from sleep import *
from path_variables import *

class Refinement:

    def __init__(self, pose, symm, ex1_t, ex2_t, flags_user = None):
        self.pose = pose
        self.symm = symm
        self.ex1_t = ex1_t
        self.ex2_t = ex2_t

        #In the next 2 if/else defines which of the 4 possible scripts to use.
        if not self.symm == 'NONE':
            self.script_path = protocols_path + 'refine_symm_'
        else:
            self.script_path = protocols_path + 'refine_'

        if self.pose.get_input_seq():
            self.script_path += 'thread_auto.xml'
        else:
            self.script_path += 'auto.xml'

        #Allows flags from user
        if flags_user:
            self.flags_path = flags_user
        else:
            self.flags_path = flags_path + 'flags_classic'

    def __str__(self):
        pass

    def define_input(self): #Executes idealization if desired and defines the input structure to be submitted to refinement.

        if self.pose.idealize_YorN() == 'No':
            self.input_path = self.pose.get_path_in() + '/'
            self.input_name = self.pose.get_pdb_renumbered()
        else:
            idealize = Idealize(self.pose)
            idealize.idealize()
            idealize.best_struct()
            self.input_path = self.pose.get_path_idealize() + '/pdbs/'
            self.input_name = idealize.get_best_name()

        return self.input_name

    def refine(self):

        self.define_input()
        job_command = 'sh ' + gen_path + 'create_job_auto.sh ' \
                        + self.pose.get_path_refinement() \
                        + ' -s ' + self.input_path + self.input_name \
                        + ' -parser:protocol ' + self.script_path \
			+ ' -parser:script_vars cst_value=' + str(self.pose.get_coorcst_value())\
                        + ' -parser:script_vars res_to_fix=' + self.pose.get_res_to_fix() \
                        + ' -parser:script_vars pdb_reference=' + self.input_path + self.input_name \
                        + ' -parser:script_vars cst_full_path=' + self.pose.get_cst_file() \
                        + ' @' + self.flags_path \
                        + ' -overwrite '

        #For large proteins having the flag ex2 might kill the run. 
        if self.pose.get_pdb_length() < self.ex1_t:
            add_to_command = '-ex1 '
            job_command = job_command + add_to_command
  
        if self.pose.get_pdb_length() < self.ex2_t:
             add_to_command = '-ex2 '
             job_command = job_command + add_to_command

        if self.pose.get_input_seq():
            add_to_command = ' -parser:script_vars input_seq=' + self.pose.get_input_seq() \
                           + ' -parser:script_vars cst_value=' + str(self.pose.get_coorcst_value())
            job_command = job_command + add_to_command

        if not self.symm == 'NONE':
            add_to_command = ' -parser:script_vars symmetry_file=' + self.pose.get_symm_full_path()
            job_command = job_command + add_to_command

        print '\nAbout to submit ' + self.input_name + ' to refinement 5 times'

        self.output_files = []
        for job in range(1,6):
            os.system(job_command + ' -out:suffix _' + str(job))
            name = self.input_name[:-4] + '_' + str(job) + '_0001.pdb'
            self.output_files.append(name)

        os.system('chmod +x ' + self.pose.get_path_refinement() + '/command')
        os.system(self.pose.get_path_refinement() + '/command')

        repeats = 30000
        sleep_and_check('Refinement', repeats, self.output_files, self.pose.get_path_refinement() + '/pdbs/')

    def best_struct(self):
        self.output_path = self.pose.get_path_refinement() + '/pdbs/'
        best = Choose_best_structure(self.output_path)
        best.sort_struct_based_on_filter_values('stability_score_full')
        best_struct_name = best.get_best_structure('stability_score_full_data')
        self.best_struct_name = best_struct_name[:7] + 'refined.pdb'
        os.system('mv ' + self.output_path + best_struct_name + ' ' + self.output_path + self.best_struct_name)
        return self.best_struct_name

    def get_best_name(self):
        return self.best_struct_name

    def get_best_path(self):
        return self.output_path

