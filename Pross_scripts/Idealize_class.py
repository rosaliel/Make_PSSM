import urllib,os,re, subprocess, sys
import time

from Thermostab_query import *
from Choose_best_structure_class import *
from sleep import *
from path_variables import *

class Idealize:

    def __init__(self, query, script_path = None, flags_user = None):
        self.query = query
        if script_path:
            self.script_path = script_path
        else:
            self.script_path = protocols_path + 'idealize.xml'

        if flags_user:
            self.flags_path = flags_user
        else:
            self.flags_path = flags_path + '/flags_classic'

    def __str__(self):
        pass

    def idealize(self):

        pdb_path = self.query.get_path_in()
        pdb_input = self.query.get_pdb_renumbered()
        ideal_path = self.query.get_path_idealize()

        job_command = 'sh ' + gen_path + 'create_job_auto.sh ' \
                        + ideal_path \
                        + ' -s ' + pdb_path + '/' + pdb_input \
                        + ' -parser:protocol ' + self.script_path \
                        + ' @' + self.flags_path \
                        + ' -overwrite '

        print '\nAbout to submit ' + pdb_input + ' to idealization'

        self.output_files = []
        for job in range(1,4):
            os.system(job_command + ' -out:suffix _' + str(job))
            name = pdb_input[:-4] + '_' + str(job) + '_0001.pdb'
            self.output_files.append(name)

        #run the command file
        os.system('chmod +x ' + ideal_path + '/command')
        os.system(ideal_path + '/command')

        repeats = 200
        sleep_and_check('Idealization', repeats, self.output_files, self.query.get_path_idealize() + '/pdbs/')


    def best_struct(self):
        self.output_path = self.query.get_path_idealize() + '/pdbs/'
        best = Choose_best_structure(self.output_path)
        best.sort_struct_based_on_filter_values('stability_pure')
        best_struct_name = best.get_best_structure('stability_pure_data')
        self.best_struct_name = best_struct_name[:7] + 'ideal.pdb'
        os.system('mv ' + self.output_path + best_struct_name + ' ' + self.output_path + self.best_struct_name)
        return self.best_struct_name


    def get_best_name(self):
        return self.best_struct_name

