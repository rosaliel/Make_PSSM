import urllib,os,re, subprocess, sys
from path_variables import *

class Choose_best_structure:

    def __init__(self, path_to_structures):
        self.path_to_structures = path_to_structures

    def sort_struct_based_on_filter_values(self, filter_name):
        #Results of the next call are stored in a report file. To get the file's name use:
        #self.get_report_file_name(filter_name). To get the file's path use: self.get_report_file_path()
        os.system('sh ' + gen_path + 'get_filter_values.sh '
                   + filter_name + ' '
                   + self.path_to_structures)

    def calc_average_filter_value(self, report_file_name, report_file_path = None):
        if report_file_path == None:
            report_file_path = self.path_to_structures

        report_file = open(report_file_path + '/' + report_file_name, 'r')
        sum = 0
        counter = 0

        for line in report_file:
            columns = line.rstrip().split()
            sum += float(columns[1])
            counter += 1

        average = sum/counter

        return average

    def get_report_file_name(self, filter_name):
        output_file_name = filter_name + '_data'
        return output_file_name

    def get_report_file_path(self):
        #relevant to call only if file was produced with default path.
        return self.path_to_structures

    def does_stability_val_make_sense(self, average_val, chain_length):

        expected_val_per_aa = -1.5  #at least
        ave_val_per_aa = float(average_val/chain_length)
        if ave_val_per_aa > expected_val_per_aa:
            print 'WARNING(Choose_best_structure_class): STABILITY VALUES SEEM TO BE TOO HIGH!'
            print 'MANUAL EXAMINATION OF THE STRUCTURES IS RECOMMENDED'
        else:
            print 'STABILITY VALUES ARE REASONABLE'

        print 'AVERAGE STABILITY PER AMINO ACID = ' + str(ave_val_per_aa)

    def does_rmsd_val_make_sense(self, rmsd_val):

        expected_val_for_coorcst_0p4 = 0.3  #Note that it is expected for coorcst of 0.4
        if rmsd_val > expected_val_for_coorcst_0p4*1.6:
            print 'WARNING(Choose_best_structure_class): RMSD VALUES SEEM TO BE TOO HIGH!'
            print 'MANUAL EXAMINATION OF THE STRUCTURES IS RECOMMENDED'
        elif rmsd_val < expected_val_for_coorcst_0p4/1.6:
            print 'WARNING(Choose_best_structure_class): RMSD VALUES SEEM TO BE TOO LOW!'
            print 'MANUAL EXAMINATION OF THE STRUCTURES IS RECOMMENDED'
        else:
            print 'RMSD VALUES ARE REASONABLE'

        print 'AVERAGE RMSD VALUE = ' + str(rmsd_val)

    def get_best_structure(self, report_file_name, ave_stability = None, report_file_path = None):
        if report_file_path == None:
            report_file_path = self.path_to_structures

        with open(report_file_path + '/' + report_file_name, 'r') as f:
            first_line_cols = f.readline().rstrip().split()
            self.best_structure_name = first_line_cols[0]
            self.best_score = float(first_line_cols[1])
            if ave_stability != None:
                self.diff_from_ave = self.best_score - ave_stability


        print 'THE TOP SCORING STRUCTURE IS ' + self.best_structure_name + ' SCORE = ' + str(self.best_score)
        if ave_stability != None:
            print 'BEST STRUCTURE\'S SCORE DIFFERENCE FROM AVERAGE = ' + str(self.diff_from_ave)

        return self.best_structure_name


# path_to_refined_structures = '/home/labs/fleishman/adig/ThermoStab_benchmark/4EY7_/refinement/pdbs/'
# Ch = Choose_best_structure(path_to_refined_structures)
# Ch.sort_struct_based_on_filter_values('stability_score_full')
# Ch.sort_struct_based_on_filter_values('stability_pure')
# Ch.sort_struct_based_on_filter_values('rmsd')
#
# full_stability_file = Ch.get_report_file_name('stability_score_full')
# pure_stability_file = Ch.get_report_file_name('stability_pure')
# rmsd_file = Ch.get_report_file_name('rmsd')
#
# stability_full_ave = Ch.calc_average_filter_value(full_stability_file, path_to_refined_structures)
# stability_pure_ave = Ch.calc_average_filter_value(pure_stability_file, path_to_refined_structures)
# rmsd_ave = Ch.calc_average_filter_value(rmsd_file, path_to_refined_structures)
#
# Ch.does_stability_val_make_sense(stability_pure_ave, 530)
# Ch.does_stability_val_make_sense(stability_pure_ave, 700)
# Ch.does_stability_val_make_sense(stability_pure_ave, 1000)
#
# Ch.does_rmsd_val_make_sense(rmsd_ave)
#
# Ch.get_best_structure(full_stability_file, stability_full_ave)

