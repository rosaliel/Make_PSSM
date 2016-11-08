from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import urllib,os,re, subprocess, sys
import linecache

class PSSM:

    def __init__(self):
        pass

    def create_PSSM_from_MSA(self, path_to_msa, msa_file_name, path_to_head_seq, head_file_name, output_path = None):

        if 'sub_MSA_for_pos_' in msa_file_name:
            self.pssm_dir = path_to_msa + '/' + 'sub_PSSMs'
            output_pssm_name = 'sub_PSSM_' + msa_file_name[12:]

        #elif ('.aln' in msa_file_name or '.fa' in msa_file_name) and 'query' in msa_file_name:
        else:
            output_pssm_name = 'query_full_pssm'
            if output_path == None:
                print 'ERROR(PSSM_class):YOU MUST PROVIDE A PATH TO THE OUTPUT PSSM FILE!!!'
                print 'OUTPUT_PATH CAN BE None ONLY WHEN THE INPUT MSA CONTAINS THE SUBSTR \"sub_MSA_for_pos\"'
            else:
                self.pssm_dir = output_path

        # else:
        #     output_pssm_name = msa_file_name + '_PSSM'
        #     if output_path == None:
        #         print 'ERROR(PSSM_class): YOU MUST PROVIDE A PATH TO THE OUTPUT PSSM FILE!!!'
        #         print 'OUTPUT_PATH CAN BE None ONLY WHEN THE INPUT MSA CONTAINS THE SUBSTR \"sub_MSA_for_pos\"'
        #     else:
        #         self.pssm_dir = output_path


        if not os.path.exists(self.pssm_dir): os.makedirs(self.pssm_dir)

        if self.pssm_dir != None:
            print ''
            print 'psiblast -subject ' + path_to_head_seq + '/' + head_file_name + ' -in_msa ' + path_to_msa + '/'\
                  + msa_file_name + ' -out_ascii_pssm ' + self.pssm_dir + '/' + output_pssm_name

            os.system('psiblast -subject '
                       + path_to_head_seq + '/' + head_file_name
                       + ' -in_msa ' + path_to_msa + '/' + msa_file_name
                       + ' -out_ascii_pssm ' + self.pssm_dir + '/' + output_pssm_name)

    def get_pssm_output_path(self):
        return self.pssm_dir
