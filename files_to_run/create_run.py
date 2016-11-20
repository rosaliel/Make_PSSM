import sys
import os
import argparse

OUTPUT = '/home/labs/fleishman/rosaliel/programs/make_pssm/output'
SCRIPT = '/home/labs/fleishman/rosaliel/programs/make_pssm/Pross_scripts/main.py'
PATH = '/home/labs/fleishman/rosaliel/programs/make_pssm'
FLAGS = '/home/labs/fleishman/rosaliel/programs/make_pssm/files_to_run/flag_template'

def parse_args():
    """"""
    desc = ('Creates jobs and command file to run Adi\'s PSSMs. '
            'Run the file command_after when all(!!!) jobs are done to get the '
            'final PSSMs to the directory pssm/')            
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('pdbs', 
                        help=('List of pdb files. Names should be: xxxxY.pdb '
                              'where xxxx is the pdbID and Y is the chain'),
                        nargs='+',
			            action='store')
    return parser.parse_args()
    
def create_dir(path):
    if not os.path.exists(path): 
        os.mkdir(path)
  
def create_dirs(path):
    """"""
    create_dir(path + '/flags/')
    create_dir(path + '/out/')
    create_dir(path + '/err/')
    create_dir(path + '/jobs/')
    create_dir(path + '/pssm/')
    
def prepare_flag(path, pdb_path, name, chain):
    """"""
    flags_file = open(FLAGS, 'r').readlines()
    flags_path = '{}/flags/{}_flags'.format(path, name)
    flag = open(flags_path, 'a')
    for line in flags_file:
        if 'pdbid' in line:
            flag.write('-pdbid={}\n'.format(pdb_path))
        elif 'chain_id' in line:
            flag.write('-chain_id={}\n'.format(chain))
        else:
            flag.write(line)  
    flag.close()   
    return flags_path

def prepare_job(path, name, flags_path):
    """"""
    job_path = '{}/jobs/{}.job'.format(path, name)
    job = open(job_path, 'a')
    job.write('#!/bin/bash\n')    
    job.write('. /usr/share/lsf/conf/profile.lsf\n')    
    job.write('cd {}\n'.format(PATH))
    job.write('python {} @{}\n'.format(SCRIPT, flags_path))
    job.close()
    os.chmod(job_path, 0770) # executable job
    return job_path    

def prepare_command(path, job_path, name):
    """"""
    command_path = path + '/command'
    if not os.path.exists(command_path):
        open(command_path, 'a').close()
        os.chmod(command_path, 0770)
    command_file = open(command_path, 'a')
    command = 'bsub -u /dev/null -R rusage[mem=2048]  -q fleishman '
    command += '-o {}/out/{}.out '.format(path, name)
    command += '-e {}/err/{}.err '.format(path, name)
    command += job_path + '\n'
    command_file.write(command)

def prepare_command_after(path, name, output_path):
    """prepares a script to run after all jobs are done, to copy the output"""
    command_after_path = path + '/command_after'
    if not os.path.exists(path): 
        open(command_after_path, 'a').close()
        os.chmod(command_after_path, 0770)
    file = open(command_after_path, 'a')
    file.write('cp {}/{}_/in/final_pssm_for_Rosetta {}/pssm/{}.pssm\n'.format(
                                                 output_path, name, path, name))
    file.write('rm -rf {}/{}_\n'.format(output_path, name))                                              
    
def main():
    args = parse_args()
    pwd = os.getcwd()
    create_dirs(pwd)
    for pdb in args.pdbs:
        pdb = os.path.abspath(pdb)
        name = os.path.basename(pdb)
        pdbID = name[:4]
        chain = name[4]
        flags_path = prepare_flag(pwd, pdb, name[:-4], chain)
        job_path = prepare_job(pwd, name[:-4], flags_path)
        prepare_command(pwd, job_path, name[:-4])
        prepare_command_after(pwd, name[:-4], OUTPUT)
        

if __name__ == '__main__':
    main()
