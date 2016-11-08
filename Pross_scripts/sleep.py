###################################################################
# the function puts the python scripts into repeating 2 seconds   #
# sleep periods to allow all Rosetta jobs to end.                 #    
# An special care for Filterscan since out file are generated     #
# at the beginning of the run so an extra condition is required.  #
###################################################################

import urllib,os,re, subprocess, sys
import time

from path_variables import *
from Thermostab_query import *

def sleep_and_check (process_name, repeats, out_files_vec, files_path):

    print "Start : %s" % time.ctime()

    for i in range(0,repeats):
        if os.path.isfile(stop_flag):
            sys.exit()

        status = 'not_ready'
        for file in out_files_vec:
            #The if takes care of refinement and combinatorial  design jobs
            if process_name != 'Filterscan':
                if not os.path.isfile(files_path + file):
                    time.sleep(2)
                    break
                else:
                    if file == out_files_vec[-1]:
                        status = 'ready'
                        break

            elif process_name == 'Filterscan':
                line = subprocess.Popen(['grep', 'reported success', files_path + file], stdout=subprocess.PIPE).stdout.read()    
                if not 'success' in line:
                    time.sleep(2)
                    break
                else:
                    if file == out_files_vec[-1]:
                        status = 'ready'
                        break
        
        if status == 'ready':
            break

    print "End : %s" % time.ctime()

    if status == 'not_ready':
        print '\n' + process_name + ' FAILED (OUTPUT FILES MISSING)'
    else:
        print '\n' + process_name + ' results are ready'
