import urllib,os,re, subprocess, sys
import time

from path_variables import *


def sleep_and_check (process_name, repeats, out_files_vec, files_path):
    #puts the script into repeating sleep periods of 2 seconds until all files are ready.
    print "Start : %s" % time.ctime()
    for i in range(0,repeats):
        exist_status = ''
        for file in out_files_vec:
            if not os.path.isfile(files_path + file):
                exist_status = 'not_ready'
        
	    if exist_status == 'not_ready':
	        time.sleep(2)
            else:
            	continue

    print "End : %s" % time.ctime()
    if exist_status == 'not_ready':
        print '\n' + process_name + ' FAILED (OUTPUT FILES MISSING)'
    else:
        print '\n' + process_name + ' results are ready'
