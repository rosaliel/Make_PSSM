#!/bin/sh
#This script is used to unify many score_log_files (output of Filterscan.cc) into a single file
#and then filter out all the mutations above a cutoff energy that is given as input by the user.  

#This script should run from the directory with all the score_log_files (usually scores)

REU_max_value=$1
path_to_scores=$2
path_to_analysis=$3

>$path_to_analysis/full_score_file
>$path_to_analysis/filtered_score_file

for i in $(ls -1 $path_to_scores)
	do
		 cat $path_to_scores$i >>$path_to_analysis/full_score_file
	done

echo "All mutations with scores lower than $REU_max_value R.E.U." >>$path_to_analysis/filtered_score_file

for i in $(awk '{print $4}' $path_to_analysis/full_score_file)
	do
                num_to_grep=${i:1:8}; #the grep function cannot grep for "-" sign so I define num_to_grep as the reu without the -

		if [ $(bc <<< "$i <= $REU_max_value") -eq 1 ] #REU_max_value is input by the user 
			then
      			grep "${num_to_grep}" $path_to_analysis/full_score_file >>$path_to_analysis/filtered_score_file
        	fi
	done
               
rm -rf $path_to_analysis/full_score_file 
