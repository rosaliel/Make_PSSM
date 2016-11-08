#!/bin/bash

##############################################################################
# The following script removes lines from a PSSM file according to input     #
# The script gets 3 inputs:                                                  # 
# 1. a string containing residue numbers to remove seperated by commas       #
# 2. a string containing single letter identities matching the 1st input     #
# 3. a full path to a PSSM file from which to remove the residues            #
##############################################################################

res_num_str=$1
res_id_str=$2
pssm_path=$3

echo "#!/bin/bash"

cmd_start="cat "
add_to_cmd=$cmd_start$pssm_path
general_addition=" | grep -v "
quot="\""
space="    "
counter=0

for i in $(echo $res_num_str | sed 's/,/ /g');
do
    counter=$(($counter+1))
    res_id=$(echo $res_id_str | sed 's/,/ /g' | awk -v num=$counter '{print $num}')
    add_to_cmd=$add_to_cmd$general_addition$quot$space$(echo $i $res_id)$quot
done

echo $add_to_cmd
