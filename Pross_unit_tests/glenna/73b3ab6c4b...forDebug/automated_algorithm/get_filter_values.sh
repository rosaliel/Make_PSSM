#!/bin/bash

parameter=$1
path_to_structs=$2

temp_file="temp_file"
>$path_to_structs$temp_file
file_name_end="_data"
output_file_name=$parameter$file_name_end

for i in $(ls -1 $path_to_structs | grep .pdb)
do
    file_name=$i
    stability_score_full=$(grep $parameter $path_to_structs$file_name | awk '{print $2}')
    echo $file_name $stability_score_full >>$path_to_structs$temp_file
done

sort -nk 2 $path_to_structs$temp_file >$path_to_structs$output_file_name

rm -rf $path_to_structs$temp_file
