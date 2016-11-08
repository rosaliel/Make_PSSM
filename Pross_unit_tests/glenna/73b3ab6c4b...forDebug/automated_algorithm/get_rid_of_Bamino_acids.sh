#! /bin/bash

path_to_pdb_file=$1
output_pdb_file_path=$2

for i in $path_to_pdb_file
    do
        egrep -v '^.{16}B' $i | egrep -v '^.{16}C' | awk '{print substr($0,1,16),substr($0,18,60)}'; done > $output_pdb_file_path
