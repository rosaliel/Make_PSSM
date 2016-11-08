#!/bin/sh

# ================================================================================================== #
#    The following script is used to detect fiterscan output mutations that are potential artifacts.
#    
#    Potential artifacts are mutations for which the main contribution to the total energy comes 
#    from energy terms like omega, fa_dun and maybe others. 
#
#    Script inputs, actions and outputs:	
# 1  Gets an input file of score_log file format.
# 2  User must enter a numerical cutoff value - mutations with total delta energy scores above
#    the cutoff are not analysed. 
# 3  For each mutation the script generates the names of the matching pdb files that were dumped 
#    during the filterscan run (wt file and mutant file)
# 4  Extract all the total energy terms of the pose (for the wt and mut pdb files)
# 5  Substract macthing terms (mutant - wt) to get the DELTA in that energy term  
# 6  2 output files are generated in this process:
#    a. report_file - contains all delta values for all the mutations below cutoff (2).
#    b. potential_artifacts - contains all mutations that have delta values below -1 that 
#       are not in fa_atr, fa_rep, fa_sol, fa_elec. The latter 4 are indicative of "real"
#	mutations so I filter them out to detect more easily the potential artifacts. 
# 7  The delta value used in (6) is defined in the variable $suspect_artifact.
# ================================================================================================== #

REU_max_value=$1
path_to_analysis=$2
path_to_fs=$3
pdb_prefix=$4

declare -i num_of_lines=$(cat $path_to_analysis/filtered_score_file | wc -l)
suspect_artifact=-1

>$path_to_analysis/report_file
>$path_to_analysis/potential_artifacts

echo -e "All mutations that have terms with Delta(mut-wt) that is lower than $suspect_artifact are shown in this file unless the terms are:" >>$path_to_analysis/potential_artifacts
echo -e "fa_atr, fa_rep, fa_elec, fa_sol, pssm\n" >>$path_to_analysis/potential_artifacts

declare -A aa_dict_1to3=( ["A"]="ALA" ["C"]="CYS" ["D"]="ASP" ["E"]="GLU" ["F"]="PHE" ["G"]="GLY" ["H"]="HIS" ["I"]="ILE" ["K"]="LYS" ["L"]="LEU" ["M"]="MET" ["N"]="ASN" ["P"]="PRO" ["Q"]="GLN" ["R"]="ARG" ["S"]="SER" ["T"]="THR" ["V"]="VAL" ["W"]="TRP" ["Y"]="TYR" )

for ((i=1; i!=$num_of_lines+1; i++))
	do
		delta_score=$(sed -n -e ${i}p $path_to_analysis/filtered_score_file | awk '{print $4}')
		if [ $(bc <<< "$delta_score<=$REU_max_value") -eq 1 ] #REU_max_value is input by the user
			then 
			res_num=$(sed -n -e ${i}p $path_to_analysis/filtered_score_file | awk '{print $1}')
			wt_id=$(sed -n -e ${i}p $path_to_analysis/filtered_score_file | awk '{print $2}')
			mut_id=$(sed -n -e ${i}p $path_to_analysis/filtered_score_file | awk '{print $3}')

			ref_pdb=$(echo "$pdb_prefix ${aa_dict_1to3[$wt_id]} $res_num ${aa_dict_1to3[$wt_id]} .pdb" | sed 's/ //g')
			mut_pdb=$(echo "$pdb_prefix ${aa_dict_1to3[$wt_id]} $res_num ${aa_dict_1to3[$mut_id]} .pdb" | sed 's/ //g')

			full_score_wt=$(grep TOTAL_SCORE $path_to_fs/pdbs/$ref_pdb | awk '{print $2}')
                        full_score_mut=$(grep TOTAL_SCORE $path_to_fs/pdbs/$mut_pdb | awk '{print $2}')
			delta_full_score=$(echo $full_score_mut"-("$full_score_wt")" | bc)
			
			echo "##################### #########" >>$path_to_analysis/report_file
			echo $wt_id$res_num$mut_id $delta_score >>$path_to_analysis/report_file
			echo "##################### #########" >>$path_to_analysis/report_file
			echo total $delta_full_score >>$path_to_analysis/report_file

			sed -n -e ${i}p $path_to_analysis/filtered_score_file >>$path_to_analysis/potential_artifacts
			
			for ((j=1; j<39; j++))
       				 do
      				 	if [[ $((j % 2)) = 0 ]]
                        			then
                        			energy_term=$(grep TOTAL_WTD $path_to_fs/pdbs/$ref_pdb | awk '{print $'$j'}')
						j=$(( j+1 )) #increment j by 1 to get an odd number
						term_value_wt=$(grep TOTAL_WTD $path_to_fs/pdbs/$ref_pdb | awk '{print $'$j'}')
                                                term_value_mut=$(grep TOTAL_WTD $path_to_fs/pdbs/$mut_pdb | awk '{print $'$j'}')
						delta_term=$(echo $term_value_mut"-("$term_value_wt")" | bc)
						echo -e "$energy_term \t$delta_term  " >>$path_to_analysis/report_file
						
						if [ $(bc <<< "$delta_term<=$suspect_artifact") -eq 1 ] # && [ $energy_term != "fa_atr:" ] && [ $energy_term != "fa_rep:" ] && [ $energy_term != "fa_sol:" ] && [ $energy_term != "fa_elec:" ] && [ $energy_term != "res_type_constraint:" ]
							then
							echo -n "$energy_term $delta_term  " >>$path_to_analysis/potential_artifacts
						fi
               				fi
			        done

				echo -e "\n" >>$path_to_analysis/report_file
				echo -e "\n" >>$path_to_analysis/potential_artifacts
		fi
	done

sed 's/coordinate_constraint/coorcst/g' $path_to_analysis/report_file | sed 's/res_type_constraint/pssm/g' >$path_to_analysis/temp1
column -t $path_to_analysis/temp1 >$path_to_analysis/report_file
sed -i '1s/^/All mutations from filterscan with all the energy_terms extracted from their pssm and substratced (mut-wt)"\n/' $path_to_analysis/report_file #add text at the 1st line of the file

rm -rf $path_to_analysis/temp1
