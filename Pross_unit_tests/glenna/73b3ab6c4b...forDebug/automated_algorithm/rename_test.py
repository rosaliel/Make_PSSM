import urllib,os,re, subprocess, sys

gen_path = '/home/labs/fleishman/adig/ThermoStab_benchmark/targets_auto/4I5V_min28/design/'
output_files = ['4I5V_A_refined_1_0001.pdb', '4I5V_A_refined_2_0001.pdb', '4I5V_A_refined_3_0001.pdb']
pdbid = '4I5V'

print len(output_files)

all_thresholds_path = [gen_path + 'designable_aa_resfile.-0.45/pdbs', gen_path + 'designable_aa_resfile.-0.75/pdbs', gen_path + 'designable_aa_resfile.-1/pdbs']

for i in range (0, len(all_thresholds_path)):
    for j in range (0, len(output_files)):
	print 'i = ' + str(i) + ' j = ' + str(j)
        abc_vec = ['a','b','c','d','e','f','g','h']
        old_name = all_thresholds_path[i] + '/' + output_files[j]
        new_name = pdbid + '_designed_' + str(i) + '_' + abc_vec[j]
        new_name = all_thresholds_path[i] + '/' + new_name
        os.system('mv ' + old_name + ' ' + new_name)


