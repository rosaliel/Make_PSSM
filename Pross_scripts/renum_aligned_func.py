import urllib,os,re, subprocess, sys

def renum_aligned_seqs(seq_with_gaps):
    ori_num_2_new_num = {}
    counter = 0
    for res_num in range (0, len(seq_with_gaps)):
        key = res_num
        if seq_with_gaps[res_num] == '-':
            counter += 1
            ori_num_2_new_num[key] = '-'
        else:
            ori_num_2_new_num[key] = res_num - counter

    return ori_num_2_new_num

