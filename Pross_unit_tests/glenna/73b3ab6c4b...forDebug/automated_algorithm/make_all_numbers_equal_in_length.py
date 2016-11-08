import urllib,os,re, subprocess, sys

def make_all_numbers_equal_in_length(num_to_handle, longest_num):
    num_to_str = str(num_to_handle)
    num_to_str_len = len(num_to_str)
    longest_num_len = len(str(longest_num))
    zeros_to_add_before = ''

    for i in range (num_to_str_len, longest_num_len):
        zeros_to_add_before += '0'

    num_to_str = zeros_to_add_before + num_to_str

    return num_to_str

