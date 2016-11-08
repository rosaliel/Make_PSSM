from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import urllib,os,re, subprocess, sys

from Align2seq_class import *
from path_variables import *

# seq1 = 'MIRIKKKLILTIIYIHLFILNRLSFENAIKKTKNQENNLTLLPIKSTEEEKDDIKNGKDKKEIDNDKENIKTNNAKDHSTYIKSYLNTNVNDGLKYLFIPSHNSFIKKYSVFNQINDGMLLNEKNDVKNNEDYKNVDYKNVNFLQYHFKELSNYNIANSIDILQEKEGHLDFVIIPHYTFLDYYKHLSYNSIYHKSSTYGKCIAVDAFIKKINETYDKVKSKCNDIKNDLIATIKKLEHPYDINNKNDDSYRYDISEEIDDKSEETDDETEEVEDSIQDTDSNHTPSNKKKNDLMNRTFKKMMDEYNTKKKKLIKCIKNHENDFNKICMDMKNYGTNLFEQLSCYNNNFCNTNGIRYHYDEYIHKLILSVKSKNLNKDLSDMTNILQQSELLLTNLNKKMGSYIYIDTIKFIHKEMKHIFNRIEYHTKIINDKTKIIQDKIKLNIWRTFQKDELLKRILDMSNEYSLFITSDHLRQMLYNTFYSKEKHLNNIFHHLIYVLQMKFNDVPIKMEYFQTYKKNKPLTQ'
# seq2 = 'DKSSITKINEDIEKFNEEIIKNEEQCLVGGKTDFDNLLIVLENAEKANVRKTLFDNTFNDYKNKKSSFYNCLKNKKNDYDKKIKNIKNEITKLLKNIESTGNMCKTESYVMNNNLYLLRVNEVKSTPIDLYLNRAKELLESSSKLVNPIKMKLGDNKNMYSIGYIHDEIKDIIKRYNFHLKHIEKGKEYIKRITQANNIADKMKKDELIKKIFESSKHFASFKYSNEMISKLDSLFIKNEQ'
#
# pairwise = Align2seq(seq1, seq2)
#
# pairwise.align2seq()
# aligned_1 = pairwise.get_aligned_seq1()
# aligned_2 = pairwise.get_aligned_seq2()
#
# print pairwise.get_identity_percent(aligned_1, aligned_2)


def read_input(input_path):
    all_lines = open(input_path, 'r')
    return all_lines

all_lines = read_input(gen_path + 'blast_results_edited')

name2seq_map = {}
keys_vec = []

for line in all_lines:
    if line.startswith('>'):
        key = line.strip()
        keys_vec.append(key)
    elif key in name2seq_map:
        name2seq_map[key] += line.strip()
    else:
        name2seq_map[key] = line.strip()


for i in range(0, len(keys_vec)):
    print keys_vec[i]

seqs_to_remove = []

print 'Length of dictionary = ' + str(len(name2seq_map))

for i in range(0, len(keys_vec)):
    for j in range(i + 1, len(keys_vec)):
        if keys_vec[i]: #not in seqs_to_remove:
            print keys_vec[i][0:20], keys_vec[j][0:20]
            pairwise = Align2seq(name2seq_map[keys_vec[i]], name2seq_map[keys_vec[j]])
            pairwise.align2seq()
            aligned_1 = pairwise.get_aligned_seq1()
            aligned_2 = pairwise.get_aligned_seq2()
            identity = pairwise.get_identity_percent(aligned_1, aligned_2)
            if identity > 99.7 and keys_vec[j] not in seqs_to_remove:
                seqs_to_remove.append(keys_vec[j])

print len(seqs_to_remove)

for seq in seqs_to_remove:
    print seq[0:20]
    print name2seq_map[seq]
    del name2seq_map[seq]
    print 'I am here'

print len(name2seq_map)


            # print aligned_1[0:200]
            # print aligned_2[0:200]
            # print ''
            # print aligned_1[200:400]
            # print aligned_2[200:400]
            # print ''
            # print aligned_1[400:]
            # print aligned_2[400:]
            # print ''
