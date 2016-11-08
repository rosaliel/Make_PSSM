from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from prody import *
from pylab import *
import urllib,os,re, subprocess, sys

class Align2seq:

    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.aligned_seq1 = ''
        self.aligned_seq2 = ''
        self.size = 0

    def align2seq(self):
        matrix = matlist.blosum62
        gap_open = -11
        gap_extend = -1

        alignment = pairwise2.align.globalds(self.seq1, self.seq2, matrix, gap_open, gap_extend)
        seq1_vs_seq2 = alignment[0]

        self.aligned_seq1, self.aligned_seq2, self.score, begin, self.size = seq1_vs_seq2

    def find_gaps(self, seq):
        self.gaps_vec = []
        for char_idx in range (0, len(seq)):
            if seq[char_idx] == '-':
                self.gaps_vec.append(char_idx)

        return self.gaps_vec

    def remove_gaps(self, seq, gaps_vec):
        counter = 0
        for gap in gaps_vec:
            seq = seq[:gap - counter] + seq[gap -counter + 1:]
            counter += 1

        return seq

    def get_identity_percent(self, aligned_1, aligned_2):
        if aligned_1[0:1] is '-':
            seq_to_scan = aligned_1
        elif aligned_2[0:1] is '-':
            seq_to_scan = aligned_2
        else:
            seq_to_scan = None

        counter = 0
        if seq_to_scan:
            for char in seq_to_scan:
                if char is '-':
                    counter += 1
                else:
                    break

        aligned_1 = aligned_1[counter:]
        aligned_2 = aligned_2[counter:]

        if aligned_1[-1:-2] is '-':
            seq_to_scan = aligned_1
        elif aligned_2[-1:-2] is '-':
            seq_to_scan = aligned_2
        else:
            seq_to_scan = None

        counter = 0  #Should stay before the if statement to reset counter back to zero
        if seq_to_scan:
            for char in reversed(seq_to_scan):
                if char is '-':
                    counter += 1
                else:
                    break

        if counter != 0:
            aligned_1 = aligned_1[:-counter]
            aligned_2 = aligned_2[:-counter]

        #Count the number of identities
        counter = 0.000
        for i in range (0, len(aligned_1)):
            if aligned_1[i] == aligned_2[i]:
                counter += 1.000

        identities = (counter/len(aligned_1))*100
        print 'Identities = ' + str(round(identities, 2)) + '%'
        return round(identities, 2)

    def get_aligned_seq1(self):
        return self.aligned_seq1

    def get_aligned_seq2(self):
        return self.aligned_seq2

    def get_alignment_length(self):
        return len(self.aligned_seq1)

    def get_gaps_vec(self):
        return self.gaps_vec

    def get_score(self):
        return self.score

