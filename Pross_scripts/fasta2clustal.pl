#!/usr/bin/perl


#       fasta2clustal.pl
#       
#       Copyright 2011 Benjamin Tovar <scenesfromamemory4@gmail.com>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#
################################################################################
#
# DATE: 27/July/2011
# AUTHOR: Benjamin Tovar
#
# This program converts aligned and gapless FASTA file to
# a gapless CLUSTAL file, ideal for convert motifs to PWM.
#
################################################################################

use warnings;
use strict;

###### Print USAGE if no argument is given
my $usage = "\nUSAGE: fasta2clustal.pl <input_fasta_file>
EXAMPLE: fasta2clustal.pl dna_sequence.fa\n\n";

###### Get the arguments:  
my $input_fasta_file = shift or die $usage;

###### Let the computer decide:

	my $dna_sequence ='';
	my $sequence_name='';
	my $output_file_name='';

	open(INPUT_FILE,$input_fasta_file);
	my @file = <INPUT_FILE>;
	close INPUT_FILE;

		# ClUSTAL HEADER
		print "CLUSTAL W 2.0 multiple sequence alignment\n\n";

	foreach my $line(@file){
		
		# Discard empty lines
		if($line =~ /^\s*$/){
		next;
		}
		
		else{
				# Extract sequence names
				if($line =~ /^>/) {
					$sequence_name = $line;
					$sequence_name =~ s/\s//g;
					$sequence_name =~ s/\>//g;
					$sequence_name =~ s/\Start_position/_Start_position/g;
					$sequence_name =~ s/\:/_/g;
					chomp $sequence_name;
				}

				# Extract DNA sequence
				else{
				$dna_sequence = $line;
				chomp $dna_sequence;
		# Print results
                printf("%-50s %10s \n", "$sequence_name","$dna_sequence");
				} 
		}
	}

# Powred by #!CrunchBang Linux
# Benjamin Tovar

exit;
