import urllib,os,re, subprocess, sys
from string import maketrans   # Required to call maketrans function.

from path_variables import *

#This class can be used for the following:
# 1. Run blast and get output as xml
# 2. Cluster an xml output with cd-hit
# 3. Make MSA from a fasta formatted file
# 4. Pairwise alignement -  for this option provide the 5th argument to the function and as full path

class Blast:

    def __init__(self, infile_name, infile_path, evalue, min_id, max_targets, coverage, pairwise_file = 'No_file'):
        self.infile_name = infile_name
        self.infile_path = infile_path
        self.infile_full = infile_path + infile_name
        self.pairwise_file = pairwise_file
        self.note_file = self.infile_path + 'blast_notes'
        self.num_of_seqs_file = self.infile_path + 'how_many_seqs'
        with open(self.note_file, 'a') as notes:
            notes.write('Notes from Blast run\n')
            notes.close()

        if coverage is 'default':
            self.min_coverage = 0.6
        else:
            self.min_coverage = coverage

        if min_id is 'default':
            self.min_id = 34.0 #at 30% sequence identity we can still assume very similar structures.
        else:
            self.min_id = min_id

        if evalue is 'default':
            self.evalue = str(0.0001)
        else:
            self.evalue = str(evalue)

        if max_targets is 'default':
            self.max_targets = str(1500)
        else:
            self.max_targets = str(max_targets)

    def run_blast(self):
        self.outfile =  self.infile_full + '_out_xml'
        print 'About to run blastp, for more details see output file blast_notes'
        with open(self.note_file, 'a') as notes:
            notes.write('About to run blastp with the following parameters:\n')
            notes.write('evalue = ' + self.evalue + \
                        '  max_targets = ' + self.max_targets + \
                        '  min_id = ' + str(self.min_id) + \
                        '  coverage = ' + str(self.min_coverage) + '\n')
            notes.close()

        os.system('/apps/RH6U4/ncbi-blast+/2.2.29/bin/blastp -query ' + self.infile_full \
                  + ' -db /shareDB/nr/nr' \
                  + ' -evalue ' + self.evalue \
                  + ' -max_target_seqs ' + self.max_targets \
                  + ' -outfmt 5'\
                  + ' -out ' + self.outfile)

        return self.outfile

    def run_psiblast(self):
        self.outfile =  self.infile_full + '_out_xml'
        print 'About to run psiblast'
        os.system('/apps/RH6U4/ncbi-blast+/2.2.29/bin/psiblast -query ' + self.infile_full \
                  + ' -db /shareDB/nr/nr' \
                  + ' -evalue ' + self.evalue \
                  + ' -max_target_seqs ' + self.max_targets \
                  + ' -num_iterations 5'
                  + ' -outfmt 5'\
                  + ' -out ' + self.outfile)

        return self.outfile

    def blast2seq_p(self):
        if self.pairwise_file == 'No_file':
            print 'ERROR: CAN\'T RUN BLAST2SEQ WITH NO FASTA FORMATTED PAIRWISE FILE CONTAINING SEQS TO ALIGN WITH QUERY'
        else:
            self.bl2seq_outfile =  self.infile_full + '_bl2seq_xml'
            os.system('/apps/RH6U4/ncbi-blast+/2.2.29/bin/blastp -query ' + self.infile_full \
                  + ' -subject ' + self.pairwise_file \
                  + ' -outfmt 5'\
                  + ' -out ' + self.bl2seq_outfile)

        return self.bl2seq_outfile

    def parse_blast_output(self, blast_xml_path):
        print 'About to parse blast xml output'
        with open(blast_xml_path, "r") as file:
            blast_data=file.read()

        hit_blocks = blast_data.split('BlastOutput_iterations>')
        if len(hit_blocks) != 3:
            print 'ERROR (Blast): XML FORMAT IS UNUSUAL, NUM OF BLOCK IS NOT 3'
        self.hits = hit_blocks[1].split('<Hit>')

        if not self.hits[0].rstrip('\n'):  #first element is an empty row with new line at the end.
            self.hits.pop(0)

        return self.hits

    def parse_xml_hit(self, hits_list):
        self.all_hits_data = []
        for hit in hits_list:
            hit_data = []
            hit_lines = hit.split('\n')
            for line in hit_lines:
                line = line.split('>')
                if '<Hit_id' in line[0] or '<Hsp_evalue' in line[0] or '<Hsp_identity' in line[0] \
                    or '<Hsp_align-len' in line[0] or '<Hsp_hseq' in line[0]:
                    if '<Hsp_hseq' in line[0]:
                        chars_to_remove = ['-']
                        line[1] = line[1].translate(None, ''.join(chars_to_remove))
                    hit_data.append(line[1].split('<')[0])

            #calculate the % identity by dividing <Hsp_identity> in the length of the hit sequence (see denominator) and *100
            if len(hit_data) == 5:
                hit_identity_percent = (float(hit_data[2])/float(hit_data[3]))*100
                hit_data.append(hit_identity_percent)

            #Checking that the XML format is correct: we expect 8 elements in each hit_data vector, see loop above
            if len(hit_data) == 6:
                self.all_hits_data.append(hit_data)
            elif len(hit_data) == 0:   #special condition for the 1st parsed element in xml output of blast2seq
                continue
            else:
                print 'WARNING (Blast): HIT IS UNUSUAL, NUM OF ELEMENTS IN HIT IS NOT 8 BUT ' + str(len(hit_data)) + ' ' + str(hit_data[0])
                print 'THIS MIGHT RESULT FROM FEW DIFFERENT ALIGNMENT FOR THE SAME SEQUENCE (MORE THAN 1 HSP PER HIT)'

        return self.all_hits_data

    def extract_query_seq(self):
        self.query_seq = ''
        with open (self.infile_full, "r") as file:
            infile_lines = file.read().split('\n')
            for line in infile_lines:
                if not '>' in line:
                    self.query_seq += line

        return self.query_seq

    def output_in_fasta(self):
        blast_data = self.run_blast()
        hit_blocks = self.parse_blast_output(blast_data)
        parsed_hits = self.parse_xml_hit(hit_blocks)
        query_seq = self.extract_query_seq()
        query_len = len(query_seq)
        coverage = self.min_coverage*query_len
        with open(self.note_file, 'a') as notes:
            notes.write('query_length is ' + str(query_len) + '\n')
            notes.write(str(self.min_coverage*query_len) + '\n')
            notes.close()

        count_fail = 0
        self.fasta_outfile = self.infile_full + '_out_filt_fa'
        f = open(self.fasta_outfile,'w')
        f.write('>Input_pdb_SEQRES|THIS_IS_QUERY' + '\n')
        f.write(query_seq + '\n')
        for hit in parsed_hits:
            if float(hit[5]) >= float(self.min_id) and float(hit[3]) >= float(coverage):
                id_2dig = '%.2f' % hit[5]
                f.write('>' + hit[0] + '| Identities = ' + str(id_2dig) + '%' '\n')
                f.write(hit[4] + '\n')
            else:
                with open(self.note_file, 'a') as notes:
                    notes.write(str('hit failed ' + hit[0] + ' ' + str(hit[5]) + ' ' + str(hit[3])) + '\n')
                    notes.close()
                count_fail += 1 #counts the sequences that do not pass the min_id threshold.
        with open(self.note_file, 'a') as notes:
            notes.write(str(count_fail) + ' hits are excluded due to sequence identity or coverage')
            notes.close()
        f.close()

        print 'THE NUMBER OF HITS IN THIS BLAST IS ' + str(len(parsed_hits))
        print str(count_fail) + ' SEQUENCES SHARE LESS THAN ' + str(self.min_id) + ' IDENTITIES WITH QUERY'
        print 'AND ARE EXCLUDED FROM THE FASTA OUTPUT'

        return self.fasta_outfile

    def run_cd_hit(self, infile, clust_threshold):
        print 'About to run cd-hit on initial blast results'
        if clust_threshold is 'default':
            clust_threshold = str(0.97)
        else:
            clust_threshold = str(clust_threshold)

        self.cd_hit_outfile = self.infile_full + '_clustered'

        os.system('/apps/RH6U4/cd-hit/4.6.1/bin/cd-hit -i ' + infile + ' -o ' + self.cd_hit_outfile + '_temp' + ' -c ' + clust_threshold)
        os.system('mv ' + self.cd_hit_outfile + '_temp.clstr ' + self.cd_hit_outfile + '.clstr')

        if 'Input_pdb_SEQRES' in open(self.cd_hit_outfile + '_temp').read():
            os.system('mv ' + self.cd_hit_outfile + '_temp ' + self.cd_hit_outfile)

        else:
            with open(self.cd_hit_outfile + '_temp', "r") as file_r:
                lines = file_r.read().split('\n')

            self.extract_query_seq()   # If implementing directly this function this line is needed.

            file_w = open(self.cd_hit_outfile, 'w')
            file_w.write('>Input_pdb_SEQRES|THIS_IS_QUERY' + '\n')
            file_w.write(self.query_seq + '\n')
            for line in lines:
                file_w.write(line + '\n')

            os.system('rm -rf ' + self.cd_hit_outfile + '_temp')

        return self.cd_hit_outfile

    def run_muscle_for_MSA_output(self, infile):
        print 'About to run muscle to get MSA'
        self.muscle_outfile = self.infile_path + 'query_msa.fa'     # In this format it's easy to move the query sequence to the top
        os.system('/apps/RH6U4/muscle/3.8.425/bin/muscle -in ' + infile + ' -out ' + self.muscle_outfile + '_temp')

        self.muscle_clustal = self.infile_path + 'query_msa.aln' # This format is good for user inspection
        #os.system(gen_path + 'fasta2clustal.pl ' + self.muscle_outfile + '_temp ' + '>' + self.muscle_clustal + '_temp')

        #Moving the query sequence to the top of the msa
        with open (self.muscle_outfile + '_temp', "r") as file:
            msa_lines=file.read().split('\n')

            for i in range (0, len(msa_lines)):
                if 'Input_pdb_SEQRES' in msa_lines[i]:
                    beg_line = i
                    while i < len(msa_lines):
                        if not ('>' in msa_lines):
                            i+=1
                            end_line = len(msa_lines)
                        else:
                            end_line = i
                            break
                    break

            f = open(self.muscle_outfile, 'w')
            for i in range (beg_line, end_line):       #print the query name and sequence first
                f.write(msa_lines[i] + '\n')
            for i in range (0, beg_line):              #then print all sequences that were originally before the query
                f.write(msa_lines[i] +'\n')
            for i in range(end_line, len(msa_lines)):  #then print all sequences that were originally after the query
                f.write(msa_lines[i] + '\n')
            f.close()

        #Getting rid of *** in the clw format.
        #f = open(self.muscle_clustal,'w')
        #with open (self.muscle_clustal + '_temp', "r") as file:
        #    msa_clw=file.read()
        #    for line in msa_clw.split('\n'):
        #        if not ("*" in line):
        #            f.write(line + '\n')
        #    f.close()

        os.system('rm -rf ' + self.muscle_outfile + '_temp')
        #os.system('rm -rf ' + self.muscle_clustal + '_temp')

    def get_seq_data_as_MSA(self):
        self.run_muscle_for_MSA_output(self.run_cd_hit(self.output_in_fasta(), 'default'))
        self.get_num_of_seqs_in_msa()

    def get_query_len(self):
        return len(self.query_seq)

    def get_num_of_seqs_in_msa(self):
        with open (self.cd_hit_outfile, "r") as file:
            cd_hit_lines = file.read().split('\n')
            file.close()

        self.seqs_in_msa = len(cd_hit_lines)/2

        with open(self.num_of_seqs_file, 'a') as file:
            file.write(str(self.seqs_in_msa))
            file.close()

        return self.seqs_in_msa


#path = query_path + '1HZY_/msa/'

#blast = Blast('blast_infile', path, 'default', 'default', 'default', 'default')
#cd_hit_file = blast.run_cd_hit(query_path + '1HZY_/msa/blast_infile_out_filt_fa', 0.97)
#blast.run_muscle_for_MSA_output(cd_hit_file)

#print blast.get_num_of_seqs_in_msa()

