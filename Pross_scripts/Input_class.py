import urllib,os,re, subprocess, sys

class Input:

    def __init__(self):
        self.struct_args = []
        self.blast_args = []

    def get_input(self):  #make sure there are no blank lines in the flag file
        import argparse
        parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
        parser.add_argument('-pdbid', default=None, type=str)
        parser.add_argument('-uploaded_fasta', default=None, type=str)      #If user uploads his own pdb this field is obligatory
        parser.add_argument('-input_seq', default=None, type=str)
        parser.add_argument('-chain_id', default=None, type=str)
        parser.add_argument('-small_ligands', default=None, type=str)       # small molecule ligand names separated by commas
        parser.add_argument('-large_ligands', type=str, default=None)       # DNA/protein ligands chain_ids separated by commas
        parser.add_argument('-symm_interfaces', type=str, default=None)     # starts with symmetry type, then white space then specify chain/s. examples: C2, B or D3, B,C  
        parser.add_argument('-res_to_fix', type=str, default=None)          # residues numbers provided by user, separated by commas
        parser.add_argument('-idealize', type=str, default='No')
        parser.add_argument('-input_msa', type=str, default='No')           # option for user to provide his own MSA file - full path
        parser.add_argument('-refined_file', type=str, default='No')
        parser.add_argument('-cst_value', type=float, default=0.4)
        parser.add_argument('-evalue', type=float, default=0.0001)
        parser.add_argument('-min_id', type=float, default=34.0)
        parser.add_argument('-max_targets', type=int, default=1500)
        parser.add_argument('-coverage', type=float, default=0.6)

        args = vars(parser.parse_args())  #args is a map

        ligands = ['small_ligands', 'large_ligands']
        for lig in ligands:
            if args[lig]:
                args[lig] = args[lig].split(',')

        return args

    def arrange_input(self, args):
        self.struct_args.append(args['pdbid']) # 0
        self.struct_args.append(args['input_seq']) # 1
        self.struct_args.append(args['chain_id']) # 2

        # ligands = ['small_ligands', 'large_ligands']
        # for lig in ligands:
        #     if args[lig]:
        #         args[lig] = args[lig].split(',')

        self.struct_args.append(args['small_ligands']) # 3
        self.struct_args.append(args['large_ligands']) # 4
        self.struct_args.append(args['res_to_fix']) # 5
        self.struct_args.append(args['symm_interfaces']) # 6
        self.struct_args.append(args['idealize']) # 7
        self.struct_args.append(args['input_msa']) # 8
        self.struct_args.append(args['refined_file']) # 9
        self.struct_args.append(args['uploaded_fasta']) # 10
        self.struct_args.append(args['cst_value']) # 11


        return self.struct_args

    def blast_args_f(self, args):
        self.blast_args.append(args['evalue'])
        self.blast_args.append(args['min_id'])
        self.blast_args.append(args['max_targets'])
        self.blast_args.append(args['coverage'])

        return self.blast_args

#input_data = Input()
#print input_data.arrange_input(input_data.get_input())
