from Bio.SubsMat import MatrixInfo as matlist
from prody import *
from pylab import *
import urllib,os,re, subprocess, sys

from path_variables import *

ion() 

class Get_input:
    #############################################################################################
    # The class gets the following input from a user that wants to thermo-stabilize a protein:  #
    # pdb structure                                                                             #
    # protein sequence (in case the pdb structure is of a homoloug and not the query itself)    #
    # Chain of interest - in the future choose it according to molprobity                       #
    # Small molecules or atoms that are in proximity to the chain of interest                   #
    # Other chains that interact with chain of interest. These could be protein or dna chains.  #
    #                                                                                           #
    # Need to add - specific residues that the user wants to fix                                #
    #############################################################################################

    def __init__(self):
        self.all_info = [None, None, None, None, None, None]

    def get_user_input(self):
        #This function gets from the user the pdb code, protein_seq (if homoloug) and chain to design.
        #Calls the helper function get_ligands() to get a list of all ligands that interact with the chain of interest
        struct_or_homoloug = raw_input('\nTo stabilize a protein that has a pdb structure press Y else press N: ')

        if struct_or_homoloug == 'Y':
            pdbid = raw_input('\nPlease enter a 4 letters pdb code (for example 4EY7): ')
            input_seq = None
        elif struct_or_homoloug == 'N':
            input_seq = raw_input('\nPlease enter the protein sequence (remove tags): ')
            pdbid = raw_input('\nPlease enter a 4 letters pdb code of a homologous structure: ')
        else:
            print 'Wrong input\n'
            return 0

        chain_id = raw_input('\nWhich chain would you like to design? ')

        ligands = self.get_ligands()
        res_to_fix_str = self.get_specific_residues()

        self.all_info[0] = pdbid
        self.all_info[1] = input_seq
        self.all_info[2] = chain_id
        self.all_info[3] = ligands[0]
        self.all_info[4] = ligands[1]
        self.all_info[5] = res_to_fix_str

    def get_ligands(self):
        ligand_types = {'small_molecule' : [], 'dna_or_protein' : []}
        ligands = []

        for lig_type in ligand_types:
            lig_check = raw_input('\nIs there a ' + lig_type + ' substrate that binds the chain of interest (Y/N)? ')
            while lig_check == 'Y':
                ligand_types[lig_type].append(raw_input('Please enter the ' + lig_type + ' name/chain as named in the pdb: '))
                lig_check = raw_input('\nIs there another ' + lig_type + ' substrate that binds the chain of interest (Y/N)? ')

            if lig_check != 'N':
                print 'Wrong input\n'
                return 0
        small_molecules = ligand_types['small_molecule']
        dna_or_pro = ligand_types['dna_or_protein']

        ligands.append(small_molecules)
        ligands.append(dna_or_pro)
        return ligands

    def get_specific_residues(self):
        res_check = raw_input('\nAre there any specific residues you would like to fix (Y/N)? ')
        if res_check == 'Y':
            res_to_fix_str = raw_input('Please enter a list of all residues to fix separated by commas (no spaces): ')
            return res_to_fix_str

    def get_all_data(self):
        return self.all_info #The different elements will be the input for the Thermostab_query class.

