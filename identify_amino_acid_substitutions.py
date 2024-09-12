#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 15:19:57 2022

@author: Lizzie (elizabeth.robbins@biology.ox.ac.uk)

This script counts amino acid mutations that have occured throughout a phylogenetic tree.

File requirements: Species tree, gene tree (generated from IQ-TREE ASR), ancestral state reconstruction alignment (FASTA) and multiple sequence alignment (FASTA).


"""
#IMPORT LIBRARIES
import os
import dendropy
import pandas as pd
import numpy as np
from tqdm import tqdm
import ast
from pathlib import Path

#GENETIC CODES
genetic_code_11 = {'F':['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                   'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
                   'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'P': ['CCT','CCC', 'CCA', 'CCG'],
                   'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A':['GCT', 'GCC','GCA', 'GCG'], 
                   'Y': ['TAT', 'TAC'], '*': ['TAA', 'TAG', 'TGA'], 'C': ['TGT', 'TGC'],
                   'W': ['TGG'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                   'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
                   'G': ['GGT', 'GGC', 'GGA', 'GGG']}

GC_11_codons = {'TTT': 'F', 'TTC': 'F',
                'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG':'L',
                'ATT':'I', 'ATC':'I', 'ATA':'I',
                'ATG':'M',
                'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
                'CCT':'P','CCC':'P', 'CCA':'P', 'CCG':'P',
                'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                'GCT':'A', 'GCC':'A','GCA':'A', 'GCG':'A',
                'TAT':'Y', 'TAC':'Y',
                'TAA':'*', 'TAG':'*', 'TGA':'*',
                'TGT':'C', 'TGC':'C',
                'TGG':'W',
                'CAT':'H', 'CAC':'H',
                'CAA':'Q', 'CAG':'Q',
                'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
                'AAT':'N', 'AAC':'N',
                'AAA':'K', 'AAG':'K',
                'GAT':'D', 'GAC':'D',
                'GAA':'E', 'GAG':'E',
                'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

###FUNCTIONS
def get_files(working_directory):
    '''
    INPUT: working_directory - directory containing alignment files and the output of the IQ-TREE ASR. 
    
    RETURNS: str_treefile_pathlist, tr_statefile_pathlist and str_alignment_pathlist - lists of file paths strings sorted alphabetically.
    Should be able to index the list of strings to get the files associated with the various genes. 
    '''
    # read in files
    treefile_pathlist = Path(working_directory).rglob('*.aln.treefile')
    statefile_pathlist = Path(working_directory).rglob('*ancestral_seqs.nuc.aln')
    alignment_pathlist = Path(working_directory).rglob('*.nuc.aln')
    
    # make lists of tree files, ASR files and alignment files.
    str_treefile_pathlist = []
    str_statefile_pathlist = []
    str_alignment_pathlist = []
    for path in treefile_pathlist:
        path_str = str(path)
        str_treefile_pathlist.append(path_str)
    str_treefile_pathlist.sort()
    for path in statefile_pathlist:
        path_str = str(path)
        str_statefile_pathlist.append(path_str)
    str_statefile_pathlist.sort()
    for path in alignment_pathlist:
        path_str = str(path)
        if '_' not in path_str.split('/')[-1]:
            str_alignment_pathlist.append(path_str)
    str_alignment_pathlist.sort()
    
    return str_treefile_pathlist, str_statefile_pathlist, str_alignment_pathlist

def get_child_parent_relationships(treefile, outgroup_mrca, mrca_of_interest):
    '''INPUT: treefile - phylogenetic tree in newick format
              outgroup_mrca - list of species in outgroup (required to root the tree)
              mrca_of_interest - the mrca of the list of species is the root of the subtree of interest
       RETURNS: list of child -> parent relationships within the tree
    '''
    # read in phylogenetic tree
    t = dendropy.Tree.get(file=open(treefile, 'r'), schema="newick")
    # identify the most recent common ancestor (mrca) of the list of species in the outgroup_mrca list
    outgroup_mrca = t.mrca(taxon_labels=outgroup_mrca)
    # get the name of the mrca node
    outgroup_label = outgroup_mrca.label
    # re-root the tree at the node edge of the outgroup
    t.reroot_at_edge(outgroup_mrca.edge, update_bipartitions=False)
    outgroup_mrca.label = outgroup_label
    
    #initiate the list of child, parents of each branch in the tree 
    child_parent_list = []
    root = ''
    internal_node_count = 0
    taxon_count = 0
    
    # get the most recent common ancestor of the tree we want to count (note: there needs to be a different way of doing this)
    tree_mrca = t.mrca(taxon_labels=mrca_of_interest)
    
    #iterate through the nodes in the clade defined by mrca_of_interest
    for nd in tree_mrca.postorder_iter():
        if nd.label == tree_mrca.label:
            continue
        
        internal_node_count +=1
        child_parent = []
        
        #identify the root (has no child -> parent relationship)
        if nd.parent_node == None:
            root = nd.label
            continue
        
        #get child-parent relationships for internal nodes and leaves
        elif nd.parent_node != None: #confirming a true child -> parent relationship
            
            #child is an internal node (use nd.label for child)
            if nd.label != None: #internal nodes have a node label (leaves seemingly do not)
                child = nd.label
                child_parent.append(child)
                child_parent.append(nd.parent_node.label) #appeding parent label
            
            #child is a leaf (use nd.taxon.label for child)
            else:
                taxon_count += 1
                child = nd.taxon.label #get child -> parent relationships for leaves
                child_parent.append(child)
                child_parent.append(nd.parent_node.label)
        # add the child, parent relationship for the branch 
        child_parent_list.append(child_parent)
    print(f"The number of nodes in the tree analysed is {internal_node_count}") # = (n-1), where n = number of species 
    print(f"The number of leaves in the tree analysed is {taxon_count}") # number of species analysed
    return child_parent_list, root

def get_sequences(alignment_file):
    '''
    INPUT: alignment_file - contains a FASTA file of extanct sequences
    RETURNS: alignment dictionary - dictionary of species sequences. Key, species; value, sequence
    '''
    alignment_dict ={}
    with open(alignment_file, 'r') as f:
        whole_file = f.read()
        for entry in whole_file.split(">")[1:]:
            data = entry.split()
            species_name = data[0]
            species_name = species_name.replace('_', ' ')
            sequence_lines = data[1:]
            sequence_lines = ''.join(sequence_lines)
            alignment_dict[species_name] = sequence_lines
    return alignment_dict

def translate(seq):
    '''
    INPUT: nucleotide sequence
    OUTPUT: protein sequence
    '''
    #NCBI's genetic code 11 (for plant plastids)
    table = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 
                    'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
                    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 
                    'AGT':'S', 'AGC':'S', 'CCT':'P', 'CCC':'P', 
                    'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 
                    'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A',
                    'GCA':'A', 'GCG':'A', 'TAT':'Y', 'TAC':'Y',
                    'TAA':'*', 'TAG':'*', 'TGA':'*', 'TGT':'C', 
                    'TGC':'C', 'TGG':'W', 'CAT':'H', 'CAC':'H',
                    'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 
                    'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
                    'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
                    'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
                    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon not in table.keys():
                protein += '-'
            else:
                protein+= table[codon]
    return protein

def count_mutations(child_parent_list, sequence_dict, residue):
        '''
        INPUT: 
        RETURNS: site matrix (20x20) with all substitutions at a specific site
        '''
        
        #intitialise a mutation matrix for that sites
        alignment_residue = residue
        
        residues = ['C', 'S', 'T', 'A', 'G', 'P', 'D', 'E', 'Q', 'N', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'F', 'Y', 'W']
        m = pd.DataFrame(np.zeros(shape = (20,20)), index = residues, columns = residues)

        #iterate over every branch in the tree
        for edge in child_parent_list:
            identities = ['x', 'x']
            child = edge[0]
            parent = edge[1]
            
            child_ident = sequence_dict[child][alignment_residue]
            identities[0] = child_ident
            parent_ident = sequence_dict[parent][alignment_residue]
            identities[1] = parent_ident
            
            if identities[0] in residues and identities[1] in residues:
                m.loc[parent_ident, child_ident] += 1 #parent is the row and child is the column
            else:
                print(edge, identities)
                continue
                
        return m  
    
if __name__ == "__main__":
    
    # USER INPUTS
    # define outgroup for rooting the gene tree on an outgroup
    outgroup_mrca = ['Nuphar japonica','Nymphaea jamesoniana','Cabomba caroliniana',
             'Victoria cruziana','Nymphaea capensis','Nuphar advena',
             'Cabomba furcata','Cabomba aquatica','Euryale ferox','Nymphaea ampla','Nymphaea lotus']
    
    #defines species in the tree of interest (defines root node of subtree of interest that excludes the outgroup species)
    mrca_of_interest = ['Talipariti hamabo', 'Aegilops geniculata', 'Acorus americanus']
    
    #define directory for input files
    file_path= '/Users/Lizzie/Documents/Oxford_University/DTP/DPhil/Writing/Selection_paper/Post_TPC_review/Photosystem_adaptation_zenodo/test/'
    output_file_path = file_path
    
    # get lists of the files to analyse
    str_treefile_pathlist, str_statefile_pathlist, str_alignment_pathlist = get_files(file_path)
    
    for i in range(0, len(str_treefile_pathlist)): 
        file1 = (os.path.basename(str_treefile_pathlist[i])).split('.')[0]
        file2 = (os.path.basename(str_statefile_pathlist[i])).split('_')[0]
        file3 = (os.path.basename(str_alignment_pathlist[i])).split('.')[0]
        if file1 != file2 or file1 != file3:
            print(file1,file2,file3)
            print('Error: Gene files have not been correctly processed')
            exit()
        else:
            #print(file1,file2,file3)
            continue
      
    #iterate through genes
    for i in range(0, len(str_treefile_pathlist)):
        
                   #access the treefile, statefile and protein alignment
                   treefile = str_treefile_pathlist[i]
                   statefile = str_statefile_pathlist[i]
                   alignment_file = str_alignment_pathlist[i]
                   
                   gene = (os.path.basename(treefile)).split('.')[0]
                   print(f"Analysing {gene}")
                   
                   #compile a list of all the child --> parent relationships in the tree
                   child_parent_list, root = get_child_parent_relationships(treefile, outgroup_mrca, mrca_of_interest)
                   
                   #get sequences for extant sequences
                   alignment_dict = get_sequences(alignment_file)
                   alignment_dict_len = len(alignment_dict)
                   alignment_len = len(list(alignment_dict.values())[1])
                   print(f"Number of extant species found: {alignment_dict_len}")
                   print(f"Length of the {gene} alignment: {alignment_len}")
                   
                   
                   #get the ancestral state reconstruction sequences for each internal node
                   nuc_node_seq_dict = get_sequences(statefile)
                   nuc_node_len = len(nuc_node_seq_dict)
                   print(f"Number of internal node dequences found: {nuc_node_len}")
                   
                   #merge alignment dictionary and state dictionary (have all sequences in one place)
                   seq_dict = {**nuc_node_seq_dict, **alignment_dict}
                   #print(len(seq_dict))
                   
                   
                   #else you need to translate the nucleotide sequence into a protein sequence
                   sequence_dict = {}
                   for node, seq in seq_dict.items():
                       prot_seq = translate(seq)
                       if len(prot_seq) != len(seq)/3: #hard coded - need changing
                           print('error')
                       sequence_dict[node] = prot_seq

                   alignment_len = len(list(sequence_dict.values())[1])
                
                   print('Parent residues are represented by the rows and child residues are represented by the columns')
                   #initalise accumulative mutation matrix
                   residues = ['C', 'S', 'T', 'A', 'G', 'P', 'D', 'E', 'Q', 'N', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'F', 'Y', 'W']

                   
                   gene_mutation_dict = {}
                   matrix_sum_list = []
                   #count mutations for all residues in the protein
                   for res in tqdm(range(0,alignment_len),total=alignment_len):
                       r = res+1 #for the gene_mutation_dict the keys need to match to human indexing. Therefore an r of 1 will be the first position in the alignmnent
                       residue_mut_matrix = count_mutations(child_parent_list, sequence_dict, res)
                       matrix_for_sum = residue_mut_matrix.to_numpy()
                       matrix_sum = np.sum(matrix_for_sum)
                       matrix_sum_list.append(matrix_sum)

                       residue_mutation_list = residue_mut_matrix.values.tolist()
                       gene_mutation_dict[r] = residue_mutation_list
                   
                   #print(matrix_sum_list)
                   #save residue specific mutation matrices
                   residue_matrices_filepath = output_file_path + gene + '_individual_residue_matrices.txt'
                   with open(residue_matrices_filepath, 'w') as f_out:
                       for key, value in gene_mutation_dict.items():
                           f_out.write(f">{key}\n{value}\n")
                   
                   mutation_list = []
                   for site, matrix_value in gene_mutation_dict.items():
                       matrix = ast.literal_eval(str(matrix_value))
                       for i, row in enumerate(matrix):
                           for j, value in enumerate(row):
                               if i != j and value != 0:  # Ignore the diagonal (i == j) and zero entries
                                   parent = residues[i]
                                   child = residues[j]
                                   mutation_list.append([gene, site, parent, child, value])
                   mutation_df = pd.DataFrame(mutation_list, columns=['Gene', 'Aln_site', 'Parent', 'Child', 'Recurrence'])
                   output_file2_path = output_file_path + gene + '_mutation_list.cav'
                   mutation_df.to_csv(output_file2_path)

                           

                       
                
