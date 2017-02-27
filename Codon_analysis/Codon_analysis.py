# -*- coding: utf-8 -*-
"""
Created on Sun Apr 03 16:37:20 2016

Author: Justin W. Kenney
Email: jkenney9a@gmail.com

Code for analyzing codon frequency and usage for genes

Parameters:


Codon analysis code

Entrez.email=...
"""

import sys, os, glob
from Bio import Entrez
from Bio import SeqIO
import csv
import string
import numpy as np



def get_gene_list(filename):
    """
    Input: File with list of genes or list of gene output names from proteomics
    Output: List of gene names in file
    """
    
    gene_list = []
    
    if filename.split('.')[-1] == 'csv':
        with open(filename, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if len(row) > 0:
                    if "GN=" in row:
                        gene = get_gene_name(row)
                    elif ';' in row[0]:
                        for gene in row[0].split(';'):
                            gene = gene.strip()
                    else:
                        gene = row[0]
                    
                    gene_list.append(gene)
    
    elif filename.split('.')[-1] == 'txt':
        f = open(filename)
        for line in f:
            if "GN=" in line:
                gene = get_gene_name(line)
                gene_list.append(gene)
            else:
                
                gene_list.append(line.rstrip("\n"))
        
        f.close()
    
    return gene_list
    
def get_gene_name(line):
    """
    Input: A line read in from a txt or csv file from some proteomic data
    that contains a 'GN=' part before the gene name
    
    Output: The gene name pulled out of the line
    """
    gene = ""
    start = line.find("GN=")
    while line[start+3] != " ":
        gene += line[start+3]
        start += 1
    
    return gene
    

def get_CDS_fasta_files(gene, organism):
    """
    Input: gene name and organism
    Output: fasta file with CDS sequence from NCBI in fasta directory
    If the file already exists does not download from NCBI
    """

    if not os.path.exists("fasta"):
        os.mkdir("fasta")
    
    #Checks to see if the file has already been downloaded previously
    #If it has not, then it will be downloaded, otherwise it will use the 
    #previously downloaded file. This speeds up processing time and minimizes
    #the strain on NCBI resources
    filename = gene + "_" + organism
    if len(glob.glob("fasta\\" + filename + "*.fasta")) == 0:
                   
        
         #Get gene id from Gene database; exclude predicted sequences        
        search_term = gene + "[Gene Name] AND " + organism + \
        "[Organism] AND mRNA[Filter] AND RefSeq[Filter] NOT PREDICTED[Title]"
        
        search_handle = Entrez.esearch(db="nucleotide", term = search_term)
        
        #Parse the resulting xml file into a dictionary and get gene ID numbers
        #associated with specific gene records
        search_record = Entrez.read(search_handle)
        gene_ids = search_record["IdList"]
        search_handle.close()
        
        count = 1
        
        #Gets the CDS file for each gene from Entrez and creates a fasta file
        #for each gene entry
        for g in gene_ids:
            handle = Entrez.efetch(db="nucleotide", id=g, rettype="fasta_cds_na",\
            retmode = "text")
            record = SeqIO.read(handle, format="fasta")
            SeqIO.write(record, "fasta\\" + filename + "_" + str(count)\
            + ".fasta", "fasta")
            count += 1


def genefile_to_seq(gene, organism):
    """
    Input: Gene and organism name for which a file exists
    Output: List of sequence objects associated with that gene
    """
    filename = gene + "_" + organism
    filenames = glob.glob("fasta\\" + filename + "*.fasta")
    
    sequences = []
    for f in filenames:
        sequences.append(SeqIO.read(f, "fasta"))
    
    return sequences


def normalized_codon_mapping(codon_mapping):
    """
    Input: Codon frequency mapping dictionary in form: {'codon':['AA',Freq])
    Output: Normalized codon frequency map relative to max frequency for each 
    codon {'codon':'rel freq'}
    Note: AA's in codon map must be encoded as single letters
    """

    AA_List = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P',\
    'S','T','W','Y','V','*']

    #Creates a dictionary mapping AA letter to frequencies associated with AA
    AA_Freq_Dict = {x:[codon_mapping[y][1] 
    for y in codon_mapping if codon_mapping[y][0] == x] for x in AA_List}

    #Dictionary mapping amino acid to maximum frequency value in codon_mapping
    AA_Max_Freq = {x:max(AA_Freq_Dict[x]) for x in AA_Freq_Dict}
    #Deal with zeros at stop codons:
    for x in AA_Max_Freq:
        if AA_Max_Freq[x] == 0:
            AA_Max_Freq[x] = 1

    #Dictionary mapping codon to relative codon frequency
    return {x:[codon_mapping[x][0],(codon_mapping[x][1] / AA_Max_Freq[codon_mapping[x][0]])] 
    for x in codon_mapping}


def sequence_to_codons(sequence):
    """
    Input: Biopython sequence object
    Output: list containing codons in order
    """
    codons = []
    if len(sequence.seq) % 3 == 0:
         for x in range(len(sequence.seq)/3):
             codons.append(str(sequence.seq[x:x+3]))
             x += 3
    else:
        return "Error, not a coding sequence"
            
    return codons


def codon_map_gene(gene, organism, codon_map):
    """
    Input: Gene name, organism and codon map
    Output: A list of codon mappings for given gene
    """
    sequences = genefile_to_seq(gene, organism)
    gene_codons = []
    for g in sequences:
        gene_codons.append(sequence_to_codons(g))
    
    gene_codon_maps = []    
    
    for x in gene_codons:
        gene_codon_maps.append([codon_map[y][1] for y in x])
    
    return gene_codon_maps


def import_codon_map(filename):
    """
    Input: filename (csv w/ comma delimiter) containing a codon map where 
    the first column is codon, second column is amino acid abbreviation (one 
    letter) and third column is value (e.g, time to decode or relative 
    abundance etc.)
    
    Output: A codon mapping dictionary with the format of
    {codon: [AA, value]}
    
    """
    codon_map = {}
 
    with open(filename, "r") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            row[0] = string.replace(row[0], "U", "T")
            codon_map[row[0]] = [row[1], float(row[2])]    
    
    return codon_map
    


def gene_codon_analysis(gene, organism, codon_map):
    """
    Input: gene name, organism, and codon mapping
    
    This function assumes the mRNA fasta files have already been downloaded
    
    Output: A dictionary of the format: {gene_name: name, avg: mappings divided by length
    of transcript, total: sum of mappings, normalized: normalized mappings relative
    to maximum value for each AA, range_ratio: a ratio of the range of differences
    for protein lengths relative to maximum protein length, length_range:range of 
    protein lengths from downloaded transcripts}
    """
    
    gene_codon_maps = codon_map_gene(gene, organism, codon_map)
    codon_map_norm = normalized_codon_mapping(codon_map)
    gene_codon_maps_norm = codon_map_gene(gene, organism, 
                                                   codon_map_norm)
                                                  
    if len(gene_codon_maps) == 0:
        return {"gene_name": gene}
    #Calculate sums, lengths and averages for given codon mapping
    codon_sums = []
    protein_lengths = []
    codon_first_25 = []
    codon_first_third = []
    
    for mapping in gene_codon_maps:
        codon_sums.append(sum(mapping))
        protein_lengths.append(len(mapping))
        codon_first_25.append(sum(mapping[0:24]))
        codon_first_third.append(sum(mapping[0:int(len(mapping)/3)]))
    
    codon_averages = [float(sums)/lengths
                      for sums,lengths in zip(codon_sums, protein_lengths)]
    
    codon_first_25_avgs = [float(sums)/25 for sums in codon_first_25]
    
    codon_first_third_avgs = [float(sums)/int(lengths / 3) for sums, lengths in 
                        zip(codon_first_third, protein_lengths)]
    
    #Calculate normalized codon usage (i.e, relative to theoretical maximum)
    codon_sums_norm = []
    
    for mapping in gene_codon_maps_norm:
        codon_sums_norm.append(sum(mapping) - 1) #-1 to deal with stop codon
        
    codon_averages_norm = [float(sums)/lengths 
                      for sums, lengths in zip(codon_sums_norm, protein_lengths)]
    
    protein_length_min = min(protein_lengths)
    protein_length_max = max(protein_lengths)
    protein_length_median = np.median(protein_lengths)
    
    output_avg = sum(codon_averages)/len(codon_averages)
    output_total = sum(codon_sums)/len(codon_sums)
    output_norm = sum(codon_averages_norm)/len(codon_averages_norm)
    output_range_ratio = (protein_length_max - protein_length_min) / float(protein_length_max)
    output_protein_size_range = str(protein_length_min) + " - " + str(protein_length_max)
    output_first_25 = sum(codon_first_25_avgs)/len(codon_first_25_avgs)
    output_first_third = sum(codon_first_third_avgs)/(len(codon_first_third_avgs))
        
    output_dict = {"gene_name":gene, "avg":output_avg, "total":output_total,
                   "normalized":output_norm, "range_ratio":output_range_ratio,
                   "length_range":output_protein_size_range, 
                   "median_protein_length":protein_length_median, "first_25_avg":
                       output_first_25, "first_third_avg": output_first_third}
    
    return output_dict
    
        
    
    
    


if __name__ == "__main__":
    
    import sys
    
    download = False
    Entrez.email = None
    
    for arg in sys.argv[1:]:
        try:
            name, value = arg.split('=', 1)
        except: print "Error parsing command line argument. No '=' found"
        
        if name.lower() == "--gene_list" or name.lower() == "--list":
            gene_list_filename = value
        
        elif name.lower() == "--email":
            Entrez.email = value
            
        elif name.lower() == "--output":
            output_filename = value
            
        elif name.lower() == "--codon_map" or name.lower() == "--map":
            codon_map_filename = value
        
        elif name.lower() == "--download":
            download = value
            
        elif name.lower() == "--organism":
            organism = value
        
    gene_list = get_gene_list(gene_list_filename)
    codon_map = import_codon_map(codon_map_filename)
    
    if download == True or download == "T":
        if Entrez.email == None:
            print "Error. If downloading from NCBI need to provide \
            email address using --email=<your email address here>"
            sys.exit()
        
        for gene in gene_list:
            get_CDS_fasta_files(gene, organism)
    
    
    with open(output_filename, 'wb') as csvfile:
        fieldnames = ['gene_name', 'avg', 'total', 'normalized', 'range_ratio',
                      'length_range', 'median_protein_length', 'first_25_avg', 
                      'first_third_avg']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, restval="ERROR")
        writer.writeheader()
        
        counter=0
        for gene in gene_list:
            writer.writerow(gene_codon_analysis(gene, organism, codon_map))
            counter += 1
            print counter, "of", len(gene_list),"genes analyzed!    \r",
        
        