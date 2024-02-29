#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 10:57:17 2024

@author: alexaburchak
"""
import os
import re 

# Directory containing the TSV files
directory = "/home/inf-40-2023/malaria/Complete_Only_Tables"

# Output file to save common genes with their associated sequence IDs
output_file = "common_genes_with_species.txt"

# Dictionary to store gene counts
common_genes = {}

# Dictionary to store gene IDs and their corresponding sequence IDs
gene_sequence_map = {}

# Loop through each file in the specified directory
for filename in os.listdir(directory):
    # Check if the file is a TSV file
    if filename.endswith(".tsv"):
        # Construct the full file path
        file_path = os.path.join(directory, filename)
        
        # Open each TSV file and read line by line
        with open(file_path, 'r') as file:
            # Iterate through each line in the file
            for line_number, line in enumerate(file, start=1):
                # Split the line into columns
                columns = line.strip().split(' ')
    
                # Extract gene ID and sequence ID from the columns
                gene_id = columns[0]
                sequence_id = columns[2]
                    
                # Store gene ID and sequence ID mapping
                gene_info = (gene_id, sequence_id)
                common_genes[gene_id] = common_genes.get(gene_id, set())
                common_genes[gene_id].add(sequence_id)

# Identify complete genes present in all 8 files
complete_genes = [gene for gene, sequences in common_genes.items() if len(sequences) == 8]

# Write common genes with their associated sequence IDs to the output file
with open(output_file, 'w') as out_file:
    for gene in complete_genes:
        # Write gene ID and associated sequence IDs to the output file, separated by tabs
        out_file.write(gene + "\t" + "\t".join(common_genes[gene]) + "\n")

# Print message indicating completion and the location of the output file
print("Common genes with sequence IDs saved to", output_file)
