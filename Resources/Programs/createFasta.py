#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

# Dictionary to store sequence IDs for each busco gene
busco_sequences = {}

# Read the file with the new format
with open('common_genes_with_species.txt', 'r') as file:
    for line in file:
        parts = line.strip().split()
        busco_id = parts[0]
        # Initialize a dictionary to store species and sequence IDs for the current busco ID
        species_sequence_ids = {}
        for entry in parts[1:]:
            species, sequence_id = entry.split('_', 1)
            species_sequence_ids[species] = sequence_id
        busco_sequences[busco_id] = species_sequence_ids

# Function to extract sequence from sequence ID file
def extract_sequence(sequence_id, species):
    # Construct the filename based on the species
    filename = f'new_{species}.faa'
    
    # Check if the file exists
    if os.path.exists(filename):
        # Open the sequence file
        with open(filename, 'r') as seq_file:
            sequence = ''  # Variable to store the sequence
            found_sequence = False  # Flag to indicate if the sequence is found
            
            # Iterate over each line in the sequence file
            for line in seq_file:
                # Check if the line is a header line starting with '>'
                if line.startswith('>'):
                    # If the sequence has already been found, return it
                    if found_sequence:
                        return sequence
                    # Check if the sequence ID is found in the header line
                    if sequence_id in line:
                        found_sequence = True  # Set the flag to True as the sequence is found
                    continue  # Continue to the next line if the header line doesn't contain the sequence ID
                # If the sequence is found and the line is not a header line, append it to the sequence
                if found_sequence:
                    sequence += line.strip()
            
            # If the sequence is found, return it
            if found_sequence:
                return sequence
            # If the sequence is not found, return None
            else:
                return None
    else:
        # If the file doesn't exist, print a message and return None
        print(f"File {filename} not found for species {species}")
        return None

# Create fasta files for each busco gene ID
output_directory = 'Fasta_files'
os.makedirs(output_directory, exist_ok=True)

# Loop through each busco gene ID and its associated sequence IDs
for busco_id, species_sequence_ids in busco_sequences.items():
    # Define the filename for the fasta file corresponding to the current busco gene ID
    fasta_filename = os.path.join(output_directory, f'{busco_id}.fasta')
    
    # Open the fasta file for writing
    with open(fasta_filename, 'w') as fasta_file:
        # Iterate over each species and sequence ID for the current busco gene ID
        for species, sequence_id in species_sequence_ids.items():
            # Extract the sequence for the current sequence ID and species
            sequence = extract_sequence(sequence_id, species)
            
            # If sequence is found for the current species, write it to the fasta file
            if sequence is not None:
                fasta_file.write(f'>{species}\n{sequence}\n')
            else:
                print(f"Skipping sequence ID {sequence_id} for busco ID {busco_id}. Sequence not found in species {species}")

print("Fasta files created successfully.")
