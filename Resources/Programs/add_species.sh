#!/bin/bash

# Loop through all files in the current directory
for table_file in *_table_*.tsv; do
    # Check if the file is a table file
    if [[ -f "$table_file" ]]; then
        echo "Processing file: $table_file"
        # Get the species from the filename
        species=$(echo "$table_file" | awk -F'_' '{print $NF}' | cut -d '.' -f 1)
        # Create a temporary file to store modified lines
        temp_file=$(mktemp)
        # Read the table file line by line
        while IFS= read -r line; do
            # Extract the columns using whitespace delimiter
            read -r col1 col2 col3 <<< "$line"
            # Modify the third column to include the species
            col3="${species}_${col3}"
            # Print the modified line to the temporary file
            echo "$col1 $col2 $col3" >> "$temp_file"
        done < "$table_file"
        # Overwrite the original table file with the modified lines
        mv "$temp_file" "$table_file"
        echo "File $table_file processed successfully."
    fi
done
