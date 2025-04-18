#!/bin/bash

name_start=$(pwd | awk -F / '{print $NF}')

echo "$name_start"

echo "count of line:"
read count

mkdir ${name_start}_1 ${name_start}_2 ${name_start}_3

# Initialize a counter for numbering files
counter=1

# Loop through each .extxyz file
for file in *.extxyz; do
    # Generate merged content
    merged_content=$(python3 series.py "$file" $count)
    new_name=$(python3 series.py "$file" $count | tr -d '\n' | xargs -n1)

    # Formulate the new file name with a sequential number
    new_file="${name_start}_${new_name}_${counter}_${file}"
    
    # Copy the file with the new name
    cp "$file" "$new_file"
    
    # Count the number of lines in merged content
    num_lines=$(echo "$merged_content" | wc -l)
    
    # Determine the destination folder based on the number of lines
    if [ "$num_lines" -eq 1 ]; then
        destination="${name_start}_1"
    elif [ "$num_lines" -eq 2 ]; then
        destination="${name_start}_2"
    elif [ "$num_lines" -eq 3 ]; then
        destination="${name_start}_3"
    else
        echo "Unexpected number of lines: $num_lines"
        continue
    fi
    
    # Move the original file to the appropriate folder
    mv "$new_file" "$destination/"
    
    # Increment the counter for the next file
    (( counter++ ))
done
