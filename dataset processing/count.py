import os
import sys
import re
from collections import defaultdict
import openpyxl

# Get the current folder path
current_folder = './'

# Check if the command-line arguments are correct
if len(sys.argv) < 2:
    print("Usage: python3 111.py <number>")
    sys.exit(1)

# Get the number from the command-line argument
target_number = sys.argv[1]

# Find folders ending with the target number
folders_ending_with_target = [folder for folder in os.listdir(current_folder) if folder.endswith(target_number)]

# Initialize the dictionary for the total letter sequence occurrences
total_letter_sequence_count = defaultdict(int)

# Initialize the dictionary for the letter sequence occurrences in each folder
folder_letter_sequence_counts = {}

# Iterate through folders ending with the target number
for folder in folders_ending_with_target:
    folder_path = os.path.join(current_folder, folder)
    if os.path.isdir(folder_path):
        # Initialize the dictionary for the letter sequence occurrences in the current folder
        letter_sequence_count = defaultdict(int)
        for filename in os.listdir(folder_path):
            # Find the content after the first underscore in the filename
            after_first_underscore = filename.split('_', 1)[-1]
            # Use regular expressions to find all letter sequences before numbers in this part
            matches = re.findall(r'([a-zA-Z]+)(?=\d)', after_first_underscore)
            # Combine the found letter sequences into one for statistics
            combined_sequence = ''.join(matches)
            # Count the occurrences of each combined letter sequence
            letter_sequence_count[combined_sequence] += 1
            # Update the total letter sequence occurrence statistics
            total_letter_sequence_count[combined_sequence] += 1

        # Save the statistics of the current folder in the dictionary
        folder_letter_sequence_counts[folder] = letter_sequence_count

# Create a new Excel workbook object
wb = openpyxl.Workbook()

# Create a sheet for the total statistics
total_sheet = wb.active
total_sheet.title = 'total'

# Write the header for the total statistics
total_sheet.append(['series', 'times'])
# Write the total statistics
for sequence, count in total_letter_sequence_count.items():
    total_sheet.append([sequence, count])

# Create sheets for the statistics of each folder
for folder, letter_sequence_count in folder_letter_sequence_counts.items():
    folder_sheet = wb.create_sheet(title=folder)
    folder_sheet.append(['series', 'times'])
    for sequence, count in letter_sequence_count.items():
        folder_sheet.append([sequence, count])

# Generate the filename for the output file
output_file = f'output_combined_{target_number}.xlsx'

# Save the Excel file
wb.save(output_file)

print(f"The statistics have been saved to the file '{output_file}'.")    