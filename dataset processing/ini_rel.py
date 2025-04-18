import os
import shutil
import re
import concurrent.futures

# Use the current folder as the input directory
input_directory = os.getcwd()
initial_output_directory = os.path.join(input_directory, 'initial')
relaxed_output_directory = os.path.join(input_directory, 'relaxed')

def process_file(file_path, initial_output_directory, relaxed_output_directory):
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    # Find lines that start with a number
    digit_lines = [i for i, line in enumerate(lines) if re.match(r'^\d', line)]

    if len(digit_lines) < 2:
        print(f"The file {file_path} has less than two lines starting with a number, cannot process.")
        return

    # Get the file name and generate new file paths
    base_filename = os.path.basename(file_path)
    initial_filename = base_filename.replace('.extxyz', '.extxyz')
    relaxed_filename = base_filename.replace('.extxyz', '.extxyz')
    initial_output_path = os.path.join(initial_output_directory, initial_filename)
    relaxed_output_path = os.path.join(relaxed_output_directory, relaxed_filename)

    # Get the first block of content starting from the first line with a number until the second line with a number
    initial_content = lines[digit_lines[0]: digit_lines[1]]
    with open(initial_output_path, 'w', encoding='utf-8') as initial_file:
        initial_file.writelines(initial_content)

    # Get the content from the last line starting with a number to the end of the file
    relaxed_content = lines[digit_lines[-1]:]
    with open(relaxed_output_path, 'w', encoding='utf-8') as relaxed_file:
        relaxed_file.writelines(relaxed_content)

def process_directory(input_directory, initial_output_directory, relaxed_output_directory):
    with concurrent.futures.ProcessPoolExecutor(max_workers=52) as executor:
        futures = []
        for root, _, files in os.walk(input_directory):
            for file in files:
                if file.endswith('.extxyz'):
                    file_path = os.path.join(root, file)
                    futures.append(executor.submit(process_file, file_path, initial_output_directory, relaxed_output_directory))
        
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"An error occurred while processing the file: {e}")

# Create output directories
os.makedirs(initial_output_directory, exist_ok=True)
os.makedirs(relaxed_output_directory, exist_ok=True)

# Process all files in the input directory
process_directory(input_directory, initial_output_directory, relaxed_output_directory)