import os
import re
import csv
from joblib import Parallel, delayed

def extract_free_energy_from_files(directory):
    free_energies = {}
    pattern = re.compile(r'free_energy=([-]?\d+\.\d+)')

    # Traverse all files in the given directory
    files = [filename for filename in os.listdir(directory) if filename.endswith('.extxyz')]
    
    # Use joblib Parallel for concurrent processing
    results = Parallel(n_jobs=-1)(delayed(extract_free_energy_from_file)(os.path.join(directory, filename), pattern) for filename in files)
    
    # Collect results into dictionary
    for filename, result in zip(files, results):
        if result is not None:
            # Remove the .txt extension for the prefix
            prefix = filename.rsplit('.extxyz', 1)[0]
            free_energies[prefix] = result
    
    return free_energies

def write_to_csv(data, output_file):
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['prefix', 'free_energy'])
        for prefix, free_energy in data.items():
            writer.writerow([prefix, free_energy])

# Example usage
directory = os.getcwd() 
output_file = os.path.join(directory, 'id_prop.csv')
free_energies = extract_free_energy_from_files(directory)
write_to_csv(free_energies, output_file)

print(f"Data has been written to {output_file}")
