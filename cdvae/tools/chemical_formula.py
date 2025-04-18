import os
from pymatgen.io.cif import CifParser

import warnings
warnings.filterwarnings("ignore")

# Specify the folder path where the CIF files are located
folder_path = r'E:\oc20\oc\all\cif\initial'  
output_path = r'E:\oc20\oc\all\cif\initial\111.txt'
# Open the output txt file for writing data
with open(output_path, 'w') as output_file:
    # Iterate through all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.cif'):
            file_path = os.path.join(folder_path, filename)
            # Parse the CIF file
            parser = CifParser(file_path)
            structure = parser.parse_structures(primitive=False)[0]  # Get the first structure
            
            # Extract element symbols and quantities (ensure the quantities are integers)
            element_counts = structure.composition.items()
            elements_str = ''.join([f"{element}{int(count)}" for element, count in element_counts])

            # Write the file name and the corresponding element symbols and quantities, separated by a tab
            output_file.write(f"{filename}\t{elements_str}\n")

print("Data has been successfully written to the 111.txt file!")