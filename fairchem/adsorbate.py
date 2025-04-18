import os
import shutil
from ase.io import read
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

def move_cif_if_needed(cif_path, target_dir, element_to_find, check_radius, count_radius, required_element=None, required_count=None):
    try:
        # Read the CIF file
        atoms = read(cif_path)
        
        # Get the scaled positions
        scaled_positions = atoms.get_scaled_positions()  # Get the scaled positions
        
        # Find all atoms of the specified element
        selected_atoms = [atom for atom in atoms if atom.symbol == element_to_find]
        
        if not selected_atoms:
            return
        
        # Find the atom with the highest z-coordinate
        highest_z_atom = max(selected_atoms, key=lambda atom: scaled_positions[atom.index][2])
        #highest_z_atom = min(selected_atoms, key=lambda atom: scaled_positions[atom.index][2])
        
        # Step 1: Check if there are enough atoms of the required_element type within the check_radius centered at this atom
        atoms_in_check_sphere = [atom for atom in atoms if calculate_distance_with_periodicity(highest_z_atom, atom, atoms, check_radius)]
        selected_atoms_in_check_sphere = [atom for atom in atoms_in_check_sphere if atom.symbol == required_element]
        
        if len(selected_atoms_in_check_sphere) < required_count:
            move_file(cif_path, target_dir)
            return
        
        # Step 2: Check if there are other atoms within the count_radius centered at this atom, excluding the element_to_find and required_element types
        atoms_in_count_sphere = [atom for atom in atoms if calculate_distance_with_periodicity(highest_z_atom, atom, atoms, count_radius)]
        atoms_in_count_sphere_filtered = [atom for atom in atoms_in_count_sphere if atom.symbol not in [element_to_find, required_element]]
        
        if not atoms_in_count_sphere_filtered:
            move_file(cif_path, target_dir)
            return

        # If both conditions are met, do not move the file
        return
    
    except Exception as e:
        # If an error occurs while reading or processing the CIF file, print the error message and record the file
        print(f"Error processing {cif_path}: {e}")
        with open("error_files.log", "a") as log_file:
            log_file.write(f"{cif_path}: {e}\n")

def calculate_distance_with_periodicity(atom1, atom2, atoms, radius):
    # Use the `ase` method to calculate the periodic distance between two atoms
    dist = atoms.get_distance(atom1.index, atom2.index, mic=True)  # mic=True indicates considering periodic boundary conditions
    return dist <= radius

def move_file(cif_path, target_dir):
    # Move the file to the target directory
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    shutil.move(cif_path, os.path.join(target_dir, os.path.basename(cif_path)))

def process_cifs_in_directory(source_dir, target_dir, element_to_find, check_radius, count_radius, required_element=None, required_count=None):
    # Get all CIF files
    cif_files = [os.path.join(source_dir, f) for f in os.listdir(source_dir) if f.endswith('.cif')]
    
    # Use ThreadPoolExecutor to process CIF files in parallel
    with ThreadPoolExecutor(max_workers=12) as executor:
        # Use the tqdm package to add a progress bar to the file list
        list(tqdm(executor.map(lambda cif: move_cif_if_needed(cif, target_dir, element_to_find, check_radius, count_radius, required_element, required_count), cif_files), total=len(cif_files), desc="Processing CIF files"))

# Example usage
source_directory = r"/home/sdengning/fairchem/test/h2o_up/"
target_directory = r"/home/sdengning/fairchem/test/h2o_up/none"

'''
count_radius = 4  # Radius for the second step (unit: Å)
required_count = 3  # Required number of the specified element
element_to_find = "N"  # Specified element
check_radius = 1.2  # Radius for the first step (unit: Å)
required_element = "H"  # Element required within the sphere
'''

'''
count_radius = 4  # Radius for the second step (unit: Å)
required_count = 3  # Required number of the specified element
element_to_find = "N"  # Specified element
check_radius = 1.7  # Radius for the first step (unit: Å)
required_element = "O"  # Element required within the sphere
'''


count_radius = 4  # Radius for the second step (unit: Å)
required_count = 2  # Required number of the specified element
element_to_find = "O"  # Specified element
check_radius = 1.1  # Radius for the first step (unit: Å)
required_element = "H"  # Element required within the sphere


process_cifs_in_directory(source_directory, target_directory, element_to_find, check_radius, count_radius, required_element, required_count)