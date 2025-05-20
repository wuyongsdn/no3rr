import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.build import molecule
from ase.optimize import LBFGS, BFGS, BFGSLineSearch, LBFGSLineSearch, GPMin, MDMin, FIRE
from ase.visualize.plot import plot_atoms
from ase.constraints import FixAtoms
import logging

from fairchem.core.common.relaxation.ase_utils import OCPCalculator
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import torch

import warnings
warnings.filterwarnings("ignore")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    filename='adsorbate_placement.log'
)

# Define molecule creation functions
def create_h2o(o_position, bond_length_o_h=0.957, angle_hoh=104.5):
    """Create a water molecule at the specified oxygen position."""
    adsorbate = Atoms('H2O')
    adsorbate.positions[0] = o_position  # Oxygen atom position
    # Hydrogen 1
    adsorbate.positions[1] = o_position + np.array([
        bond_length_o_h * np.cos(np.radians(angle_hoh / 2)),
        bond_length_o_h * np.sin(np.radians(angle_hoh / 2)),
        0
    ])
    # Hydrogen 2
    adsorbate.positions[2] = o_position + np.array([
        bond_length_o_h * np.cos(np.radians(-angle_hoh / 2)),
        bond_length_o_h * np.sin(np.radians(-angle_hoh / 2)),
        0
    ])
    return adsorbate

def create_nh3(n_position, bond_length_n_h=1.01, angle_hnh=107):
    """Create an ammonia molecule at the specified nitrogen position."""
    adsorbate = Atoms('NH3')
    adsorbate.positions[0] = n_position  # Nitrogen atom position
    # Hydrogen 1
    adsorbate.positions[1] = n_position + np.array([
        bond_length_n_h * np.cos(np.radians(angle_hnh - 90)),
        bond_length_n_h * np.sin(np.radians(angle_hnh - 90)),
        0
    ])
    # Hydrogen 2
    adsorbate.positions[2] = n_position + np.array([
        bond_length_n_h * np.cos(np.radians(angle_hnh - 90 + 120)),
        bond_length_n_h * np.sin(np.radians(angle_hnh - 90 + 120)),
        0
    ])
    # Hydrogen 3
    adsorbate.positions[3] = n_position + np.array([
        bond_length_n_h * np.cos(np.radians(angle_hnh - 90 + 240)),
        bond_length_n_h * np.sin(np.radians(angle_hnh - 90 + 240)),
        0
    ])
    return adsorbate

def create_h(position):
    """Create a hydrogen atom at the specified position."""
    adsorbate = Atoms('H')
    adsorbate.positions[0] = position
    return adsorbate

def create_no3(n_position, bond_length_n_o=1.22, angle_ono=120):
    """Create a nitrate ion at the specified nitrogen position."""
    adsorbate = Atoms('NO3')
    adsorbate.positions[0] = n_position  # Nitrogen atom position
    # Oxygen 1
    adsorbate.positions[1] = n_position + np.array([
        bond_length_n_o,
        0,
        0
    ])
    # Oxygen 2
    adsorbate.positions[2] = n_position + np.array([
        bond_length_n_o * np.cos(np.radians(120)),
        bond_length_n_o * np.sin(np.radians(120)),
        0
    ])
    # Oxygen 3
    adsorbate.positions[3] = n_position + np.array([
        bond_length_n_o * np.cos(np.radians(240)),
        bond_length_n_o * np.sin(np.radians(240)),
        0
    ])
    return adsorbate

# Create adsorbate function mapping
adsorbate_functions = {
    'nh3': create_nh3,
    'no3': create_no3,
    'h': create_h,
    'h2o': create_h2o
}

# Input and output path dictionary, keys are system names (customizable), 
# values are (input folder path, output folder path)
input_dir = r'/home/sdengning/method/fairchem/predict/first/cal'
input_output_dirs = {
    'nh3': (input_dir, r'/home/sdengning/method/fairchem/predict/first/nh3'),
    'no3': (input_dir, r'/home/sdengning/method/fairchem/predict/first/no3'),
    'h': (input_dir, r'/home/sdengning/method/fairchem/predict/first/h'),
    'h2o': (input_dir, r'/home/sdengning/method/fairchem/predict/first/h2o'),
}

# Setup OCPCalculator
calc = OCPCalculator(
        model_name="DimeNet++-IS2RE-OC20-All",
        local_cache="/home/sdengning/method/fairchem/predict/first",
        cpu=False
    )

# Process each adsorbate's corresponding folder
for adsorbate_name, (input_dir, output_dir) in input_output_dirs.items():
    try:
        # Check if the output folder exists, create it if not
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            logging.info(f"Created directory: {output_dir}")
            
        # Get all CIF files
        cif_files = [f for f in os.listdir(input_dir) if f.endswith('.cif')]
        logging.info(f"Found {len(cif_files)} CIF files in {input_dir}")

        # Process each CIF file
        for cif_file in tqdm(cif_files, desc=f"Processing {adsorbate_name} CIF files"):
            try:
                slab = read(os.path.join(input_dir, cif_file))  # Read structure from CIF file
                logging.info(f"Successfully read file: {cif_file}")

                # Get lattice parameters
                lattice_param_x = slab.get_cell()[0, 0]
                lattice_param_y = slab.get_cell()[1, 1]
                lattice_param_z = slab.get_cell()[2, 2]  # Get z-direction lattice parameter

                # Calculate half of the lattice height
                surface_x_max = slab.get_positions()[:, 0].max()
                surface_x_min = slab.get_positions()[:, 0].min()
                surface_y_max = slab.get_positions()[:, 1].max()
                surface_y_min = slab.get_positions()[:, 1].min()
                half_z_height = (slab.get_positions()[:, 2].max() + slab.get_positions()[:, 2].min()) / 2

                # Find atoms below half_z_height and fix them
                fixed_indices = [i for i, pos in enumerate(slab.get_positions()) if pos[2] < half_z_height]
                constraint = FixAtoms(indices=fixed_indices)
                slab.set_constraint(constraint)

                # Set the height of the adsorbate above the surface
                height_above_atoms = 0.05 * lattice_param_z  # Set height as a fraction of lattice parameter z

                # Generate grid positions on the surface
                num_points = 100
                grid_x, grid_y = np.meshgrid(
                    np.linspace(surface_x_min, surface_x_max, int(np.sqrt(num_points))),
                    np.linspace(surface_y_min, surface_y_max, int(np.sqrt(num_points)))
                )
                grid_positions = np.column_stack([grid_x.flatten(), grid_y.flatten()])

                # Loop over each grid position to place the adsorbate molecule
                for i, (x, y) in enumerate(grid_positions):
                    try:
                        # Copy the slab
                        slab_with_adsorbate_up = slab.copy()

                        # Determine the position above the surface where the adsorbate will be placed
                        anchor_position = np.array([x, y, slab.get_positions()[:, 2].max() + height_above_atoms])

                        # Create the adsorbate at the specified position
                        adsorbate = adsorbate_functions[adsorbate_name](anchor_position)

                        # Add the adsorbate to the slab
                        slab_with_adsorbate_up += adsorbate
                        slab_with_adsorbate_up.calc = calc
                        
                        # Calculate potential energy
                        final_energy_up = slab_with_adsorbate_up.get_potential_energy()

                        # Save the structure without optimization
                        up_file_path = rf'{output_dir}/{os.path.splitext(cif_file)[0]}_{i}_E_{final_energy_up:.6f}.cif'
                        write(up_file_path, slab_with_adsorbate_up)
                        
                        logging.info(f"Successfully saved structure: {up_file_path}")
                        
                    except Exception as e:
                        logging.error(f"Error processing grid position {i}: {str(e)}")
                        continue
                        
            except Exception as e:
                logging.error(f"Error processing CIF file {cif_file}: {str(e)}")
                continue
                
    except Exception as e:
        logging.error(f"Error processing adsorbate {adsorbate_name}: {str(e)}")
        continue

logging.info("All adsorbate placement calculations completed")    
