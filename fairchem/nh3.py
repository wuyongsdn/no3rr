import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.build import molecule
from ase.optimize import LBFGS, BFGS, BFGSLineSearch, LBFGSLineSearch, GPMin, MDMin, FIRE
from ase.visualize.plot import plot_atoms
from ase.constraints import FixAtoms

from fairchem.core.common.relaxation.ase_utils import OCPCalculator
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import sys  # Added for command-line arguments

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")


def create_nh3o(n_position, bond_length_n_h=1.01, bond_length_n_o=1.40, angle_hno=107, angle_hnh=107):
    """Create NH3O adsorbate molecule."""
    adsorbate = Atoms('NH3O')
    adsorbate.positions[0] = n_position
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
    # Oxygen
    adsorbate.positions[4] = n_position + np.array([
        bond_length_n_o * np.cos(np.radians(angle_hno)),
        bond_length_n_o * np.sin(np.radians(angle_hno)),
        0
    ])
    return adsorbate


def create_nh2(n_position, bond_length_n_h=1.01, angle_hnh=104.5):
    """Create NH2 adsorbate molecule."""
    adsorbate = Atoms('NH2')
    adsorbate.positions[0] = n_position
    # Hydrogen 1
    adsorbate.positions[1] = n_position + np.array([
        bond_length_n_h * np.cos(np.radians(angle_hnh / 2)),
        bond_length_n_h * np.sin(np.radians(angle_hnh / 2)),
        0
    ])
    # Hydrogen 2
    adsorbate.positions[2] = n_position + np.array([
        bond_length_n_h * np.cos(np.radians(-angle_hnh / 2)),
        bond_length_n_h * np.sin(np.radians(-angle_hnh / 2)),
        0
    ])
    return adsorbate


# Set up OCPCalculator
calc = OCPCalculator(
    trainer="energy",
    model_name="EquiformerV2-153M-S2EF-OC20-All+MD",
    local_cache="",# model path
    cpu=False,
)

# Dictionary of adsorbate creation functions
adsorbate_functions = {
    'nh3o': lambda n_position: create_nh3o(n_position),
    'nh2': lambda n_position: create_nh2(n_position),
    # ... (all other entries as before) ...
    'nho2': lambda n_position: create_nho2(n_position)
}

# Dictionary of input/output directories: key=adsorbate name, value=(input_dir, output_dir)
input_output_dirs = {
    'nh3': (input_dir, output_dir),
    'nh3o': (input_dir, output_dir),
    # ... (all other entries as before) ...
    'nho2': (input_dir, output_dir)
}

def main():
    # Get adsorbate to process from command-line argument (e.g., 'no3'), or process all if none specified
    target_adsorbate = sys.argv[1].lower() if len(sys.argv) > 1 else None
    
    for adsorbate_name, (input_dir, output_dir) in input_output_dirs.items():
        # Skip if target adsorbate is specified and doesn't match current name
        if target_adsorbate is not None and adsorbate_name != target_adsorbate:
            continue
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Get all CIF files in input directory
        cif_files = [f for f in os.listdir(input_dir) if f.endswith('.cif')]
        
        for cif_file in tqdm(cif_files, desc=f"Processing {adsorbate_name} CIF files"):
            slab = read(os.path.join(input_dir, cif_file))  # Read slab structure from CIF
            
            # Get lattice parameters (assuming 3D cell, focusing on z-direction for height)
            lattice_param_z = slab.get_cell()[2, 2]
            surface_positions = slab.get_positions()
            surface_z_max = surface_positions[:, 2].max()
            surface_z_min = surface_positions[:, 2].min()
            half_z_height = (surface_z_max + surface_z_min) / 2
            
            # Fix atoms below half_z_height
            fixed_indices = [i for i, pos in enumerate(surface_positions) if pos[2] < half_z_height]
            slab.set_constraint(FixAtoms(indices=fixed_indices))
            
            # Define adsorption height above surface (5% of lattice z parameter)
            height_above_surface = 0.05 * lattice_param_z
            n_grid_points = 100
            surface_x_min = surface_positions[:, 0].min()
            surface_x_max = surface_positions[:, 0].max()
            surface_y_min = surface_positions[:, 1].min()
            surface_y_max = surface_positions[:, 1].max()
            
            # Generate grid positions for adsorption
            grid_x, grid_y = np.meshgrid(
                np.linspace(surface_x_min, surface_x_max, int(np.sqrt(n_grid_points))),
                np.linspace(surface_y_min, surface_y_max, int(np.sqrt(n_grid_points)))
            )
            grid_positions = np.column_stack([grid_x.flatten(), grid_y.flatten()])
            
            for idx, (x, y) in enumerate(grid_positions):
                slab_with_adsorbate = slab.copy()
                n_position = np.array([x, y, surface_z_max + height_above_surface])
                
                # Create adsorbate molecule
                adsorbate = adsorbate_functions[adsorbate_name](n_position)
                slab_with_adsorbate += adsorbate
                slab_with_adsorbate.calc = calc
                
                # Structure optimization
                optimizer = LBFGS(slab_with_adsorbate, logfile=None)
                optimizer.run(fmax=0.05, steps=100)
                
                # Get final energy and save structure
                final_energy = slab_with_adsorbate.get_potential_energy()
                output_filename = f"{output_dir}/{os.path.splitext(cif_file)[0]}_{idx}_E_{final_energy:.6f}.cif"
                write(output_filename, slab_with_adsorbate)

if __name__ == "__main__":
    main()