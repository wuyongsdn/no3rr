import os
import pandas as pd
import warnings
from pymatgen.core import Structure
from tqdm import tqdm

warnings.filterwarnings("ignore", category=UserWarning, module='pymatgen')

def generate_custom_format(structure):
    lattice = structure.lattice
    atoms = structure.sites

    custom_format = []

    # basic information
    custom_format.append("data_C")
    custom_format.append("_symmetry_space_group_name_H-M   'P 1'")
    custom_format.append(f"_cell_length_a   {lattice.a:.8f}")
    custom_format.append(f"_cell_length_b   {lattice.b:.8f}")
    custom_format.append(f"_cell_length_c   {lattice.c:.8f}")
    custom_format.append(f"_cell_angle_alpha   {lattice.alpha:.8f}")
    custom_format.append(f"_cell_angle_beta   {lattice.beta:.8f}")
    custom_format.append(f"_cell_angle_gamma   {lattice.gamma:.8f}")
    custom_format.append("_symmetry_Int_Tables_number   1")
    custom_format.append(f"_chemical_formula_structural   {structure.formula}")
    custom_format.append(f"_chemical_formula_sum   {structure.composition.formula.replace(' ', '')}")
    custom_format.append(f"_cell_volume   {lattice.volume:.8f}")
    custom_format.append(f"_cell_formula_units_Z   {structure.composition.num_atoms}")

    # symmetry information
    custom_format.append("loop_")
    custom_format.append(" _symmetry_equiv_pos_site_id")
    custom_format.append(" _symmetry_equiv_pos_as_xyz")
    custom_format.append("  1  'x, y, z'")

    # atomic information
    custom_format.append("loop_")
    custom_format.append(" _atom_site_type_symbol")
    custom_format.append(" _atom_site_label")
    custom_format.append(" _atom_site_symmetry_multiplicity")
    custom_format.append(" _atom_site_fract_x")
    custom_format.append(" _atom_site_fract_y")
    custom_format.append(" _atom_site_fract_z")
    custom_format.append(" _atom_site_occupancy")

    for i, atom in enumerate(atoms):
        custom_format.append(f"  {atom.specie}  {atom.specie}{i}  1  {atom.frac_coords[0]:.8f}  {atom.frac_coords[1]:.8f}  {atom.frac_coords[2]:.8f}  1")

    return "\n".join(custom_format)

def process_files_in_folder(folder):
    data = []
    files = [f for f in os.listdir(folder) if f.endswith('.cif') or f.endswith('.vasp')]
    error_files = []

    for file_name in tqdm(files, desc="Processing files"):
        file_path = os.path.join(folder, file_name)
        try:
            structure = Structure.from_file(file_path)
            custom_format = generate_custom_format(structure)

            file_base_name = os.path.splitext(file_name)[0]
            data.append([file_base_name, custom_format])
        except Exception as e:
            error_files.append(file_name)
            print(f"Error processing {file_name}: {e}")

    if error_files:
        print(f"\nThe following files could not be processed: {', '.join(error_files)}")

    return data

def save_to_xlsx(data, output_file):
    df = pd.DataFrame(data, columns=['FileName', 'Content'])
    df.to_excel(output_file, index=False)

def main():
    folder = r'E:\oc20\oc\no3\all'  # input dir
    output_file = r'E:\oc20\oc\no3\all\output111.xlsx' # output

    data = process_files_in_folder(folder)
    save_to_xlsx(data, output_file)

if __name__ == "__main__":
    main()
