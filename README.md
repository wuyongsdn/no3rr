# Materials Science Automation Toolkit

This repository contains Python scripts and configuration files for materials science applications, including dataset processing, crystal structure analysis, and energy calculations.

## 📁 Script Descriptions

### Core Processing Scripts
- **file.sh**  
  Collaborates with `series.py` to organize files into directories based on substrate type (ternary/quaternary alloys)

- **count.py**  
  Processes extxyz files to:
  - Retain one instance per substrate category
  - Generate spreadsheet with duplicate counts

- **energy.py**  
  Extracts free energy values from extxyz files

- **ini_rel.py**  
  Extracts initial/final configurations from trajectory files

### Structure Analysis
- **cdvae_cif.py**  
  Converts crystal files (.cif/.vasp) to custom format and saves to Excel

- **chemical_formula.py**  
  Generates text file mapping CIF filenames to chemical formulas

- **adsorbate.py**  
  Filters CIF files based on atomic conditions and moves qualifying files

### Data Management
- **dict.py**  
  `python dict.py filename.extxyz`  
  Extracts files containing 'N' in adsorbate species

- **copy_file.py**  
  Copies files listed in Excel spreadsheet's first column to target location

### Energy Calculations
- **cal_energy.py**  
  Performs adsorption simulations by:
  1. Creating adsorbate molecules (NH3O, NH2, etc.)
  2. Running slab adsorption
  3. Optimizing structures
  4. Saving energy/structure data

### Utilities
- **series.py**  
  Analyzes file content frequencies and prints statistics

## ⚙️ Configuration Files

**no3rr.yaml**  
Crystal data module configuration containing:
- Dataset paths
- Preprocessing parameters
- Training settings
- Property definitions

## 🙏 Acknowledgments

All contributors for their valuable input

Open-source library maintainers

Research institutions supporting materials informatics
