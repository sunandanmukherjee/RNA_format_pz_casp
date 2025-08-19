# RNA_format_pz_casp
A command-line pipeline for formatting 3D RNA models to meet the strict submission guidelines of the RNA-Puzzles and CASP experiments.<p align="center">
<img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
<img src="https://img.shields.io/badge/python-3.7+-blue.svg" alt="Python 3.7+">
<img src="https://img.shields.io/badge/dependency-Biopython-green.svg" alt="Dependency: Biopython">
</p>

Dependencies
Biopython

How to use:
1. Download the template file from the official website of RNA-Puzzles or CASP for a target
2. Split the files with "--model" option using the following command:

   ```
   python split_PDB.py <filename.pdb> --model
   ```
3. Use format_PDB.py to adjust the file based on the official template:

   ```
   python format_PDB.py --template <filename_model1.pdb> --model <pz_model1.pdb>
   ```

   
