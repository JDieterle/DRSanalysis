# DRSanalysis
Skripts for extracting the optical properties from differential reflectance spectra (DRS) by kramers-kronig constrained variational analysis. This program is able to
 iteratively account for significant absorptions at energies above the measured range.

Usage:
required files (examples can be found in the folder 'example'):
- ascii file with tabulated optical constants (e1 e2 versus energy in eV) of the substrate: 1. col.: Energy, 2. col.: e1, 3. col.: e2 
- ascii file with tabulatd DRS data: 1.-4. row: thickness [nm], 5+.row: data, 1. col.: Energy, 2. col.: wvl, 3+. col.: DRS data
- optional: ascii file containg initial parameter: 1. col.: Energy, 2. col.: e1, 3. col.: e2 

To start the graphical user interface run kk_analysis_drs_GUI.m and follow the instructions. The results are automatically saved in ascii files.
