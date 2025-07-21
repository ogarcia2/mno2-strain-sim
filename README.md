# MnO₂ Phase Transition Strain Simulations

This repository contains code and data for computing 2D strain fields and deformation gradients during phase transitions between MnO₂ polymorphs. It supports simulations using input from first-principles calculations and image-based measurements.

## Structure

- `src/`: Python scripts for computing deformation gradients (F) and Green-Lagrange strain tensors (E₀)
- `data/`: Atom position data (e.g., from POSCARs)
- `results/`: Output files including tidy F and E₀ values
- `figures/`: Figures used in the report

## Scripts

- `peng_pixel_strain_analysis.py`: Calculates F and E₀ from pixel measurements (Peng et al.)
- `ma_f_e0_calculator.py`: Computes 3D F and E₀ from POSCAR-based atomic positions
- `ma_projection_to_2D.py`: Projects 3D structures into 2D (x–z plane) and recalculates F and E₀

## Requirements

- Python 3.11+
- numpy, pandas

Install:
```bash
pip install -r requirements.txt
Authors
Oskar K. Garcia — San Francisco State University
Dr. Nicole Adelstein — LLNL
Dr. Bo Wang — LLNL
