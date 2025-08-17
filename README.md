# Atomic Mass Center Finder (DIPOL tag)
Automatically find the atomic mass center (DIPOL tag) from POSCAR and POTCAR.

## Overview

This project provides a Python-based tool to **calculate the DIPOL mass center** in VASP.

It automatically reads `POSCAR` and `POTCAR` files from the **current working directory** and computes the atomic mass center while properly handling periodic boundary conditions (PBC).

## Features

- Reads `POSCAR` and `POTCAR` directly from the current directory.
- Calculates the atomic mass center (DIPOL tag) with **periodic boundary handling**:
  * Boundary atoms are excluded, since their equivalents exist on the opposite side of the periodic cell.

## Usage

1. Download the `dipol.py` file.

2. Place your `POSCAR` and `POTCAR` files in the working directory.

3. Run the notebook to calculate the DIPOL mass center.


## License

MIT License
