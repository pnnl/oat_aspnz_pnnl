# oat_aspnz_pnnl

This is the code repository for the ASPNZ software module.  ASPNZ is an extension of the Open Applied Toplogy (OAT) library. 

The purpose of ASPNZ is to support network signal processing by accelerating zigzag homology computation.

## Deprecation warning

This repository is not maintained. Functionality from this project will be integrated directly into future versions of the Open Applied Topology software.

## Installation

All packages are installed and used in Python.

- OAT with Lamellar: cd into `aspnz_benchmarks/software/oat_aspnz_lamellar/oat_python_ZIGZAG`. Follow the installation instructions in the README.md contained in that folder.
- OAT without Lamellar: cd into `aspnz_benchmarks/software/oat_aspnz_serial/oat_python_ZIGZAG`. Follow the installation instructions in the README.md contained in that folder. **If you are using a local machine, such as a laptop, open `aspnz_benchmarks/software/oat_aspnz_serial/oat_rust_ZIGZAG` and delete `, features = ["enable-rofi-shared"]` from the line for Lamellar under the list of dependencies.**