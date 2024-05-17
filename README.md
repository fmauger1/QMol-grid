# QMol-grid: A MATLAB package for quantum-mechanical simulations in atomic and molecular systems

The QMol-grid package provides a suite of routines for performing quantum-mechanical simulations in atomic and molecular systems. All simulations use an underlying Cartesian-grid discretization scheme. The QMol-grid package is written using MATLAB's object-oriented features and handle classes.

The QMol-grid package provides support for the following types of computations:
- `DFT` Ground-state density-functional theory, both using a Cartesian-grid or basis-set discretization
- `HF` Ground-state Hartree Fock, both using a Cartesian-grid or basis-set discretization
- `SE` Ground-state Schrodinger equation, both using a Cartesian-grid or basis-set discretization
- `TDDFT` Real-time time-dependent density-functional theory, both using a Cartesian-grid or basis-set discretization
- `TDHF` Time-dependent Hartree Fock, for basis-set discretization only
- `TDSE` Time-dependent Schrodinger equation, both using a Cartesian-grid or basis-set discretization


## Installation
- dowload zipped file in the release, and unzip it in a folder of the type *User*/Documents/MATLAB/QMol-grid
- add permanently the folder *User*/Documents/MATLAB/QMol-grid to the MATLAB path

#### External components
The QMol-grid package uses the following two external components (included in the Tools toolbox, version 01.00)
- `fourierTool`, used with fast-Fourier transforms (e.g., when computing derivatives or convolutions)
- 'convertUnits', used to perform miscellaneous unit conversions

## Example 1:


## Example 2:

## References
- F. Mauger *et al*, *QMol-grid: A MATLAB package for quantum-mechanical simulations in atomic and molecular systems*, [arXiv:24](https://arxiv.org/abs/24)
```bibtex
@unpublished{mauger2024,
  title = {QMol-grid: A MATLAB package for quantum-mechanical simulations in atomic and molecular systems},
  author = {Mauger, F. et al},
  year = {2024},
  URL = {https://arxiv.org/abs/24}
}
```
For more information: <fmauger@lsu.edu>


## Acknowledgments
The original development of the 1D QMol-grid package was supported by the U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences, under Award No. DE-SC0012462.
