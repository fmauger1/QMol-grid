# QMol-grid: A MATLAB package for quantum-mechanical simulations in atomic and molecular systems

The QMol-grid package provides a suite of routines for performing quantum-mechanical simulations in atomic and molecular systems. All simulations use an underlying Cartesian-grid discretization scheme. The QMol-grid package is written using MATLAB's object-oriented features and handle classes.

The QMol-grid package provides support for the following types of computations:
- `DFT` Ground-state density-functional theory, both using a Cartesian-grid or basis-set discretization
- `HF` Ground-state Hartree Fock, both using a Cartesian-grid or basis-set discretization
- `SE` Ground-state Schrodinger equation, both using a Cartesian-grid or basis-set discretization
- `TDDFT` Real-time time-dependent density-functional theory, using a Cartesian-grid
- `TDSE` Time-dependent Schrodinger equation, using a Cartesian-grid


## Installation
- dowload zipped file in the release, and unzip it in a folder of the type *User*/Documents/MATLAB/QMol-grid
- add permanently the folder *User*/Documents/MATLAB/QMol-grid to the MATLAB path

### External components
The QMol-grid package uses the following two external components (included in the Tools toolbox, version 01.00)
- `fourierTool`, used with fast-Fourier transforms (e.g., when computing derivatives or convolutions)
- 'convertUnits', used to perform miscellaneous unit conversions

### Tests

Display the list of components in the QMol-grid package
```Matlab
QMol_doc.showComponents;
```

The QMol-grid package provides a suite of tests that performs low-level checks on the various properties that are defined throughout its components. The test suite can be used to check basic functions after installing the package or modifying its components.

Run all available unit tests using
```Matlab
QMol_test.test;
```

### Limitations

The current release is available for one-dimensional computations. 
Time-dependent Hartree Fock, for basis-set discretization only is not currently available
Time propagation on basis sets is not available 

## Example 1: Schrödinger Equation ground state 

This tutorial walks through the process of setting up and calculating the Schrödinger-equation ground state of a one-dimensional hydrogen atom model. The Schrödinger-equation ground-state corresponds to the lowest energy solution to the eigenvalue problem $\hat{H}\psi(x)=E\psi(x)$, where $\hat{\mathcal{H}}$ is the Schrödinger-equation Hamiltonian operator, $\psi$ is the wave function, and  $E$ is its associated energy. In atomic units, the Hamiltonian operator is $\hat{\mathcal{H}} = -\frac{\Delta}{2} + \hat{\mathcal{V}}.$

In this tutorial, we illustrate how to use the QMol-grid package to calculate the ground-state wave function of a one-dimensional hydrogen-like atom. Specifically, it walks through defining (i) the domain and grid-discretization over which the Schrodinger-equation and wave function are calculated, the (ii) atomic potential and (iii) Schrödinger-equation model, and calculating (iv) the ground state associated with these properties.

We model the 1D hydrogen model atom using a soft-Coulomb potential with
```Matlab
H = QMol_Va_softCoulomb('softeningParameter',sqrt(2));
```
where  we choose the softening parameter $\sqrt{2}$ to match H's ground state energy. By default, the atom is located at the origin $x=0$.
Note that H only corresponds to the atomic model, which is shared with molecular systems and various quantum frameworks. Thus, it must be turned into a valid Schrodinger-equation potential, using
```Matlab
V = QMol_SE_V('atom',H);
```

The simulation domain must be a Cartesian grid -- with all increasing, equally spaced discretization points -- and should be wide and with small enough of a discretization step to properly capture the system's wave function. In our case we select a domain ranging -15 to 15 a.u., with a discretization steps of 0.1 a.u.
```Matlab
x = -15:.1:15;
```

We now have all the elements to define a Schrodinger-equation model object with the potential and domain defined above
```Matlab
SE = QMol_SE(                        ...
         'xspan',                x,  ...
         'potential',            V);
```

With the Schrodinger-equation object defined above, we next move to calculating its associated ground-state wave function and energy using the two commands
```Matlab
GSS = QMol_SE_eigs;
GSS.computeGroundState(SE);
```

At the end of the calcultation, the ground-state wave function is stored in the input Schrodinger-equation object `SE`, together with relevant information such as the domain discretization. For instance, solely relying on `SE`, one can plot the ground-state wave function with
```Matlab
figure
    plot(SE.xspan,SE.waveFunction.waveFunction,'-','LineWidth',2)
    set(gca,'box','on','FontSize',12,'LineWidth',2)
    xlabel('x (a.u.)')
    ylabel('wave function (a.u.)')
    xlim(SE.xspan([1 end]))
```
producing (note that the ground-state calculation start from a random seed and thus the resulting wave function is defined with an arbitrary sign that can change from calculation to calculation)

<p align="center">
  <img src="https://github.com/fmauger1/QMol-grid/blob/main/GS__T01.png" alt="Example 1" width="300"/>
</p>

From the `plot` command line, we see that the domain-discretization grid may be recovered using the xspan property in the object `SE` (using the standard object-oriented dot notation `SE.xspan`). On the other hans, the wave function is nested inside another object, which explains the consecutive dots `SE.waveFunction.waveFunction`. 

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
