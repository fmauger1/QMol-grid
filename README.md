# QMol-grid: A MATLAB package for quantum-mechanical simulations in atomic and molecular systems

The `QMol-grid` package provides a suite of routines for performing quantum-mechanical simulations in atomic and molecular systems. All simulations use an underlying Cartesian-grid discretization scheme. The `QMol-grid` package is written using MATLAB's object-oriented features and handle classes.

The `QMol-grid` package provides support for the following types of computations:
- `DFT` Ground-state density-functional theory, both using a Cartesian grid or basis-set discretization
- `HF` Ground-state Hartree Fock, both using a Cartesian grid or basis-set discretization
- `SE` Ground-state Schrödinger equation, both using a Cartesian grid or basis-set discretization
- `TDDFT` Real-time time-dependent density-functional theory, using a Cartesian grid
- `TDSE` Time-dependent Schrödinger equation, using a Cartesian grid

___
## Table of Contents
  * [Installation](#installation)
  * [Example 1: Schrödinger-equation ground state ](#example1)
  * [Example 2: ](#example2)
  * [Reference](#reference)
  
___

## Installation
- download the zipped file in the release, and unzip it in a folder of the type `<user>/Documents/MATLAB/QMol-grid`
- add permanently the folder `<user>/Documents/MATLAB/QMol-grid` (without its subfolders) to the [MATLAB path](https://mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html)
- after successful installation, the package documentation will be accessible in MATLAB's, in the "Supplemental Software" section

### Tests

Display the list of components in the `QMol-grid` package
```Matlab
QMol_doc.showComponents;
```

The `QMol-grid` package provides a suite of tests that performs low-level checks on the various properties that are defined throughout its components. The test suite can be used to check basic functions after installing the package or modifying its components.

Run all available unit tests using
```Matlab
QMol_test.test;
```

### Limitations

- The current release on supports one-dimensional computations 
- Time-dependent Hartree Fock is not currently available
- Time propagation on basis sets is not currently available 

[&uarr;](#table-of-contents)
___
## <a name="example1"></a>Example 1: Schrödinger-equation ground state 

Here we illustrate how to use the `QMol-grid` package to calculate the ground-state wave function of a one-dimensional hydrogen-like atom. The Schrödinger-equation ground-state corresponds to the lowest-energy solution to the eigenvalue problem $\hat{\mathcal{H}}\psi(x)=E\psi(x)$, where $\hat{\mathcal{H}}$ is the Schrödinger-equation Hamiltonian operator, $\psi$ is the wave function, and  $E$ its associated energy. In atomic units, the Hamiltonian operator is $\hat{\mathcal{H}} = -\frac{\Delta}{2} + \hat{\mathcal{V}}$.

Specifically, it walks through defining (i) the domain and grid discretization over which the Schrödinger-equation and wave function are calculated, the (ii) atomic potential and (iii) Schrödinger-equation model, and calculating (iv) the ground state associated with these properties.

We model the one-dimensional hydrogen model atom using a soft-Coulomb potential $V(x)=-1/\sqrt{x^2+a^2}$ with
```Matlab
H = QMol_Va_softCoulomb('softeningParameter',sqrt(2));
```
where  we choose the softening parameter $a=\sqrt{2}$ to match H's ground state energy. By default, the atom is located at the origin $x=0$. Note that H only corresponds to the atomic model, which is shared with molecular systems and various quantum frameworks. Thus, it must be turned into a valid Schrödinger-equation potential, using
```Matlab
V = QMol_SE_V('atom',H);
```

The simulation domain must be a Cartesian grid -- with all increasing, equally spaced discretization points -- and should be wide and with small enough of a discretization step to properly capture the wave function. In our case we select a domain ranging -15 to 15 a.u., with a discretization steps of 0.1 a.u.
```Matlab
x = -15:.1:15;
```

We now have all the elements to define a Schrödinger-equation model object with the potential and domain defined above
```Matlab
SE = QMol_SE(                        ...
         'xspan',                x,  ...
         'potential',            V);
```

We next move to calculating its associated ground-state wave function and energy using the two commands
```Matlab
GSS = QMol_SE_eigs;
GSS.computeGroundState(SE);
```

At the end of the calculation, the ground-state wave function is stored in the input Schrödinger-equation object SE, together with relevant information such as the domain discretization. For instance, solely relying on SE, one can plot the ground-state wave function with
```Matlab
figure
    plot(SE.xspan,SE.waveFunction.waveFunction,'-','LineWidth',2)
    set(gca,'box','on','FontSize',12,'LineWidth',2)
    xlabel('x (a.u.)')
    ylabel('wave function (a.u.)')
    xlim(SE.xspan([1 end]))
```
producing
<p align="center">
  <img src="https://github.com/fmauger1/QMol-grid/blob/main/GS__T01.png" alt="Example 1" width="300"/>
</p>

From the `plot` command line, we see that the domain-discretization grid may be recovered using the xspan property in the object SE (using the standard object-oriented dot notation `SE.xspan`). On the other hand, the wave function is nested inside another object, which explains the consecutive dots `SE.waveFunction.waveFunction`. 

[&uarr;](#table-of-contents)
___
## <a name="example2"></a>Example 2:

[&uarr;](#table-of-contents)
___
## Reference
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
The original development of the `QMol-grid` package was supported by the U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences, under Award No. DE-SC0012462.
Addition of the (TD)SE features was supported by the National Science Foundation under Grant No. 2207656

[&uarr;](#table-of-contents)
