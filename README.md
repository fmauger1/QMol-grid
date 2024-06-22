# QMol-grid: A MATLAB package for quantum-mechanical simulations in atomic and molecular systems

![Version](https://img.shields.io/badge/version-v1.21.000-blue)
[![Generic badge](https://img.shields.io/badge/MATLAB-R2022a-orange)](https://fr.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/license-BSD-lightgray)](https://github.com/fmauger1/QMol-grid/blob/main/LICENSE)

The `QMol-grid` package provides a suite of routines for performing quantum-mechanical simulations in atomic and molecular systems in one spatial dimension. All simulations use an underlying Cartesian-grid discretization scheme. The `QMol-grid` package is written using MATLAB's object-oriented features and handle classes.

The `QMol-grid` package provides support for the following types of computations:
- `DFT` Ground- and excited-state density-functional theory, both using a Cartesian grid or basis-set discretization
- `HF` Ground- and excited-state Hartree Fock, both using a Cartesian grid or basis-set discretization
- `SE` Ground- and excited-state Schrödinger equation, both using a Cartesian grid or basis-set discretization
- `TDDFT` Real-time time-dependent density-functional theory, using a Cartesian grid
- `TDSE` Time-dependent Schrödinger equation, using a Cartesian grid

The time-propagators are computed using symplectic-split operators (2nd, 4th, and 6th order in time, and spectral in space). They support field-free and laser-driven simulations (in the dipole approximation) with the following on-the-fly features, each specifying their own time sampling
- Checkpointing, with the creation of a restart MATLAB file that can be used to resume a calculation that was stopped before it was finished
- Calculation and storage of the dipole, dipole velocity, and dipole acceleration signals
- Calculation and storage of the wave function(s)/Kohn-Sham orbitals and Hamiltonian-component energies
- Storage of the wave function(s) (Schrödinger), and Kohn-Sham orbitals and one-body density (DFT)
- Calculation and storage of the ionization signal, keeping track of how much electronic density is absorbed at the domain boundaries
- Calculation and storage of the results of installable output functions of the wave function(s) (Schrödinger), and Kohn-Sham orbitals or one-body density (DFT)
- Saving the intermediate Schrödinger- or DFT-model objects in separate MATLAB files (`.mat`)

___
## Table of Contents
  * [Installation](#installation)
  * [Example 1: Schrödinger-equation ground state ](#example1)
  * [Example 2: Time-dependent density-functional theory](#example2)
  * [Reference](#reference)
  
___

## Installation
- from a release: download the zipped file in the release, and unzip the content of the folder `/src/`, including all its subfolders, in a folder of the type `<user>/Documents/MATLAB/QMol-grid` or
- from the github repository (development version): download the entire content of the folder `/src/`, including all its subfolders, in a folder of the type `<user>/Documents/MATLAB/QMol-grid`
- add permanently the folder `<user>/Documents/MATLAB/QMol-grid` (without its subfolders) to [MATLAB path](https://mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html)
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

- The current release only supports one-dimensional computations 
- Time-dependent Hartree Fock is currently not available
- Time propagation on basis sets is currently not available 

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

The simulation domain must be a Cartesian grid -- with all increasing, equally spaced discretization points -- and should be wide enough and with small enough of a discretization step to properly capture the wave function. In our case we select a domain ranging -15 to 15 a.u., with a discretization steps of 0.1 a.u.
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

From the `plot` command line, we see that the domain-discretization grid may be recovered using the xspan property in the object SE (using the standard object-oriented dot notation `SE.xspan`). On the other hand, the wave function is nested inside another object, which explains the consecutive dots `SE.waveFunction.waveFunction`. Other properties in the object `SE.waveFunction` are used by ground/excited-state and TDSE calculations; we refer to the `QMol_SE_wfcn` documentation page for further details.

[&uarr;](#table-of-contents)
___
## <a name="example2"></a>Example 2: Time-dependent density-functional theory

For a given set of initial Kohn-Sham orbitals $\phi_{k}$, the TDDFT dynamics is described by the nonlinear system of partial differential equations, in atomic units (a.u.)
$$i\partial_t \phi_k({\bf x};t) =\hat{\mathcal{H}}_{\rm DFT}[\{\phi_k\}_k;t]({\bf x};t)\ \phi_k({\bf x};t)  .\quad\quad\quad (2.1)$$
The `QMol-grid` package relies on the canonical Hamiltonian structure of TDDFT [Mauger 2024](https://doi.org/10.1016/j.cnsns.2023.107685) to integrate the dynamics of equation (2.1). 
In this example, we illustrate how to use the `QMol-grid` package to integrate the TDDFT dynamics of an open-shell one-dimentional molecular ion model with 3 atomic centers and 5 active electrons.

### Initial condition
In the `QMol-grid` package, TDDFT simulations are decoupled from setting up the initial condition, which must be done independently. For our example, we start by calculating the neutral-molecule ground state:
```Matlab
% Molecular model
V_1     =   QMol_Va_softCoulomb('atom','X_1','charge',2,'position',-3);
V_2     =   QMol_Va_softCoulomb('atom','X_2','charge',2,'position', 0);
V_3     =   QMol_Va_softCoulomb('atom','X_3','charge',2,'position', 3);
 
% DFT model
Vext    =   QMol_DFT_Vext('atom',{V_1,V_2,V_3});
Vh      =   QMol_DFT_Vh_conv;
Vxc     =  {QMol_DFT_Vx_LDA_soft,QMol_DFT_Vc_LDA_soft};
 
DFT     =   QMol_DFT_spinPol(                                       ...
                'xspan',                       -50:.1:50,           ...
                'occupation',                  {[1 1 1],[1 1 1]},   ...
                'externalPotential',            Vext,               ...
                'HartreePotential',             Vh,                 ...
                'exchangeCorrelationPotential', Vxc,                ...
                'selfInteractionCorrection',    'ADSIC'             );
 
% DFT ground state
SCF     =   QMol_DFT_SCF_Anderson;
SCF.solveSCF(DFT);
```
Next, we manually induce an excitation in the molecular cation by successively (i) replacing one of the Kohn-sham orbitals by a superposition of molecular-orbital states (excitation part) and (ii) removing an electron, going from 3 to 2, from the down-spin Kohn-Sham orbitals (ionization part).
```Matlab
% Induce excitation
DFT.orbital.set('orbitalDown',[DFT.KSO.KSOdw(:,1) (DFT.KSO.KSOdw(:,2)+DFT.KSO.KSOdw(:,3))/sqrt(2)]);

% Induce ionization
DFT.set('occupation',{[1 1 1],[1 1]});
```
We now have a non-stationary set of Kohn-Sham orbitals, leading to field-free dynamics under the flow of equation (2.1).

### TDDFT simulation
With the DFT molecular model, including the initial condition, in hand, we now move to integrating the subsequent field-free TDDFT dynamics. For this, we select a [fourth-order Forest Ruth](https://doi.org/10.1016/0167-2789(90)90019-L) symplectic split-operator scheme (see [Mauger 2024](https://doi.org/10.1016/j.cnsns.2023.107685) for more details).
```Matlab
TDDFT   =   QMol_TDDFT_SSO_4FR(                     ...
                'time',                 0:10:100,   ...
                'timeStep',             2e-2,       ...
                'saveDensity',          true,       ...
                'saveDensityTime',      1);
```
In our example, the TDDFT object is created with:
- The first pair of arguments specifies that the integration should start at time t=0 and end at t=100 a.u. The step of 10 a.u., is unrelated to the propagation time step and instead specifies the time intervals to use in progress display.
- The second pair of arguments specifies the (fixed) time step for the propagation.
- The third pair of arguments indicates that the one-body density should be saved periodically, with the period specified by the fourth pair of arguments, _i.e._, every 1 a.u. in our case.

Then, we launch the TDDFT integration with 
```Matlab
TDDFT.propagate(DFT);
```
At the end of the simulation, the DFT object has been updated to contain the Kohn-Sham orbitals at t=100 a.u. The time-dependent one-body density is stored in the TDDFT object itself.

### Plotting the result
Next we recover calculated observables out of the TDDFT object. Each set of observable is stored in a separate structure property in the TDDFT object, which containts (i) the exact time vector at which the quantity has been saved and (ii) the observable itself. In our case, the structure of interest is `TDDFT.outDensity` with the up- and down-spin densities respectively stored in the fields `totalUp` and `totalDown`. The densities are matrices with columns corresponding to the successive saved times.
To plot the spin density, defined as the difference between the up- and down-spin one-body densities, we use
```Matlab
figure
    imagesc(TDDFT.outDensity.time,DFT.xspan,TDDFT.outDensity.totalUp-TDDFT.outDensity.totalDown)
    set(gca,'box','on','FontSize',12,'LineWidth',2,'YDir','normal')
    xlim(TDDFT.outDensity.time([1 end]))
    ylim([-10 10])
    xlabel('time (a.u.)')
    ylabel('position (a.u.)')
    title('spin density')
    colorbar vert
```
producing 
<p align="center">
  <img src="https://github.com/fmauger1/QMol-grid/blob/main/GS__T07.png" alt="Example 2" width="300"/>
</p>

[&uarr;](#table-of-contents)
___
## Reference
- F. Mauger and C. Chandre, *QMol-grid: A MATLAB package for quantum-mechanical simulations in atomic and molecular systems*, [arXiv:24](https://arxiv.org/abs/24)
```bibtex
@unpublished{mauger2024,
  title = {QMol-grid: A MATLAB package for quantum-mechanical simulations in atomic and molecular systems},
  author = {Mauger, Fran\c{c}ois and Chandre, Cristel},
  year = {2024},
  URL = {https://arxiv.org/abs/24}
}
```
For more information: <fmauger@lsu.edu>


## Acknowledgments
The original development of the `QMol-grid` package, and its (TD)DFT features, was supported by the U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences, under Award No. DE-SC0012462.
Addition of the (TD)SE features was supported by the National Science Foundation under Grant No. 2207656

[&uarr;](#table-of-contents)
