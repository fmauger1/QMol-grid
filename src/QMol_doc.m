classdef (Abstract) QMol_doc < QMol_suite
%QMol_doc run-time documentation interface for the QMol toolbox. It defines
%   > display list of components
%   > bibliographic entries for kernel components
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release + consolidate bibliography

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access=?QMol_doc)
function showInfo
    fprintf( '  * QMol_doc:\n');
    fprintf(['      > Interface for run-time documentation\n'         ...
             '      > Citable references for kernel components\n']); 
    QMol_doc.version;
end
end
methods (Static,Access=public)
function showComponents
%showComponents displays the list of available components in the QMol
%   toolbox. Components associated with the selected discretization
%   framework are listed in the QMol_info class.
    
    % Header
    QMol_doc.showHeader;

    % Kernel components
    QMol_doc.showSection('Kernel');
    QMol_doc.showInfo;
    QMol_info.showInfo;
    QMol_test.showInfo;
    fprintf('\n')

    fprintf('  Density-functional theory (DFT)\n');
    QMol_DFT.showInfo;
    QMol_DFT_SCF_Anderson.showInfo;
    QMol_DFT_eigs.showInfo;
    QMol_DFT_eig_basis.showInfo;
    QMol_DFT_profiler.showInfo;
    fprintf('\n')

    fprintf('  Time-dependent DFT (TDDFT)\n');
    QMol_TDDFT.showInfo;
    QMol_TDDFT_sympSplitOp.showInfo;
    QMol_TDDFT_SSO_2S.showInfo;
    QMol_TDDFT_SSO_4BM.showInfo;
    QMol_TDDFT_SSO_4FR.showInfo;
    QMol_TDDFT_SSO_6BM.showInfo;
    fprintf('\n')

    fprintf('  Schrodinger equation (SE)\n');
    QMol_SEq.showInfo;
    QMol_SE_eigs.showInfo;
    QMol_SE_eig_basis.showInfo;
    QMol_SE_profiler.showInfo;
    fprintf('\n')

    fprintf('  Time-dependent Schrodinger equation (TDSE)\n');
    QMol_TDSE.showInfo;
    QMol_TDSE_sympSplitOp.showInfo;
    QMol_TDSE_SSO_2S.showInfo;
    QMol_TDSE_SSO_4BM.showInfo;
    QMol_TDSE_SSO_4FR.showInfo;
    QMol_TDSE_SSO_6BM.showInfo;
    fprintf('\n')

    fprintf('  Miscellaneous\n');
    QMol_Vmol.showInfo;
    fprintf('\n')
    
    % Specific discretization components
    QMol_info.showComponents; 

    % Test suite components
    QMol_test.showComponents;

    % References (package is cited by default)
    QMol_doc.showBibliography({});

    % Funding
    QMol_doc.showFunding;
    
    % Footer
    QMol_doc.showFooter;

end
end
%% Run-time documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_suite)
function showHeader
%showHeader displays the header for the QMol toolbox

    % QMol header
    fprintf('########################### QMol-grid package ############################\n');
    fprintf('Quantum simulation methods for atomic and molecular systems [Mauger XXXX].\n\n')

    fprintf('  * Kernel (release)\n');QMol_doc.showKernelHeader; QMol_info.showHeader; fprintf('\n');

    QMol_doc.showFooter; fprintf('\n');

    % License
    QMol_doc.showSection('License');
    QMol_doc.showLicense;
    fprintf('\n')

    % External components
    QMol_doc.showSection('External components');
    fprintf('  * convertUnit\n'); QMol_doc.showVersion(convertUnit.version,convertUnit.lastMod,convertUnit.author);
    fprintf('  * fourierTool\n'); QMol_doc.showVersion(fourierTool.version,fourierTool.lastMod,fourierTool.author);
    fprintf('\n')

end
function showLicense
%showLicense displays the license (2-Clause BSD License)

    fprintf('  Copyright (c) 2024, Francois Mauger\n  All right reserved\n\n');

    fprintf('  Redistribution and use in source and binary forms, with or without\n')
    fprintf('  modification, are permitted provided that the following conditions are\n');
    fprintf('  met:\n\n');

    fprintf('  1. Redistributions of source code must retain the above copyright\n');
    fprintf('     notice, this list of conditions and the following disclaimer.\n\n');

    fprintf('  2. Redistributions in binary form must reproduce the above copyright\n');
    fprintf('     notice, this list of conditions and the following disclaimer in the\n');
    fprintf('     documentation and/or other materials provided with the distribution.\n\n');

    fprintf('  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS\n');
    fprintf('  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED\n');
    fprintf('  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n');
    fprintf('  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n');
    fprintf('  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n');
    fprintf('  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED\n');
    fprintf('  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR\n');
    fprintf('  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF\n');
    fprintf('  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING\n');
    fprintf('  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\n');
    fprintf('  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n');
    
end
function showFunding
%showFunding displays the funding source that supported the development of
%   the toolbox
    
    QMol_doc.showSection('Funding');
    fprintf('    The original development of the QMol-grid toolbox, and its (TD)DFT\n');
    fprintf('  features, was supported by the U.S. Department of Energy, Office of\n');
    fprintf('  Science, Office of Basic Energy Sciences, under Award No. DE-SC0012462.\n');
    fprintf('    Addition of the (TD)SE features was supported by the National Science\n');
    fprintf('  Foundation under Grant No. 2207656.\n\n');
end
function showFooter
%showFooter displays the footer for the QMol toolbox
    fprintf('##########################################################################\n');
end
function showVersion(V,L,A)
%showVersion formats and displays the input version, last modification 
%   date, and author 
    fprintf('    V-%-9s (%-10s) %+45s\n',V,L,A);
end
function showSection(S)
%showSection formats and diplays a section break
    if nargin == 0   ||   isempty(S)
        fprintf('==========================================================================\n');
    else
        X='====================================================================';
        fprintf('=== %s %s\n',S,X(numel(S):end));
    end
end
function showBibliography(ref)
%showBibliography displays the list of cited references (ref)

    % Add reference to package
    if isempty(ref),    ref     =   { 'Mauger XXXX'};
    else,               ref     =   [{'Mauger XXXX'}, ref];                 end

    % Print reference list
        % Section header
        QMol_doc.showSection('References');
        
        % Sort references (and remove potential duplicates)
        ref             =   unique(lower(erase(ref,' ')));
        
        % Print bibliography
        for k = 1:numel(ref)
            QMol_info.(ref{k});
        end

        % Blank line
        fprintf('\n')
end
end
%% Bibliography %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods(Hidden,Static,Access=public,Sealed)                                            %\n (last possible position for line break)
function missingreference
    fprintf( '  [REF] Missing reference.\n');
end
function anderson1965
    fprintf(['  [Anderson 1965] D.G. Anderson, "Iterative procedures for nonlinear\n' ...
             '    integral equations," Journal of the ACM 12, 547 (1965).\n']);
end
function baker2015
    fprintf(['  [Baker 2015] T.E. Baker, E.M. Stoudenmire, L.O. Wagner, K. Burke, and\n'   ...
             '    S.R. White, "One-dimensional mimicking of electronic structure: The\n'   ...
             '    case for exponentials," Physical Review B 91, 235141 (2015).\n']);
end
function blanes2002
    fprintf(['  [Blanes 2002] S. Blanes and P.C. Moan, "Practical symplectic partitioned\n' ...
             '    Runge-Kutta and Runge-Kutta-Nystrom methods," Journal of Computational\n' ...
             '    and Applied Mathematics 142, 313 (2002).\n']);
end
function casa2006
    fprintf(['  [Casa 2006] F. Casas and A. Iserles, "Explicit Magnus expansions for\n'     ...
             '    nonlinear equations," Journal of Physics A: Mathematical and General\n'   ...
             '    39, 5445 (2006).\n']);
end
function forest1990
    fprintf(['  [Forest 1990] E. Forest and R.D. Ruth, "Fourth-order symplectic\n'           ...
             '    integration," Physica D: Nonlinear Phenomena 43, 105 (1990).\n']);
end
function helbig2011
    fprintf(['  [Helbig 2011] N. Helbig, J.I. Fuks, M. Casula, M.J. Verstraete,\n'  ...
             '    M.A.L. Marques, I.V. Tokatly, and A. Rubio, "Density functional theory\n'...
             '    beyond the linear regime: Validating an adiabatic local density\n'...
             '    approximation," Physical Review A 83, 032503 (2011).\n']);
end
function javanainen1988
    fprintf(['  [Javanainen 1988] J. Javanainen, J.H. Eberly, and Q. Su, "Numerical\n' ...
             '    simulations of multiphoton ionization and above-threshold electron\n' ...
             '    spectra," Physical Review A 38, 3430 (1988).\n']);
end
function johnson1988
    fprintf(['  [Johnson 1988] D. D. Johnson, "Modified Broyden''s method for\n' ...
             '    accelerating convergence in self-consistent calculations," Physical\n'...
             '    Review B 38, 12807 (1988).\n']);
end
function kohn1965
    fprintf(['  [Kohn 1965] W. Kohn and L.J. Sham, "Self-Consistent Equations Including\n' ...
             '    Exchange and Correlation Effects," Phys. Rev. 140, A1133 (1965).\n']);
end
function legrand2002
    fprintf(['  [Legrand 2002] C. Legrand, E. Suraud, and P.-G. Reinhard, "Comparison of\n' ...
             '    self-interaction-corrections for metal clusters," J. Phys. B 35, 1115\n' ...
             '    (2002).\n']);
end
function mauger2024
    fprintf(['  [Mauger 2024] F. Mauger, C. Chandre, M.B. Gaarde, K. Lopata, and K.J. \n' ...
             '    Schafer, "Hamiltonian  formulation and symplectic split-operator \n' ...
             '    schemes for time-dependent  density-functional-theory equations of \n' ...
             '    electron dynamics in molecules," Communications in Nonlinear Science \n' ...
             '    and Numerical Simulation 129, 107685 (2024).\n']);
end
function maugerxxxx
    fprintf(['  [Mauger XXXX] F. Mauger, et al., "QMol-grid: A MATLAB package for\n' ...
             '    quantum-mechanical simulations in atomic and molecular systems," \n' ...
             '    arXiv:xxx (XXXX).\n']);
end                                                                                    %\n (last possible position for line break)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

