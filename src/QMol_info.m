classdef(Abstract) QMol_info < QMol_doc
%QMol_info run-time documentation for implementation-specific components of
%   the QMol-grid package.
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release + consolidate bibliography
%   01.22.000   04/19/2025  F. Mauger
%       Integrale CI module (+ new version number)
%   01.22.001   04/30/2025  F. Mauger
%       Integrate CI module split
%   01.23.000   06/04/2025  F. Mauger
%       Add QMol_TDCI_abs_CAP and QMol_TDCI_abs_mask
%       Add QMol_TDCI_SSO
%   01.23.001   06/06/2025  F. Mauger
%       Add QMol_TDCI_abs_CAP_space
%   01.23.002   06/07/2025  F. Mauger
%       Add QMol_TDCI_abs_mask_space
%   01.23.003   07/22/2025  F. Mauger
%       Add QMol_TDDFT_ESSO

% 2-Clause BSD License
%
% Copyright (c) 2025, Francois Mauger, all right reserved
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.003','07/22/2025','F. Mauger')
end
end
methods (Static,Access=?QMol_doc)
function showInfo
    fprintf( '  * QMol_info:\n');
    fprintf(['      > Implementation-specific run-time documentation\n'         ...
             '      > Citable references for implementation-specific components\n']); 
    QMol_info.version;
end
end
methods (Static,Access=?QMol_suite)
function showHeader
%showHeader displays the discretization input for the header

    fprintf('  * One-dimensional model (release)\n');
    QMol_info.version;

end
end
methods (Static,Access=public)
function showComponents
%showComponents displays the list of available components in the QMol
%   toolbox. Components associated with the selected discretization
%   framework are listed in the QMol_info class.


    % Domain discretization ===============================================
    QMol_doc.showSection('Domain discretization');
    QMol_disc.showInfo;
    QMol_disc_basis.showInfo;
    %QMol_data.showInfo;
    fprintf('\n')

    % DFT =================================================================
    QMol_doc.showSection('DFT theory');
    fprintf('  Model\n');
    QMol_DFT_spinPol.showInfo;
    QMol_DFT_spinRes.showInfo;
    fprintf('\n')

    fprintf('  Components\n');
    QMol_DFT_density.showInfo;
    QMol_DFT_orbital.showInfo;
    QMol_DFT_orbital_basis.showInfo;
    QMol_DFT_Vks.showInfo;
    QMol_DFT_Vks_basis.showInfo;
    QMol_DFT_Vks_grad.showInfo;
    fprintf('\n')

    fprintf('  Functionals\n');
    QMol_DFT_Vext.showInfo;
    QMol_DFT_Vh_conv.showInfo;
    QMol_DFT_Vh_fft.showInfo;
    QMol_DFT_Vx_LDA_exp.showInfo;
    QMol_DFT_Vx_LDA_soft.showInfo;
    QMol_DFT_Vx_XX_conv.showInfo;
    QMol_DFT_Vx_XX_fft.showInfo;
    QMol_DFT_Vc_LDA_soft.showInfo;
    fprintf('\n')

    % TDDFT ===============================================================
    QMol_doc.showSection('TDDFT theory');

    fprintf('  Absorbing boundaries\n');
    QMol_TDDFT_abs_CAP.showInfo;
    QMol_TDDFT_abs_mask.showInfo;
    fprintf('\n')

    fprintf('  Propagators\n');
    QMol_TDDFT_SSO.showInfo;
    QMol_TDDFT_ESSO.showInfo;
    fprintf('\n')

    % CI ==================================================================
    QMol_doc.showSection('Configuration-interaction theory');
    fprintf('  Model\n');
    QMol_CI_conv.showInfo;
    fprintf('\n')

    % TDCI ================================================================
    QMol_doc.showSection('TDCI theory');

    fprintf('  Absorbing boundaries\n');
    QMol_TDCI_abs_CAP.showInfo;
    QMol_TDCI_abs_CAP_space.showInfo;
    QMol_TDCI_abs_mask.showInfo;
    QMol_TDCI_abs_mask_space.showInfo;
    fprintf('\n');

    fprintf('  Propagators\n');
    QMol_TDCI_SSO.showInfo;
    fprintf('\n');

    % SE ==================================================================
    QMol_doc.showSection('Schrodinger-equation theory');
    fprintf('  Model\n');
    QMol_SE.showInfo;
    fprintf('\n')

    fprintf('  Components\n');
    QMol_SE_wfcn.showInfo;
    QMol_SE_wfcn_basis.showInfo;
    QMol_SE_V.showInfo;
    fprintf('\n')

    % TDSE ================================================================
    QMol_doc.showSection('TDSE theory');

    fprintf('  Absorbing boundaries\n');
    QMol_TDSE_abs_CAP.showInfo;
    QMol_TDSE_abs_mask.showInfo;
    fprintf('\n');

    fprintf('  Propagators\n');
    QMol_TDSE_SSO.showInfo;
    fprintf('\n');

    % External fields =====================================================
    QMol_doc.showSection('External field');
    QMol_extField_dipole.showInfo;
    fprintf('\n')

    % Pseudopotentials ====================================================
    QMol_doc.showSection('Pseudopotentials');
    QMol_Va_Gaussian.showInfo;
    QMol_Va_softCoulomb.showInfo;
    fprintf('\n')

end
end
%% Bibliography %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods(Hidden,Static,Access=public,Sealed)                                            %\n (last possible position for line break)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

