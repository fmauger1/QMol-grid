classdef QMol_SE_wfcn_basis < QMol_SE_wfcn
%QMol_SE_wfcn_basis Schrodinger-equation (SE) wave function(s) with basis-
%   set discretization.
%
%   NOTES:
%   > One can change the data-equality criterion through the member
%     property eqTol. Side effects to doing so are untested, though (the
%     code is developed with eqTol = 1e-10).
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_SE_wfcn})
function showInfo
    fprintf('  * QMol_SE_wfcn_basis:\n      > SE wave function(s)\n      > Basis-set discretization\n'); 
    QMol_SE_wfcn_basis.version;
end
end
methods (Access=public)
function mem = getMemoryProfile(obj,opt,N)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 3,  N   =   1;
    if nargin < 2,  opt =   false;  end, end

    % Display components (do not trust wfcn, it might have been initialized 
    % empty; fetch size from domain and input number of wave functions)
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(size(obj.disc.v,2),'imag'); if opt
    
    QMol_SE_profiler.showMemoryFootprint('Wave function(s)',N*mem,1);                       end
    mem                 =   N*mem;
    
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public,Static)
    function b = isBasis,   b   =   true;  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

