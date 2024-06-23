classdef QMol_DFT_orbital_basis < QMol_DFT_orbital
%QMol_DFT_orbital_basis density-functional-theory (DFT) orbital(s) with
%   basis-set discretization.
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
methods (Static,Access={?QMol_doc,?QMol_DFT_orbital})
function showInfo
    fprintf('  * QMol_DFT_orbital_basis:\n      > Kohn-Sham orbital\n      > Basis-set discretization\n'); 
    QMol_DFT_orbital_basis.version;
end
end
methods (Access=public)
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components (do not trust KSO, KSOup, or KSOdw, they might
    % have been initialized empty; fetch size from domain and parent DFT)
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(size(obj.disc.v,2),'imag');

    if obj.isSpinPol,                                                                               if opt
        % Spin polarized
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham orbitals (basis set)', 0,1);
        QMol_DFT_profiler.showMemoryFootprint('spin up'     ,numel(obj.disc.QM.occ{1})*mem,2);
        QMol_DFT_profiler.showMemoryFootprint('spin down'   ,numel(obj.disc.QM.occ{2})*mem,2);      end
        
        mem             =   (numel(obj.disc.QM.occ{1}) + numel(obj.disc.QM.occ{2})) * mem;
    else,                                                                                           if opt
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham orbitals',numel(obj.disc.QM.occ)*mem,1);   end
        
        mem             =   numel(obj.disc.QM.occ)*mem;
    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public,Static)
    function b = isBasis,   b   =   true;  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

