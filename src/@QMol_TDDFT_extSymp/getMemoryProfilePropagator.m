function mem = getMemoryProfilePropagator(obj,opt)
%getMemoryProfilePropagator estimates the memory foorprint of the TDDFT
%   propagator

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nargin < 2,  opt     =   false;  end

    if opt, QMol_DFT_profiler.showMemoryFootprint('TDDFT propagator (symplectic split operator)', 0,1); end
    
    % Liouville operators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    m                   =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.DFT.disc.x),'imag');
    if obj.DFT.isSpinPol,   m   =   3*m;            else,   m   =   2*m;    end
    mem                 =   m;

    if opt, QMol_DFT_profiler.showMemoryFootprint('Liouville operators', m,2); end

    % Extended phase space ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    m                   =   2*obj.DFT.KSO.getMemoryProfile(false);
    mem                 =   mem + m;
    
    if opt, QMol_DFT_profiler.showMemoryFootprint('Extended phase space', m,2); end,                if obj.DFT.dim  == 1
    
    % Miscellaneous ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    m                   =   numel(obj.DFT.disc.x);                                                  elseif  obj.DFT.dim == 2
    m                   =   2*numel(obj.DFT.disc.x)*numel(obj.DFT.disc.y);                          elseif  obj.DFT.dim == 3
    m                   =   3*numel(obj.DFT.disc.x)*numel(obj.DFT.disc.y)*numel(obj.DFT.disc.z);    end
    m                   =   QMol_DFT_profiler.getMemoryFootprint(m,'real');
    mem                 =   mem + m;

    if opt, QMol_DFT_profiler.showMemoryFootprint('Dipole operator', m,2); end

end