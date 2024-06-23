function mem = getMemoryProfilePropagator(obj,opt)
%getMemoryProfilePropagator estimates the memory foorprint of the TDSE
%   propagator

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nargin < 2,  opt     =   false;  end

    if opt, QMol_SE_profiler.showMemoryFootprint('TDSE propagator (symplectic split operator)', 0,1); end
    
    % Liouville operators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    m                   =   2*QMol_SE_profiler.getMemoryFootprint(numel(obj.SE.disc.x),'imag');
    mem                 =   m;

    if opt, QMol_SE_profiler.showMemoryFootprint('Liouville operators', m,2); end,              if obj.SE.dim  == 1
    
    % Miscellaneous ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    m                   =   numel(obj.SE.disc.x);                                               elseif  obj.SE.dim == 2
    m                   =   2*numel(obj.SE.disc.x)*numel(obj.SE.disc.y);                        elseif  obj.SE.dim == 3
    m                   =   3*numel(obj.SE.disc.x)*numel(obj.SE.disc.y)*numel(obj.SE.disc.z);    end
    m                   =   QMol_SE_profiler.getMemoryFootprint(m,'real');
    mem                 =   mem + m;

    if opt, QMol_SE_profiler.showMemoryFootprint('Dipole operator', m,2); end

end