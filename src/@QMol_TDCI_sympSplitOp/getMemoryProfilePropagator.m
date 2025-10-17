function mem = getMemoryProfilePropagator(obj,opt)
%getMemoryProfilePropagator estimates the memory foorprint of the TDSE
%   propagator

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nargin < 2,  opt     =   false;  end

    if opt, QMol_SE_profiler.showMemoryFootprint('TDCI propagator (symplectic split operator)', 0,1); end
    
    % Liouville operators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    m                   =   2*QMol_SE_profiler.getMemoryFootprint(numel(obj.CI.CI),'imag');
    mem                 =   m;

end