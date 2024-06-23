function rho = DFT_allocateDensity(obj,rho,isFull)
%DFT_allocateDensity allocates a one-body-density (rho) object. Optionally,
%   if an old density object is given as an input, it is reinitialized
%   instead. Optionally, specify whether the density-discretization should
%   be allocated too (to zero)
    
    % Initialization
    if nargin < 3,  isFull  =   false;
    if nargin < 2,  rho     =   [];         end, end
    
    % Allocate/reinitialized density object
    if isempty(rho)
        rho             =   QMol_DFT_density('isSpinPol',obj.QM.isSpinPol);
    else
        rho.clear;
        rho.set('isSpinPol',obj.QM.isSpinPol);
    end

    % Allocate density discretization
    if isFull
        if rho.isSpinPol,   rho.rhoUp   =   zeros(numel(obj.x),1);          %#ok<ALIGN> 
                            rho.rhoDw   =   zeros(numel(obj.x),1);
        else,               rho.rho     =   zeros(numel(obj.x),1);          end
    end

    % Initialize density
    rho.initialize(obj);
    
end