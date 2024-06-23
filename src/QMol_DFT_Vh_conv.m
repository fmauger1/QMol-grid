classdef QMol_DFT_Vh_conv < QMol_suite
%QMol_DFT_Vh_conv implementation of the Hartree potential functional for
%   density-functional theory (DFT) simulations discretized on a Cartesian
%   grid. Hartree-functional components are computed using an explicit
%   convolution scheme defined over an extended domain. It is almost always
%   slower than the fast-Fourier transform implementation but does not
%   require any tuning of the extended domain.
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT_Vh_conv})
function showInfo
    fprintf( '  * QMol_DFT_Vh_conv:\n');
    fprintf(['      > Hartree potential\n',...
             '      > Explicit-convolution scheme\n']); 
    QMol_DFT_Vh_conv.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Initialization
    ref                 =   {};

    % Header
    fprintf('  * Hartree-potential functional                      explicit convolution\n');

    % SIC (if any)
    if obj.SIC == 1
        fprintf('    with average-density self-interaction correct. (ADSIC) [Legrand 2002].\n');
        ref             =   [ref, {'Legrand 2002'}];
    end

    % Electron-interaction potential
    fprintf('    Interaction pot. = %s (elec.-elec.)\n',func2str(obj.Vee));

    % Version
    obj.version;

end
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components (do not trust V, or DV)
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(2*numel(obj.DFT.disc.x)-1,'real');

    if opt
        QMol_DFT_profiler.showMemoryFootprint('Hartree functional',mem,1);
    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Electron-electron interaction potential
    Vee                 =   @(x) 1./sqrt(x.^2+2)
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    interactionPotential            % Vee
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_DFT_Vh_conv,?QMol_DFT})
    % Linked objects
    DFT                             % DFT.disc always defines the domain
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Run-time variables
    V
    SIC                             % 0 = none, 1 = ADSIC
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % Vee ~~~~~~~~~~~~~~
    function set.interactionPotential(obj,val),     obj.Vee =   val;        end
    function val = get.interactionPotential(obj),   val     =   obj.Vee;    end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.DFT             =   [];
    obj.V               =   [];
    obj.SIC             =   [];
    
    % Initialization status
    obj.isInit          =   false;
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.Vee),    obj.Vee     =   @(x) 1./sqrt(x.^2+2);           end
end
function initialize(obj,DFT,SIC)
%initialize initializes the object

    % Initialization
    obj.reset;
    
    % Set links
    obj.DFT             =   DFT;

    % Identify the flavor of SIC
    if nargin == 2,     SIC = [];       end
    if isempty(SIC),    SIC = 'none';   end

    switch lower(SIC)
        case {'none','off'}
            obj.SIC     =   0;
        case {'adsic','average density','averagedensity','average_density'}
            obj.SIC     =   1;
        otherwise
            warning('QMol:QMol_DFT_Vh_conv:SIC', ...
                    ['Unknown flavor of self-interaction correction (SIC) ' SIC '\n; ' ...
                     'SIC ignored for the Hartree potential.'])
            obj.SIC     =   0;
    end
    
    % Extended potential
    obj.V               =   obj.Vee(...
                                (0:2*numel(obj.DFT.disc.x)-2).' * obj.DFT.disc.dx ...
                               -(obj.DFT.disc.x(end)-obj.DFT.disc.x(1)) );
    
    % Miscellaneous
    obj.isInit          =   true;

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_DFT_Vh_conv';
    PropNames           =   {'Vee','interactionPotential'};
end
end
%% Get Hartree potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function V = getPotential(obj,rho,V,isAdd)
%getPotential returns the discretization of the Hartree potential
%   associated with the member properties and input one-body density (rho).
%   Optionally, provide the potential object(s), where the discretization
%   should be stored, and if the Hartree potential should be added to it.
    
    % Initialization
    if nargin < 4,  isAdd   =   false;                                      %#ok<ALIGN> 
    if nargin < 3,  V       =   [];
    if nargin < 2,  rho     =   [];             end, end, end
    
    if isempty(V)   ||   ~isAdd  
        V               =   obj.DFT.disc.DFT_allocatePotential(V);          % If any potential, reset it
    end
    
    if isempty(rho)
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);                % Use DFT one-body density storage
        rho             =   obj.DFT.rho;
    end

    % Compute Hartree potential
    if rho.isSpinPol
        Vh              =   conv(rho.rhoUp+rho.rhoDw,obj.V,'same') * obj.DFT.disc.dx;
    else
        Vh              =   conv(rho.rho            ,obj.V,'same') * obj.DFT.disc.dx;
    end
    
    % Self-interaction correction (SIC == 0 for no correction)
    if obj.SIC == 1
        % ADSIC
        N               =   sum(obj.DFT.Ntot);
        Vh              =   ( (N-1)/N ) * Vh;
    end

    % Return Hartree potential
    V.add(Vh);

end
function DV = getPotentialDerivative(obj,dim,rho,DV,isAdd)
%getPotentialDerivative returns the discretization of the derivative of the
%   Hartree potential associated with the input one-body density (rho).
%   Optionally, provide the potential object(s), where the discretization
%   should be stored, and if the exchange potential should be added to it.
%   The potential derivative is simply the Hartree potential of the
%   derivative of the density.

    % Initialization
    if nargin < 5,  isAdd   =   false;                                      %#ok<ALIGN> 
    if nargin < 4,  DV      =   [];
    if nargin < 3,  rho     =   [];     end, end, end
    
    if isempty(DV)   ||   ~isAdd  
        DV              =   obj.DFT.disc.DFT_allocatePotentialGradient(DV); % If any potential, reset it
    end
    
    if isempty(rho)
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);                % Use DFT one-body density storage
        rho             =   obj.DFT.rho;
    end

    % Compute Hartree potential gradient
    switch dim
        case 1
            % Gradient is Hartree potential of gradient of density
            if rho.isSpinPol
                DVh     =   conv(rho.D_rhoUp+rho.D_rhoDw,obj.V,'same') * obj.DFT.disc.dx;
            else
                DVh     =   conv(rho.D_rho              ,obj.V,'same') * obj.DFT.disc.dx;
            end

        otherwise
            error('QMol:QMol_DFT_Vh_conv:getPotentialDerivative', ...
                 ['Unexpected dimension (' num2str(dim) ') for Hartree-potential derivative computation.']);
    end
    
    % Self-interaction correction (SIC == 0 for no correction)
    if obj.SIC == 1
        % ADSIC
        N               =   sum(obj.DFT.Ntot);
        DVh             =   ( (N-1)/N ) * DVh;
    end

    % Return Hartree potential
    DV.add(dim,DVh);

end
function E = getEnergy(obj,rho)
%getEnergy returns the Hartree energy associated with the member properties
%   and input one-body density (rho). Empty or missing density uses the 
%   one-body-density of the member DFT object (discouraged).
    
    % Initialization
    if nargin < 2,      rho     =   [];                                     end
    if isempty(rho),    rho     =   obj.DFT.getDensity(obj.DFT.rho);        end

    % Compute Hartree energy
    if rho.isSpinPol
        E               =   0.5*sum((rho.rhoUp+rho.rhoDw).*conv(rho.rhoUp+rho.rhoDw,obj.V,'same')) * obj.DFT.disc.dx^2;
    else
        E               =   0.5*sum( rho.rho             .*conv(rho.rho            ,obj.V,'same')) * obj.DFT.disc.dx^2;
    end
    
    % Self-interaction correction (SIC == 0 for no correction)
    if obj.SIC == 1
        % ADSIC
        N               =   sum(obj.DFT.Ntot);
        E               =   ( (N-1)/N ) * E;
    end

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

