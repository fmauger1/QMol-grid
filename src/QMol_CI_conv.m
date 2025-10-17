classdef QMol_CI_conv < QMol_CI
%QMol_CI_conv configuration interaction module
%   Use QMol_CI_conv to describe a configuration-interaction (CI) model
%   associated with a spin-restricted set of orbitals, typically obtained
%   via a density-functional theory (DFT) or Hartree-Fock (HF) ground state
%   calculation. The class also provides functionalities for calculating
%   the CI matrix and associated ground- and excited-state wave functions,
%   as well as the dipole-coupling matrix, defined as < chi | r | chi' >,
%   where | chi > and | chi' > are two configuration states and r is the
%   position vector.
%
%   Internally, the class uses an explicit convolution to calculate the
%   two-electron integrals (QMol_CI_fft provides a fast-Fourier alternative
%   for faster calculation, but with additional parameters to set).
%
%   CI = QMol_CI_conv('name1','value1',___) creates a CI object with the
%     name properties set to the specified values.
%
%   CI = QMol_CI_conv(DFT,'name1','value1',___), where DFT is a spin-
%     restricted DFT object, reads off the domain, number of electrons, 
%     orbital basis, atomic or molecular potential, and electron-electron
%     interaction potential off DFT. Additional or overwriting parameters
%     can be specified as a list of name-value pair arguments.
%
%     For basis-set DFT models, the orbitals are first reconstructed over
%     the domain discretization before being stored in the CI object.
% 
%   Editable properties:
%     * Domain discretization: discretization, xspan
%     * Electronic structure: numberElectron, orbitalBasis, 
%         configurationBasis, waveFunction
%     * Potential: externalPotential, interactionPotential
%     * CI calculation: display, tolerance, CImatrix, DXmatrix, isDipole
%     * Configuration basis: type, reference, active, frozen, noDouble,
%         noEmpty
%
%   Uneditable properties:
%     * isInitialized, isSetConfigBasis, isSpinPol, dim
%
%   Methods:
%     * Changing class properties: set, reset, clear
%     * Initialization: initialize, setConfigurationBasis, 
%         cleanConfigurationBasis
%     * Run-time documentation: showDocumentation, getMemoryProfile
%     * CI calculation: computeCImatrix, sparsify, computeGroundState, 
%         getEnergy, showEnergy, analyzeCIspectrum
%     * Density: getDensityMatrix, getDensity
%
%     > The CI object must be initialized before using any of the 
%       computeCImatrix, computeGroundState, getEnergy, showEnergy, and
%       analyzeCIspectrum (see the documentation for further requirement
%       for the CI object to hold valid CImatrix and waveFunction).
%
%   See also QMol_DFT

%   Version     Date        Author
%   01.22.000   05/20/2025  F. Mauger
%       Creation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.22.000','05/20/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_CI})
function showInfo
    fprintf('  * QMol_CI_conv:\n      > Explicit convolution scheme\n'); 
    QMol_CI_conv.version;
end
end
methods (Access=?QMol_CI)
function ref = showDoc(obj)
%showDoc QMol_CI_conv specific part of the run-time documentation

    % Implementation of the CI model
    fprintf('  * Explicit convolution for calculating electron interaction integrals.\n');
    obj.version;

    ref = {};

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,Transient,GetAccess=public,SetAccess=?QMol_suite)
    x                       % Implicit definition of the domain
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    xspan                   % (x) domain grid discretization (see full documentation) [ vector (default []) ]
end
properties (Transient,Access=private)
    V2e                     % internal electron-electron interaction potential
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % x ~~~~~~~~~~~~~~~~
    function val = get.x(obj),              if obj.isInit,  val         =   obj.disc.x; % obj.x should be empty after the object initialization
                                            else,           val         =   obj.x;      end,    end
    function set.xspan(obj,val),                            obj.x       =   val;        end
    function val = get.xspan(obj),                          val         =   obj.x;      end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function initialize(obj,varargin)
%initialize initializes the CI object and all its components
    
    % Implicit definition of the domain discretization
    if ~isempty(obj.x)
        % Define/update domain
        if isempty(obj.disc)
            obj.disc    =   QMol_disc('x',obj.x);
        else
            obj.disc.set('x',obj.x);
        end
        % Clear data
        obj.x           =   [];
    end

    % Common initialization
    initialize@QMol_CI(obj,varargin{:});

end
function reset(obj)
%reset resets all temporary (transient) properties of the object
    
    % Run-time variables
    obj.V2e             =   [];
    obj.isSetConfigBasis=   false;

    % Call parent clear (will update initialzation status)
    reset@QMol_CI(obj);

end
end
methods (Access=?QMol_CI)
function isOK = initializeVee(obj)
%initializeVee initializes the internal discretization of the electron-
%   interaction potential.
    
    % Interaction potential
    obj.V2e             =   obj.Vee(...
                                (0:2*numel(obj.disc.x)-2).' * obj.disc.dx ...
                               -(obj.disc.x(end)-obj.disc.x(1)) );

    % Return status
    isOK                =   true;
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    % Parent-class components
    [~,PropNames]       =   QMol_CI.propertyNames;

    ClassName           =   'QMol_CI_conv';
    PropNames           =   [PropNames,{'x','xspan'}];
end
end
%% Core and two-electron integrals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_CI)
function h = getCoreIntegral(obj,k,l)
%getCoreIntegral calculates and returns the core integral < k | h | l >

    % Kinetic operator
    if obj.isReal,  h   =   sum(     obj.SOB(:,k) .*real(ifft(obj.disc.T .* fft(obj.SOB(:,l)))));
    else,           h   =   sum(conj(obj.SOB(:,k)).*     ifft(obj.disc.T .* fft(obj.SOB(:,l))) );   end

    % External potential
    if obj.isReal,  h   =   h + sum(     obj.SOB(:,k)  .* obj.Vext.V .* obj.SOB(:,l));
    else,           h   =   h + sum(conj(obj.SOB(:,k)) .* obj.Vext.V .* obj.SOB(:,l));              end

    % Finalize
    h                   =   h * obj.disc.dx;
end
function v = getTwoElectronIntegral(obj,k,l,m,n)
% getTwoElectronIntegral calculates and returns the two-electron integral 
%   < k l | m n>

    if obj.isReal,  v   =   sum(     obj.SOB(:,k) .*obj.SOB(:,m) .* conv(     obj.SOB(:,l) .*obj.SOB(:,n),obj.V2e,'same'));
    else,           v   =   sum(conj(obj.SOB(:,k)).*obj.SOB(:,m) .* conv(conj(obj.SOB(:,l)).*obj.SOB(:,n),obj.V2e,'same')); end

    v                   =   v * obj.disc.dx^2;
end
function d = getDipoleIntegral(obj,k,l)
%getDipoleIntegral calculates and return the dipole integral < k | x | l >

    % Dipole coupling
    if obj.isReal,  d   =   sum(     obj.SOB(:,k)  .* obj.disc.x(:) .* obj.SOB(:,l)) * obj.disc.dx;
    else,           d   =   sum(conj(obj.SOB(:,k)) .* obj.disc.x(:) .* obj.SOB(:,l)) * obj.disc.dx; end

end
end
%% One-body density and one-particle density matrices %%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function rho = getDensity(obj,varargin)
%getDensity computes the one-body density
%   Use getDensity to calculate the one-body density associated with a wave
%   function or one-particle reduced density matrices (1-RDM) of a CI
%   model.
%
%   rho = obj.getDensity() returns the one-body density of the waveFunction
%   property, as QMol_DFT_density object. Note: waveFunction must only
%   contain a single wave function.
%
%   rho = obj.getDensity(psi) retuns the one-body density rho of the wave
%   function psi, specified as a column vector of the population
%   coefficients in the configuration-state basis configurationBasis. The
%   output rho is a QMol_DFT_density object.
%
%   rho = obj.getDensity(rhoUp,rhoDown) returns the one-body density of the
%   1-RDM matrices rhoUp and rhoDown, which are typically calculated using
%   the getDensityMatrix method.
%
%   obj.getDensity(___,rho) stores the one-body density in the input 
%   QMol_DFT_density object rho.

    % Initialization
    if nargin == 1
        % rho = obj.getDensity()
        [rUp,rDw]       =   obj.getDensityMatrix();
        rho             =   obj.disc.DFT_allocateDensity;
    elseif nargin == 2,                                                     if isnumeric(varargin{1})
        % rho = obj.getDensity(psi)
        [rUp,rDw]       =   obj.getDensityMatrix(varargin{1});
        rho             =   obj.disc.DFT_allocateDensity;                   else

        % rho = obj.getDensity(rho)
        [rUp,rDw]       =   obj.getDensityMatrix();
        rho             =   varargin{1};                                    end
    elseif nargin == 3,                                                     if isnumeric(varargin{2})
        % rho = obj.getDensity(rhoUp,rhoDown)
        rUp             =   varargin{1};
        rDw             =   varargin{2};
        rho             =   obj.disc.DFT_allocateDensity;                   else

        % rho = obj.getDensity(psi,rho)
        [rUp,rDw]       =   obj.getDensityMatrix(varargin{1});
        rho             =   varargin{2};                                    end
    elseif nargin == 4
        % rho = obj.getDensity(rhoUp,rhoDown,rho)
        rUp             =   varargin{1};
        rDw             =   varargin{2};
        rho             =   varargin{3};
    else
        error('QMol:CI:getDensity','Unexpected number of input arguments. Check the documentation')
    end

    % Build the density
    rho.isSpinPol       =   true;

    nMO                 =   size(obj.SOB,2);
    rho.rhoUp           =   zeros(size(obj.SOB,1),1);
    rho.rhoDw           =   zeros(size(obj.SOB,1),1);
    rho.rho             =   [];

    for k = 1:nMO
        % Diagonal elements
        if obj.isReal,  rho.rhoUp   =   rho.rhoUp + rUp(k,k)*    obj.SOB(:,k).^2;
                        rho.rhoDw   =   rho.rhoDw + rDw(k,k)*    obj.SOB(:,k).^2;
        else,           rho.rhoUp   =   rho.rhoUp + rUp(k,k)*abs(obj.SOB(:,k)).^2;
                        rho.rhoDw   =   rho.rhoDw + rDw(k,k)*abs(obj.SOB(:,k)).^2;
        end
                                                                            for l = k+1:nMO
        % Off-diagonal elements
        if obj.isReal,  rho.rhoUp   =   rho.rhoUp + 2*real(rUp(k,l))*     obj.SOB(:,k) .*obj.SOB(:,l);
                        rho.rhoDw   =   rho.rhoDw + 2*real(rDw(k,l))*     obj.SOB(:,k) .*obj.SOB(:,l);
        else,           rho.rhoUp   =   rho.rhoUp + 2*real(rUp(k,l) *conj(obj.SOB(:,k)).*obj.SOB(:,l));
                        rho.rhoDw   =   rho.rhoDw + 2*real(rDw(k,l) *conj(obj.SOB(:,k)).*obj.SOB(:,l));
        end,                                                                end
    end

    % Initialize object (making sure the density and CI disc match)
    rho.initialize(obj.disc);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

