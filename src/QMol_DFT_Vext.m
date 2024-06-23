classdef QMol_DFT_Vext < QMol_Vmol
%QMol_DFT_Vext implementation of the external potential, for density-
%   functional theory (DFT) based computations, discretized on a Cartesian
%   grid. It
%   > overloads the QMol_Vmol class 
%   > provides support for user-defined potential (can be mixed with a list
%     of atomic center components)
%     * User-defined potential
%       Vext    function handle or domain-compatible discretization vector
%               AKA: externalPotential
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_Vmol})
function showInfo
    fprintf('  * QMol_DFT_Vext:\n      > External potential\n'); 
    QMol_DFT_Vext.version;
end
end
methods (Access=?QMol_Vmol)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the DFT model
    
    % Initialization
    ref                 =   {};

    % No external potential?
    if isempty(obj.V)   &&   obj.nbA == 0
        fprintf('  * Empty external potential\n')
        return
    end

    % Component info
    fprintf('  * External-potential functional\n');
    obj.version;

    % User-defined potential
    if ~isempty(obj.Vext)
        fprintf('  * User-defined component                                     (free form)\n');
        fprintf('    Potential = '),                            if isa(obj.Vext,'function_handle')
            fprintf('    function handle\n');                   else
            fprintf('    user-supplied discretization\n');      end
        if ~isempty(obj.DVext), fprintf('    Gradient  = ');    if isa(obj.DVext,'function_handle')
            fprintf('    function handle\n');                   else
            fprintf('    user-supplied discretization\n');      end
        end
    end
end
end
methods (Access=public)
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components (do not trust V, or DV)
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.DFT.disc.x),'real');

    if opt
        QMol_DFT_profiler.showMemoryFootprint('External functional',0,1);
        QMol_DFT_profiler.showMemoryFootprint('potential'         ,mem,2);
        QMol_DFT_profiler.showMemoryFootprint('potential gradient',mem,2);
    end

    % Return memory footprint
    mem                 =   2*mem;

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % User-defined external potential
    Vext
    DVext
end
properties (GetAccess=public,SetAccess=?QMol_suite)
    % Miscellaneous
    diffDx              =   1e-5    % For finite-difference derivative
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    externalPotential               % Vext
    externalPotentialDerivative     % DVext
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_DFT_Vext,?QMol_DFT})
    % Linked objects
    DFT                             % DFT.disc always defines the discretization domain
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Run-time variables
    V
    DV
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % Vext ~~~~~~~~~~~~~
    function set.externalPotential(obj,val),                obj.Vext    =   val;        end
    function val = get.externalPotential(obj),              val         =   obj.Vext;   end
    % DVext ~~~~~~~~~~~~
    function set.externalPotentialDerivative(obj,val),      obj.DVext   =   val;        end
    function val = get.externalPotentialDerivative(obj),    val         =   obj.DVext;  end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.DFT             =   [];
    obj.V               =   [];
    obj.DV              =   [];
    
    % Reset atomic centers
    reset@QMol_Vmol(obj);                   % Parent class updates isInit
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.diffDx),obj.diffDx   =   1e-5;       end
end
function initialize(obj,DFT)
%initialize initializes the object
    
    % Initialize atomic-center components
    initialize@QMol_Vmol(obj,DFT);          % Parent class updates isInit

    % Set links (parent class reset the object)
    obj.DFT             =   DFT;

    % Discretization of the total external potential
    x                   =   obj.DFT.disc.x(:);

    if isa(obj.Vext,'function_handle')
        obj.V           =   obj.Vext(x);
    elseif ~isempty(obj.Vext)
        obj.V           =   obj.Vext(:);
    else
        obj.V           =   zeros(size(x));
    end
    
    for k = 1:obj.nbA
        obj.V           =   obj.V + obj.atom{k}.getPotential(x);
    end

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    [~,PropNames]       =   propertyNames@QMol_Vmol;

    ClassName           =   'QMol_DFT_Vext';
    PropNames           =   [PropNames,{'Vext','DVext','diffDx',...
                                'externalPotential','externalPotentialDerivative'}];
end
end
%% Get external potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods(Access=public)
function V = getPotential(obj,V,isAdd)
%getPotential returns the discretization of the external potential
%   associated with the member properties. Optionally, provide the
%   potential object(s), where the discretization should be stored, and if
%   the external potential should be added to it. 
    
    % Initialization
    if nargin < 3,  isAdd   =   false;
    if nargin < 2,  V       =   [];         end, end

    if isempty(V)   ||   ~isAdd  
        V                   =   obj.DFT.disc.DFT_allocatePotential(V);      % If any potential, reset it
    end
    
    % Return external potential
    V.add(obj.V);
end
function DV = getPotentialDerivative(obj,dim,DV,isAdd)
%getPotentialDerivative returns the discretization of the derivative of the
%   external potential associated with the member properties. Optionally,
%   provide the potential derivative object, where the discretization 
%   should be stored, and if the external potential should be added to it.
    
    % Initialization
    if nargin < 4,  isAdd   =   false;
    if nargin < 3,  DV      =   [];         end, end

    if isempty(DV)   ||   ~isAdd  
        DV              =   obj.DFT.disc.DFT_allocatePotentialGradient(DV); % If any potential, reset it
    end

    if isempty(obj.DV),         obj.setDerivative;              end         % Compute the potential gradient (if needed)

    % Return the derivative
    switch dim
        case 1
            DV.add(1,obj.DV);
        otherwise
            error('QMol:QMol_DFT_Vext:getPotentialDerivative', ...
                 ['Unexpected dimension (' num2str(dim) ') for external-potential derivative computation.']);
    end
    
end
function E = getEnergy(obj,rho)
%getEnergy returns the external energy associated with the input density
%   object. Empty or missing density uses the one-body-density of the
%   member DFT object (discouraged).
    
    % Initialization
    if nargin < 2,      rho     =   [];                                     end
    if isempty(rho)    
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);
        rho             =   obj.DFT.rho;
    end

    % Compute energy
    if obj.DFT.isSpinPol
        E               =  [sum(obj.V.*rho.rhoUp),sum(obj.V.*rho.rhoDw)]*obj.DFT.disc.dx;
    else
        E               =   sum(obj.V.*rho.rho)*obj.DFT.disc.dx;
    end
end
end
methods (Access=private)
function setDerivative(obj)
%initializeDerivative initializes the potential-derivative component(s)
    
    % Initialization
    if ~obj.isInit
        warning('QMol:QMol_DFT_Vext:setDerivative', ...
            ['Potential derivative can only be initialized after the object itself is initialized.\n' ...
            '    No potential-derivative initialization performed.'])
        return
    end

    x                   =   obj.DFT.disc.x(:);

    % No derivative provided => finite difference approximation
    if isempty(obj.DVext)
        if isa(obj.Vext,'function_handle')
            obj.DV      =   (obj.Vext(x+.5*obj.diffDx) - obj.Vext(x-.5*obj.diffDx))/obj.diffDx;
        elseif ~isempty(obj.Vext)
            obj.DV      =   NaN(size(x));
            obj.DV(2:end-1)=.5*(obj.Vext(3:end)-obj.Vext(1:end-2))/obj.DFT.disc.dx;
            obj.DV(1)   =   (obj.Vext(2)-obj.Vext(1))/obj.DFT.disc.dx;
            obj.DV(end) =   (obj.Vext(end)-obj.Vext(end-1))/obj.DFT.disc.dx;
        else
            obj.DV      =   0;
        end
    % Derivative provided
    else
        if isa(obj.DVext,'function_handle')
            obj.DV      =   obj.DVext(x);
        else
            obj.DV      =   obj.DVext(:);
        end
    end

    % Atomic centers
    for k = 1:obj.nbA
        obj.DV          =   obj.DV + obj.atom{k}.getPotentialDerivative(1,x);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

