classdef QMol_SE_V < QMol_Vmol
%QMol_SE_V implementation of the potential, for Schrodinger-equation (SE) 
%   based computations, discretized on a Cartesian grid. It
%   > overloads the QMol_Vmol class 
%   > provides support for user-defined potential (can be mixed with a list
%     of atomic center components)
%     * User-defined potential
%       Vin     function handle or domain-compatible discretization vector
%               AKA: inputPotential
    
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
    fprintf('  * QMol_SE_V:\n      > Potential\n'); 
    QMol_SE_V.version;
end
end
methods (Access=?QMol_Vmol)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the SE model
    
    % Initialization
    ref                 =   {};

    % No external potential?
    if isempty(obj.V)   &&   obj.nbA == 0
        fprintf('  * Empty potential\n')
        return
    end

    % Component info
    fprintf('  * Potential\n');
    obj.version;

    % User-defined potential
    if ~isempty(obj.Vin)
        fprintf('  * User-defined component                                     (free form)\n');
        fprintf('    Potential = '),                            if isa(obj.Vin,'function_handle')
            fprintf('    function handle\n');                   else
            fprintf('    user-supplied discretization\n');      end
        if ~isempty(obj.DVin), fprintf('    Gradient  = ');    if isa(obj.DVin,'function_handle')
            fprintf('    function handle\n');                   else
            fprintf('    user-supplied discretization\n');      end
        end
    end
end
end
methods (Access=public)
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the SE object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components (do not trust V, or DV)
    mem                 =   QMol_SE_profiler.getMemoryFootprint(numel(obj.SE.disc.x),'real');   if obj.isBasis
    mem                 =   mem + QMol_SE_profiler.getMemoryFootprint(size(obj.SE.disc.v,2),'real'); end

    if opt
        QMol_SE_profiler.showMemoryFootprint('Potential'         ,0,1);
        QMol_SE_profiler.showMemoryFootprint('potential'         ,mem,2);
        QMol_SE_profiler.showMemoryFootprint('potential gradient',mem,2);
    end

    % Return memory footprint
    mem                 =   2*mem;

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % User-defined external potential
    Vin
    DVin
end
properties (GetAccess=public,SetAccess=?QMol_suite)
    % Miscellaneous
    diffDx              =   1e-5    % For finite-difference derivative
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    inputPotential                  % Vin
    inputPotentialDerivative        % DVin
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_SE_V,?QMol_SEq})
    % Linked objects
    SE                              % SE.disc always defines the discretization domain
    % Discretization
    isBasis
end
properties (Dependent,GetAccess=public,SetAccess=private)
    potential                       % V
    potentialDerivative             % DV
    potentialMatrix                 % mV
    potentialDerivativeMatrix       % mDV
end
properties (Transient,Hidden,GetAccess=public,SetAccess=private)
    % Run-time variables
    V
    DV
    mV
    mDV
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % Vin ~~~~~~~~~~~~~~
    function set.inputPotential(obj,val),               obj.Vin     =   val;        end
    function val = get.inputPotential(obj),             val         =   obj.Vin;    end
    % DVin ~~~~~~~~~~~~~
    function set.inputPotentialDerivative(obj,val),     obj.DVin    =   val;        end
    function val = get.inputPotentialDerivative(obj),   val         =   obj.DVin;   end
    % V ~~~~~~~~~~~~~~~~
    function val = get.potential(obj),                  val         =   obj.V;      end
    % DV ~~~~~~~~~~~~~~~
    function val = get.potentialDerivative(obj),        val         =   obj.DV;     end
    % mV ~~~~~~~~~~~~~~~
    function val = get.potentialMatrix(obj),            val         =   obj.mV;     end
    % mDV ~~~~~~~~~~~~~~
    function val = get.potentialDerivativeMatrix(obj),  val         =   obj.mDV;    end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.SE              =   [];
    obj.isBasis         =   [];
    obj.V               =   [];
    obj.DV              =   [];
    obj.mV              =   [];
    obj.mDV             =   [];
    
    % Reset atomic centers
    reset@QMol_Vmol(obj);                   % Parent class updates isInit
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.diffDx),obj.diffDT   =   1e-5;       end
end
function initialize(obj,SE)
%initialize initializes the object
    
    % Initialize atomic-center components
    initialize@QMol_Vmol(obj,SE);           % Parent class updates isInit

    % Set links (parent class reset the object)
    obj.SE              =   SE;
    obj.isBasis         =   SE.disc.isBasis;

    % Discretization of the total external potential
    x                   =   obj.SE.disc.x(:);

    if isa(obj.Vin,'function_handle')
        obj.V           =   obj.Vin(x);
    elseif ~isempty(obj.Vin)
        obj.V           =   obj.Vin(:);
    else
        obj.V           =   zeros(size(x));
    end
    
    for k = 1:obj.nbA
        obj.V           =   obj.V + obj.atom{k}.getPotential(x);
    end
                                                                            if obj.isBasis
    % Potential matrix
    obj.mV              =   NaN(SE.disc.nV,SE.disc.nV);

    if isreal(SE.disc.v), for k = 1:SE.disc.nV, for l = 1:k                 %#ok<ALIGN> 
        obj.mV(k,l)     =   sum(     SE.disc.v(:,k) .* obj.V.*SE.disc.v(:,l)) * SE.disc.dx;    obj.mV(l,k) =   obj.mV(k,l);
    end, end, else ,      for k = 1:SE.disc.nV, for l = 1:k                 %#ok<ALIGN> 
        obj.mV(k,l)     =   sum(conj(SE.disc.v(:,k)).* obj.V.*SE.disc.v(:,l)) * SE.disc.dx;    obj.mV(l,k) =   conj(obj.mV(k,l));
    end, end, end
                                                                            end
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    [~,PropNames]       =   propertyNames@QMol_Vmol;

    ClassName           =   'QMol_SE_V';
    PropNames           =   [PropNames,{'Vin','DVin','diffDx',...
                                'inputPotential','inputPotentialDerivative'}];
end
end
%% Get external potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods(Access=public)
function E = getEnergy(obj,wfcn)
%getEnergy returns the energy associated with the input wave function
%   object. Empty or missing density uses the wave function of the member 
%   SE object (discouraged).
    
    % Initialization
    if nargin < 2,      wfcn    =   [];                                     end
    if isempty(wfcn),   wfcn    =   obj.SE.wfcn;                            end

    % Compute energy
    E                   =   0;
    if ~obj.isBasis,                                                        for k = 1:size(wfcn.wfcn,2)
        E               =   E + sum(obj.V.*abs(wfcn.wfcn(:,k)).^2);         end
        E               =   E*obj.SE.disc.dx;
    else,                                                                   for k = 1:size(wfcn.wfcn,2)
        E               =   E + wfcn.wfcn(:,k)'*obj.mV*wfcn.wfcn(:,k);      end
        E               =   real(E);                                        % should aleardy be real
    end

end
end
methods (Access=public)
function setDerivative(obj)
%initializeDerivative initializes the potential-derivative component(s)
    
    % Initialization
    if ~obj.isInit
        warning('QMol:QMol_SE_V:setDerivative', ...
            ['Potential derivative can only be initialized after the object itself is initialized.\n' ...
            '    No potential-derivative initialization performed.'])
        return
    end

    x                   =   obj.SE.disc.x(:);

    % No derivative provided => finite difference approximation
    if isempty(obj.DVin)
        if isa(obj.Vin,'function_handle')
            obj.DV      =   (obj.Vin(x+.5*obj.diffDx) - obj.Vin(x-.5*obj.diffDx))/obj.diffDx;
        elseif ~isempty(obj.Vin)
            obj.DV      =   NaN(size(x));
            obj.DV(2:end-1)=.5*(obj.Vin(3:end)-obj.Vin(1:end-2))/obj.SE.disc.dx;
            obj.DV(1)   =   (obj.Vin(2)-obj.Vin(1))/obj.SE.disc.dx;
            obj.DV(end) =   (obj.Vin(end)-obj.Vin(end-1))/obj.SE.disc.dx;
        else
            obj.DV      =   0;
        end
    % Derivative provided
    else
        if isa(obj.DVin,'function_handle')
            obj.DV      =   obj.DVin(x);
        else
            obj.DV      =   obj.DVin(:);
        end
    end

    % Atomic centers
    for k = 1:obj.nbA
        obj.DV          =   obj.DV + obj.atom{k}.getPotentialDerivative(1,x);
    end
                                                                            if obj.isBasis
    % Matrix representation of the potential
    obj.mDV             =   NaN(obj.SE.disc.nV,obj.SE.disc.nV);

    if isreal(obj.SE.disc.v), for k = 1:obj.SE.disc.nV, for l = 1:k         %#ok<ALIGN> 
        obj.mDV(k,l)    =   sum(     obj.SE.disc.v(:,k) .* obj.DV.*obj.SE.disc.v(:,l)) * obj.SE.disc.dx;    obj.mDV(l,k) =   obj.mDV(k,l);
    end, end, else ,      for k = 1:obj.SE.disc.nV, for l = 1:k             %#ok<ALIGN> 
        obj.mDV(k,l)    =   sum(conj(obj.SE.disc.v(:,k)).* obj.DV.*obj.SE.disc.v(:,l)) * obj.SE.disc.dx;    obj.mDV(l,k) =   conj(obj.mDV(k,l));
    end, end, end
                                                                            end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

