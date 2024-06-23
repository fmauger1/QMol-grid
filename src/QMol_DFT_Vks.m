classdef QMol_DFT_Vks < QMol_suite
%QMol_DFT_Vks Kohn-Sham (KS) potential operator discretized on a Cartesian-
%   grid domain. It defines
%   > Common interface for orbitals-dependent and independent KS potentials
%   > The potential discretization itself
%    
%   NOTES (for developers)
%   > The DFT object is accessible through the (transient) disc property,
%     if needed. This is for consistency with KS orbitals and density
%     objects.
%   > Only the explicit (with respect to the one-body density) part of the
%     KS potential is discretized. Implicit components (e.g., exact
%     exchange) only exist as a (linear) operator and cannot be
%     discretized on the grid.
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT_Vks})
function showInfo
    fprintf('  * QMol_DFT_Vks:\n      > Kohn-Sham potential\n'); 
    QMol_DFT_Vks.version;
end
end
methods (Access=public)
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components (do not trust V, Vup, or Vdw, they might
    % have been initialized empty; fetch size from domain and parent DFT)
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.disc.x),'real');

    if opt,                                                                     if obj.isSpinPol
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham potential',2*mem,1);   else
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham potential',  mem,1);   end
    end

    % Return memory footprint
    if obj.isSpinPol,   mem     =   2*mem;      end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Transient,Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Linked objects
    disc                            % discretization object (for possible future development)
end
properties (GetAccess=public,SetAccess=?QMol_suite)
    isSpinPol
end
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    V
    Vup
    Vdw
    Vimp                =   {}
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    potential                       % V
    potentialUp                     % Vup
    potentialDown                   % Vdown
    potentialImplicit               % Vimp
end
methods (Access=public,Static)
    function b = isBasis,   b   =   false;  end
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % V ~~~~~~~~~~~~~~~~
    function val = get.potential(obj),              val         =   obj.V;      end
    function set.potential(obj,val),                obj.V       =   val;        end
    % Vup ~~~~~~~~~~~~~~
    function val = get.potentialUp(obj),            val         =   obj.Vup;    end
    function set.potentialUp(obj,val),              obj.Vup     =   val;        end
    % Vdw ~~~~~~~~~~~~~~
    function val = get.potentialDown(obj),          val         =   obj.Vdw;    end
    function set.potentialDown(obj,val),            obj.Vdw     =   val;        end
    % Vimp ~~~~~~~~~~~~~
    function val = get.potentialImplicit(obj),      val         =   obj.Vimp;   end
    function set.potentialImplicit(obj,val),        obj.Vimp    =   val;        end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function clear(obj,varargin)
%clear resets all or selected member properties
    
    % streamlined clear all
    if nargin == 1
        % Clear all properties
        obj.reset;

        obj.disc        =   [];
        obj.isSpinPol   =   [];
        obj.V           =   [];
        obj.Vup         =   [];
        obj.Vdw         =   [];
        obj.Vimp        =   {};

    % parent-class selected clear
    else
        clear@QMol_suite(obj,varargin{:});
    end
end
function initialize(obj,disc)
%initialize initializes the object
    
    % Reinitialize any component
    obj.reset;

    % Update discretization
    if nargin == 1,     obj.disc    =   [];
    else,               obj.disc    =   disc;                               end
    
    % Nothing to do -- update isInit
    obj.isInit          =   true;
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_DFT_Vks';
    PropNames           =  {'isSpinPol','V','Vup','Vdw','Vimp',...
                            'potential','potentialUp','potentialDown','potentialImplicit'};
end
end
%% Arithmetic on potentials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function add(obj,varargin)
%add adds the input potential component to the object
    
    % Add potential component
    if ~isnumeric(varargin{1})
        % Implicit potential
        if isempty(obj.Vimp),   obj.Vimp=  varargin(1);
        else,                   obj.Vimp=  [obj.Vimp varargin(1)];          end

    elseif obj.isSpinPol
        % Spin-polarized DFT
        if nargin == 2
            % add(obj,V), with same V to add to spin up and down components
            if isempty(obj.Vup),obj.Vup =   varargin{1};                    %#ok<ALIGN> 
                                obj.Vdw =   varargin{1};
            else,               obj.Vup =   obj.Vup + varargin{1};
                                obj.Vdw =   obj.Vdw + varargin{1};          end
        elseif nargin == 3
            % add(obj,V_up,V_down), with different spin up and down components
            if isempty(obj.Vup),obj.Vup =   varargin{1};                    %#ok<ALIGN> 
                                obj.Vdw =   varargin{2};
            else,               obj.Vup =   obj.Vup + varargin{1};
                                obj.Vdw =   obj.Vdw + varargin{2};          end
        else
            % unexpected case, throw warning 
            warning('QMolGrid:DFT_Vks:add',...
                   ['Unknown addition type for spin-polarized DFT-potential discretization.\n' ...
                    'No potential component added.']);
        end
    else
        % Spin-restricted DFT
        if nargin == 2
            % add(obj,V)
            if isempty(obj.V),  obj.V   =   varargin{1};
            else,               obj.V   =   obj.V + varargin{1};    end
        else
            % unexpected case, throw warning
            warning('QMolGrid:DFT_Vks:add',...
                   ['Unknown addition type for spin-restricted DFT-potential discretization.\n' ...
                    'No potential component added.']);
        end
    end

    % Adding a component undo any initialization
    obj.isInit          =   false;
end
function Hp = applyPotential(obj,p,isUp)
%applyPotential potential operator V | psi >
%   The potential is assumed to have been properly initialized, but no
%   check for that is performed (try applying it no matter what).
    
    % Explicit part
    if obj.isSpinPol                                                        %#ok<ALIGN> 
        if isUp,    if isempty(obj.Vup),    Hp  =   0*p;    else,   Hp  =   obj.Vup .* p;   end
        else,       if isempty(obj.Vdw),    Hp  =   0*p;    else,   Hp  =   obj.Vdw .* p;   end, end
    else,           if isempty(obj.V),      Hp  =   0*p;    else,   Hp  =   obj.V   .* p;   end, end

    % Implicit part
    for k = 1:numel(obj.Vimp)
        if obj.isSpinPol,       Hp  =   Hp + obj.Vimp{k}(p,isUp);
        else,                   Hp  =   Hp + obj.Vimp{k}(p);                end
    end

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

