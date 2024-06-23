classdef QMol_DFT_Vks_grad < QMol_suite
%QMol_DFT_Vks_grad gradient of the Kohn-Sham (KS) potential operator 
%   discretized on a Cartesian-grid domain. It defines
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
methods (Static,Access={?QMol_doc,?QMol_DFT_Vks_grad})
function showInfo
    fprintf('  * QMol_DFT_Vks_grad:\n      > Kohn-Sham potential gradient\n'); 
    QMol_DFT_Vks_grad.version;
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

    if opt,                                                                             if obj.isSpinPol
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham potential gradient',2*mem,1);  else
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham potential gradient',  mem,1);  end
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
    DV
    DVup
    DVdw
    DVimp               =   {};
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    potentialGradient               % DV
    potentialGradientUp             % DVup
    potentialGradientDown           % DVdown
    potentialGradientImplicit       % DVimp
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % V ~~~~~~~~~~~~~~~~
    function val = get.potentialGradient(obj),          val         =   obj.DV;     end
    function set.potentialGradient(obj,val),            obj.DV      =   val;        end
    % Vup ~~~~~~~~~~~~~~
    function val = get.potentialGradientUp(obj),        val         =   obj.DVup;   end
    function set.potentialGradientUp(obj,val),          obj.DVup    =   val;        end
    % Vdw ~~~~~~~~~~~~~~
    function val = get.potentialGradientDown(obj),      val         =   obj.DVdw;   end
    function set.potentialGradientDown(obj,val),        obj.DVdw    =   val;        end
    % Vimp ~~~~~~~~~~~~~
    function val = get.potentialGradientImplicit(obj),  val         =   obj.DVimp;  end
    function set.potentialGradientImplicit(obj,val),    obj.DVimp   =   val;        end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % streamlined clear all
    if nargin == 1
        % Clear all properties
        obj.reset;

        obj.disc        =   [];
        obj.isSpinPol   =   [];
        obj.DV          =   [];
        obj.DVup        =   [];
        obj.DVdw        =   [];
        obj.DVimp       =   {};

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
    
    ClassName           =   'QMol_DFT_Vks_grad';
    PropNames           =  {'isSpinPol','DV','DVup','DVdw','DVimp',...
                            'potentialGradient','potentialGradientUp',...
                            'potentialGradientDown','potentialGradientImplicit'};
end
end
%% Arithmetic on potentials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function add(obj,dim,varargin)
%add adds the input potential component to the object

    % Initialization
    if numel(dim) ~= 1   ||   dim ~= 1
        error('QMolGrid:DFT_Vks_grad:addDim',...
            ['Unexpected dimension (' num2str(dim) ') when adding a component to the potential gradient.'])
    end
    
    % Add potential component
    if ~isnumeric(varargin{1})
        % Implicit potential
        if isempty(obj.DVimp),      obj.DVimp   =  varargin(1);
        else,                       obj.DVimp   =  [obj.DVimp varargin(1)];     end

    elseif obj.isSpinPol
        % Spin-polarized DFT
        if nargin == 3
            % add(obj,1,V), with same V to add to spin up and down components
            if isempty(obj.DVup),   obj.DVup    =   varargin{1};            %#ok<ALIGN> 
                                    obj.DVdw    =   varargin{1};
            else,                   obj.DVup    =   obj.DVup + varargin{1};
                                    obj.DVdw    =   obj.DVdw + varargin{1}; end
        elseif nargin == 4
            % add(obj,1,V_up,V_down), with different spin up and down components
            if isempty(obj.DVup),   obj.DVup    =   varargin{1};            %#ok<ALIGN> 
                                    obj.DVdw    =   varargin{2};
            else,                   obj.DVup    =   obj.DVup + varargin{1};
                                    obj.DVdw    =   obj.DVdw + varargin{2}; end
        else
            % unexpected case, throw warning 
            warning('QMolGrid:DFT_Vks_grad:add',...
                   ['Unknown addition type for spin-polarized DFT-potential gradient discretization.\n' ...
                    'No potential gradient component added.']);
        end
    else
        % Spin-restricted DFT
        if nargin == 3
            % add(obj,1,V)
            if isempty(obj.DV),     obj.DV      =   varargin{1};
            else,                   obj.DV      =   obj.DV + varargin{1};   end
        else
            % unexpected case, throw warning
            warning('QMolGrid:DFT_Vks_grad:add',...
                   ['Unknown addition type for spin-restricted DFT-potential gradient discretization.\n' ...
                    'No potential gradient component added.']);
        end
    end

    % Adding a component undo any initialization
    obj.isInit          =   false;
end
function Hp = applyPotentialGradient(obj,dim,p,isUp)
%applyPotential potential operator V | psi >
%   The potential is assumed to have been properly initialized, but no
%   check for that is performed (try applying it no matter what).
    
    % Initialization
    if numel(dim) ~= 1   ||   dim ~= 1
        error('QMolGrid:DFT_Vks_grad:applyPotentialGradientDim',...
            ['Unexpected dimension (' num2str(dim) ') when applying the potential gradient.'])
    end

    % Explicit part
    if obj.isSpinPol                                                        %#ok<ALIGN> 
        if isUp,    if isempty(obj.DVup),   Hp  =   0*p;    else,   Hp  =   obj.DVup .* p;  end
        else,       if isempty(obj.DVdw),   Hp  =   0*p;    else,   Hp  =   obj.DVdw .* p;   end, end
    else,           if isempty(obj.DV),     Hp  =   0*p;    else,   Hp  =   obj.DV   .* p;   end, end

    % Implicit part
    for k = 1:numel(obj.DVimp)
        if obj.isSpinPol,       Hp  =   Hp + obj.DVimp{k}(p,isUp);
        else,                   Hp  =   Hp + obj.DVimp{k}(p);            end
    end

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

