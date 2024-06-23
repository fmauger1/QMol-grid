classdef QMol_DFT_spinPol < QMol_DFT
%QMol_DFT_spinPol implementation of spin-polarized (SP) density-functional
%   theory (DFT) discretization on a Cartesian grid. It 
%   > overloads the QMol_DFT class to implement the SP-DFT model
%   > implements algebraic operations for DFT and related computation
%   > defines specific components for grid-based DFT computations
%     * Discretization
%       disc   discretization domain (QMol_grid class)
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT})
function showInfo
    fprintf('  * QMol_DFT_spinPol:\n      > 1D spin-polarized DFT\n'); 
    QMol_DFT_spinPol.version;
end
end
methods (Access=?QMol_DFT)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the DFT model

    % Implementation of the DFT model
    fprintf('  * Density-functional theory (DFT)                    spin polarized (SP)\n');
    fprintf('    SP-DFT level of theory, using the Kohn-Sham formalism [Kohn 1695].\n');
    obj.version;

    ref                 =   {'Kohn 1965'};

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,Transient,GetAccess=public,SetAccess=?QMol_suite)
    x                       % Implicit definition of the domain
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    xspan                   % x
end
properties (Constant,Access=public)
    isSpinPol           =   true
    dim                 =   1
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % x ~~~~~~~~~~~~~~~~
    function val = get.x(obj),  if obj.isInit,  val     =   obj.disc.x;     % obj.x should be empty after the object initialization
                                else,           val     =   obj.x;          end,    end
    function set.xspan(obj,val),                obj.x   =   val;            end
    function val = get.xspan(obj),              val     =   obj.x;          end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reset is handled by the parent (abstract) DFT class
methods (Access=public)
function initialize(obj,varargin)
%initialize initializes implementation-specific features
    
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
    initialize@QMol_DFT(obj,varargin{:});

end
end
methods (Hidden,Access = ?QMol_DFT)
function disc = getDiscCopy(obj)
%getDiscCopy returns a streamlined copy of the discretization (object)
    
    % Create discretization object
    if isempty(obj.disc)
        % Get copy from x (implicit definition)
        disc            =   QMol_disc('x',obj.x);
    elseif obj.disc.isBasis
        % Basis-set model
        disc            =   QMol_disc_basis('x',obj.disc.x,'v',obj.disc.v);
    else
        % Regular grid discretization
        disc            =   QMol_disc('x',obj.disc.x);
    end

    % Set QM component
    disc.QM             =   obj;

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    % Parent-class components
    [~,PropNames]       =   QMol_DFT.propertyNames;

    ClassName           =   'QMol_DFT_spinPol';
    PropNames           =   [PropNames,{'x','xspan'}];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

