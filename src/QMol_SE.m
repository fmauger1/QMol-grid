classdef QMol_SE < QMol_SEq
%QMol_SE implementation of the Schrodinger-equation (SE) discretization on 
%   a Cartesian grid. It 
%   > overloads the QMol_SEq class 
%   > implements algebraic operations for SE and related computation
%   > defines specific components for grid-based SE computations
%     * Discretization
%       disc   discretization domain (QMol_grid)
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_SEq})
function showInfo
    fprintf('  * QMol_SE:\n      > 1D SE\n'); 
    QMol_SE.version;
end
end
methods (Access=?QMol_SEq)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the DFT model

    % Implementation of the DFT model
    fprintf('  * Schrodinger equation (SE)\n');
    obj.version;

    ref                 =   {};

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
    initialize@QMol_SEq(obj,varargin{:});

end
end
methods (Hidden,Access = {?QMol_SE,?QMol_SEq})
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
    [~,PropNames]       =   QMol_SEq.propertyNames;

    ClassName           =   'QMol_SE';
    PropNames           =   [PropNames,{'x','xspan'}];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

