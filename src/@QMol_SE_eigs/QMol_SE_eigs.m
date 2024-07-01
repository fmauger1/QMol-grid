classdef QMol_SE_eigs < QMol_suite
%QMol_SE_eigs eigen-state solver for Schrodinger-equation (SE) level of 
%   theory.
%
%   NOTES:
%     > QMol_SE_eigs uses a private random number stream to ensure 
%       reproducibility across runs
%     > QMol_SE_eigs imposes eigen state to be real valued
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release
%   01.21.001   07/01/2024  F. Mauger
%       Add (missing) reference and funcing information to ground-state
%       run time documentation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.001','07/01/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_ES_eigs})
function showInfo
    fprintf('  * QMol_SE_eigs:\n');
    fprintf('      > SE eigen-state solver\n      > MATLAB eigs function\n'); 
    QMol_SE_eigs.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % General options
    fprintf('  * Eigen-state solver for SE Hamiltonians            MATLAB eigs function\n')
    if ~obj.issym
        fprintf('    WARNING: The Hamiltonian operator is not assumed to be symmetric')
    end
    fprintf('    Tolerance  = %-6.2g\n', obj.tol);
    fprintf('    Max. iter. = %i\n', obj.maxit);
    fprintf('    Basis dim. = %i\n', obj.p);

    % Specific symmetry options
    obj.showSymmetry;
    
    % Version
    obj.version;
    ref                 =   {};
end
end
methods (Access=private)
    showSymmetry(obj)
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % EIGS options
    issym               =   true
    tol                 =   1e-12
    maxit               =   300
    p                   =   100
    disp                =   'final'
    sym
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    IsFunctionSymmetric     % issym
    Tolerance               % tol
    MaxIterations           % maxit
    SubspaceDimension       % p
    Display                 % disp
    Symmetry                % sym
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Linked objects
    SE
    % Parsed symmetry option
    symState
    % Level of display
    dispIter
    dispGS
    % Currently running a GS calculation
    isRun               =   false
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % issym ~~~~~~~~~~~~
    function set.IsFunctionSymmetric(obj,val),      obj.issym   =   val;        end
    function val = get.IsFunctionSymmetric(obj),    val         =   obj.issym;  end
    % tol ~~~~~~~~~~~~~~
    function set.Tolerance(obj,val),                obj.tol     =   val;        end
    function val = get.Tolerance(obj),              val         =   obj.tol;    end
    % maxit ~~~~~~~~~~~~
    function set.MaxIterations(obj,val),            obj.maxit   =   val;        end
    function val = get.MaxIterations(obj),          val         =   obj.maxit;  end
    % Vext ~~~~~~~~~~~~~
    function set.SubspaceDimension(obj,val),        obj.p       =   val;        end
    function val = get.SubspaceDimension(obj),      val         =   obj.p;      end
    % disp ~~~~~~~~~~~~~
    function set.Display(obj,val),                  obj.disp    =   val;        end
    function val = get.Display(obj),                val         =   obj.disp;   end
    % sym ~~~~~~~~~~~~~~
    function set.Symmetry(obj,val),                 obj.sym     =   val;        end
    function val = get.Symmetry(obj),               val         =   obj.sym;    end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object
    
    % Run-time variables
    obj.SE              =   [];
    obj.symState        =   [];
    obj.dispIter        =   [];
    obj.dispGS          =   [];
    obj.isRun           =   false;
    
    % Initialization status
    obj.isInit          =   false;
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.issym),  obj.issym   =   true;       end
    if isempty(obj.tol),    obj.tol     =   1e-12;      end
    if isempty(obj.maxit),  obj.maxit   =   300;        end
    if isempty(obj.p),      obj.p       =   100;        end
    if isempty(obj.disp),   obj.disp    =   'final';    end
end
end
methods (Access=public)
function initialize(obj,SE)
%initialize initializes all the atom (center) components. For compatibility
%   with user-defined molecular potential, it supports having an empty atom
%   list
    
    % Initialization
    obj.reset;

    % Set links
    obj.SE              =   SE;

    % Parse symmetry configuration
    obj.initializeSymmetry;

    % Level of display
    switch lower(obj.disp)
        case {'all','iter','iteration'}
            obj.dispIter=   true;
            obj.dispGS  =   true;
        case {'final','summary','on'}
            obj.dispIter=   false;
            obj.dispGS  =   true;
        case {'none','off'}
            obj.dispIter=   false;
            obj.dispGS  =   false;
        otherwise   % Unknown case (use default summary mode)
            obj.dispIter=   false;
            obj.dispGS  =   true;
    end

    % House keeping
    obj.isInit          =   true;
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_SE_eigs';
    PropNames           =  {'issym','tol','maxit','p','disp','sym', ...
                            'IsFunctionSymmetric','Tolerance','MaxIterations', ...
                            'SubspaceDimension','Display','Symmetry'};
end
end
methods (Access=private)
    % Symmetry configuration
    initializeSymmetry(obj)
    S = parseSymmetry_1D(~,symList)
end
%% Eigenstate computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
    computeGroundState(obj,SE)
end
methods (Access=private)
    computeEigenstates(obj)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

