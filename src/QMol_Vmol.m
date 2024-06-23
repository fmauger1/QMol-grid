classdef (Abstract) QMol_Vmol < QMol_suite
%QMol_Vmol molecular potential, as a collection of centers. It defines
%   > (abstract) interface class
%   > common documentation interface
%   > common components for atom-based molecular potentials
%     * Atomic centers
%       atom    cell array of the list of atomic centers
%     * run-time variables
%       nbA     number of atomic centers in the molecular potential
%               AKA: nbAtom
    
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
    fprintf( '  * QMol_VMol:\n');
    fprintf(['      > Interface for molecular potentials\n' ...
             '      > Defined as sum of atom-like centers\n']); 
    QMol_Vmol.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Implementation-specific documentation
    ref                 =   obj.showDoc;

    % Any atomic centers?
    if obj.nbA == 0,    return;         end

    % Parse list of (atomic) centers
    fprintf('  * Atom-like center(s)\n')
    ID                  =   cell(1,obj.nbA);
    for k = 1:obj.nbA
        ID{k}           =   obj.atom{k}.showDocumentation('param');
    end

    % Parse unique potential types
    [~,ind]             =   unique(ID);
    for k = 1:numel(ind)
        ref             =   [ref, obj.atom{ind(k)}.showDocumentation('pot')]; %#ok<AGROW> 
    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (GetAccess=public,SetAccess=?QMol_suite)
    % Atomic centers
    atom
end
properties (Hidden,Transient,GetAccess=public,SetAccess=?QMol_DFT)
    % Number of atomic centers
    nbA
end
properties (Dependent,GetAccess=public)
    numberAtom                      % nbA
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % nbA ~~~~~~~~~~~~~~
    function val = get.numberAtom(obj),         val         =   obj.nbA;    end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.nbA         =   [];

    % Also reset all (atomic) center components 
    % (cross-dependencies might be broken);
    if iscell(obj.atom)
        for k = 1:numel(obj.atom), obj.atom{k}.reset; end
    elseif ~isempty(obj.atom)
        obj.atom.reset;
    end
    
    % Initialization status
    obj.isInit          =   false;
end
function initialize(obj,QM)
%initialize initializes all the atom (center) components. For compatibility
%   with user-defined molecular potential, it supports having an empty atom
%   list
    
    % Initialization
    if nargin < 2,          QM  =   [];             end
    obj.reset;

    % Discretization
    if isempty(obj.atom)
        % Empty list of atom (for user-defined potential)
        obj.nbA         =   0;
    elseif iscell(obj.atom)
        % Initialize each atomic center individually
        obj.nbA         =   numel(obj.atom);
        for k = 1:obj.nbA, obj.atom{k}.initialize(QM); end
    else
        % Single atomic center (put back in a cell for support)
        obj.atom        =   {obj.atom};
        obj.nbA         =   1;
        obj.atom{1}.initialize(QM);
    end

    % House keeping
    obj.isInit       =   true;
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_Vmol';
    PropNames           =  {'atom'};
end
end
%% Interface methods (to be overloaded) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Abstract,Access=?QMol_Vmol)
    showDoc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

