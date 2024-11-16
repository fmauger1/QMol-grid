classdef QMol_SE_eig_basis < QMol_suite
%QMol_SE_eig_basis eigen-state solver for Schrodinger-equation (SE)
%   level of theory discretized on a basis.
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release
%   01.21.001   07/01/2024  F. Mauger
%       Add (missing) reference and funding information to ground-state
%       run time documentation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.001','07/01/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_SE_eig_basis})
function showInfo
    fprintf('  * QMol_SE_eig_basis:\n'); 
    fprintf('      > SE eigen-state solver for basis-set models\n      > MATLAB eig function\n'); 
    QMol_SE_eig_basis.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % General options
    fprintf('  * Eigen-state solver for SE Hamiltonians             MATLAB eig function\n')
    fprintf('    using a direct diagonalization of the Hamiltonian matrix.\n');
    
    % Version
    obj.version;
    ref                 =   {};
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    disp                =   true
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    display                 % disp
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Linked objects
    SE
end
properties (Transient,Access=private)
    % Currently running an eigen-state calculation
    isRun               =   false
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % disp ~~~~~~~~~~~~~
    function set.display(obj,val),              obj.disp    =   val;        end
    function val = get.display(obj),            val         =   obj.disp;   end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object
    
    % Run-time variables
    obj.SE              =   [];
    
    % Initialization status
    obj.isInit          =   false;
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.disp),   obj.disp    =   true;       end
end
end
methods (Access=public)
function initialize(obj,SE)
%initialize initializes the eigen solver
    
    % Initialization
    obj.reset;

    % Set links
    obj.SE             =   SE;

    % House keeping
    obj.isInit       =   true;
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_SE_eig_basis';
    PropNames           =   {'disp','display'};
end
end
%% Eigenstate computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function computeGroundState(obj,SE)
%computeGroundState computes the ground state(s) for the input Schrodinger-
%   equation model
    
    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.initialize(SE);                                                     % resets the object first
    obj.isRun   =   true;
    
    if obj.disp,    QMol_doc.showHeader;        end
    
    % (Re)link everything ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.disp
        QMol_doc.showSection('Build Schrodinger-equation (SE) model');
        SE.initialize(2);
    else
        SE.initialize(0);
    end
    
    if obj.disp                                                             % doc, ref, and funding info
        ref     =   obj.showDocumentation;  fprintf('\n');
        QMol_doc.showBibliography(ref);
        QMol_doc.showFunding;
    end

    % Compute eigen states ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.computeEigenstates;

    % Finalization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.disp
        
        % Show wave function and Schrodinger-equation energies
        QMol_doc.showSection('Wave-function energies');
        SE.showEnergy('wfcn');      fprintf('\n');
        
        QMol_doc.showSection('Schrodinger-equation-component energies');
        SE.showEnergy('SE');        fprintf('\n');
        
        % Footer
        QMol_doc.showFooter;
    end
    
    % Clear member temporary variables
    obj.isRun           =   false;
end
end
methods (Access=private)
function computeEigenstates(obj)
%computeEigenstates computes the eigenstates, and corresponding energies,
%   associated with the member SE object. The result eigenstates are stored
%   in the member SE's wave function.

    % Diagonalize Hamiltonian matrix
    [VV,E]              =   eig(obj.SE.disc.mT + obj.SE.V.mV);

    % Select and normalize wave functions
    [~,ind]             =   sort(real(diag(E)));
    ind                 =   ind(1:obj.SE.N);

    obj.SE.wfcn.wfcn    =   obj.SE.disc.SE_normalizeWaveFunction(VV(:,ind));
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

