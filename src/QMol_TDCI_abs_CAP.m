classdef QMol_TDCI_abs_CAP < QMol_suite
%QMol_TDCI_abs_CAP CAP-type absorbing-boundry damping matrix for TDCI
%   Use QMol_TDCI_abs_CAP to define the damping matrix modeling complex-
%   absorbing potential (CAP) type absorbing-boundary conditions in TDCI
%   calculations. The damping matrix shoud be specified in the same
%   configuration-state basis as for the (TD)CI model to be propagated.
%   Physical damping matrices are symmetric real, but the class does not
%   check or enforce this.
%
%   The CAP damping matrix is added to the CI matrix and thus automatically
%   scales with the propagation time step. The damping rate should be a
%   positive number, with 0 corresponding to no damping being applied.
%
%   Editable properties:
%     * dampingMatrix
%
%   Methods:
%     * Changing class properties: set, reset, clear
%     * Run-time documentation: showDocumentation, getMemoryProfile
%     * Damping matrix: getDampingMatrix
%
%   ABS = QMol_TDCI_abs_CAP('name1','value1',___) creates a damping matrix
%     with the name properties set to the specified values.
%
%   See also QMol_TDCI

%   Version     Date        Author
%   01.23.000   06/04/2025  F. Mauger
%       Creation
%   01.23.001   06/05/2025  F. Mauger
%       Add getDampingMatrix (for uniform interface with CAP that need to
%       build the damping matrix)
%   01.23.002   06/07/2025  F. Mauger
%       Fix isCAP flag

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.002','06/07/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDCI_abs_CAP})
function showInfo
    fprintf('  * QMol_TDCI_abs_CAP:\n      > Damping matrix for CAP-type absorbing boundary condition\n'); 
    QMol_TDCI_abs_CAP.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Absorber properties
    fprintf('  * Damping matrix (absorbing boundaries)      complex-absorbing potential\n');
    
    ref                 =   {};

    % Version
    obj.version;
 
end
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the absorber with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end
 
    % Evaluate (and display) estimate of memory footprint
    if isreal(obj.M),   mem     =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.M),'real');
    else,               mem     =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.M),'complex');   end

    if opt, QMol_DFT_profiler.showMemoryFootprint('Damping matrix (absorbing boundaries)', mem,1);  end
 
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    M                       % Damping matrix
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    dampingMatrix           % (M) damping matrix [ matrix (default []) ]
end
properties (Constant)
    isCAP               =   true
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % M ~~~~~~~~~~~~~~~~
    function set.dampingMatrix(obj,val),    obj.M       =   val;            end
    function val = get.dampingMatrix(obj),  val         =   obj.M;          end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function initialize(obj,~)
%initialize initializes the object

    % Update initialization status
    obj.isInit          =   true;

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_TDCI_abs_CAP';
    PropNames           =  {'M','dampingMatrix'};
end
end
%% Damping matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function M = getDampingMatrix(obj)
%getDampingMatrix get the damping matrix
%
%   M = obj.getDampingMatrix() returns the damping matrix.

    M                   =   obj.M;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

