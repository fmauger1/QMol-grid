classdef QMol_TDDFT_sympSplitOp < QMol_TDDFT
%QMol_TDDFT_sympSplitOp symplectic split-operator TDDFT propagator (for
%   local, i.e., functionally dependent of the one-body density and its
%   derivatives, only, not of the Kohn-Sham orbitals)

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDDFT})
function showInfo
    fprintf('  * QMol_TDDFT_sympSplitOp:\n');
    fprintf('      > TDDFT propagator\n'); 
    fprintf('      > Symplectic split-operator interface\n'); 
    QMol_TDDFT_sympSplitOp.version;
end
end
methods (Access=protected)
    ref = showDoc(obj)
    mem = getMemoryProfileOrbitalDensity(obj,opt)
    mem = getMemoryProfilePropagator(obj,opt)
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (GetAccess=public,SetAccess=?QMol_suite)
    % Split motif
    splitMotif              % (default) VTV or TVT
    % External field
    diffDT              =   1e-5
end
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % External field
    EFG
end
properties (Transient,GetAccess=public,SetAccess=?QMol_suite)
    % External fields (to be repurposed in EF)
    externalFieldGauge      % EFG ('none'/'off', 'length', or 'velocity')
    potentialVector         %
    electricField           %
    electricFieldDerivative %
end
properties (Access=protected)
    % Time propagation
    isVTV                   % whether to use the VTV scheme
    t_                      % internal time variable
    dt_                     % current time step in propagator
    % External fields
    FG                      % Field gauge (0 = ignore any field, 1 = length, 2 = velocity)
    FA                      % actual potential vector
    FE                      % actual electric field
    FDE                     % actual electric field derivative
    uA                      % use EF.potentialVector
    uE                      % use EF.electricField
    uDE                     % use EF.electricFieldDerivative
    % Miscellaenous
    xi                      % energy coordinate
end
properties (Transient,GetAccess=public,SetAccess=?protected)
    % Time propagation
    c_                      % [T,V] coefficients in the split
    expT                    % exponentiated kinetic-term propagator
    expV                    % exponentiated potential-term propagator
    expVup
    expVdw
    % Miscellaneous
    nKSO                    % Number of orbitals
    dv                      % Symplex volume (dx in 1D, dx*dy in 2D, dx*dy*dz in 3D)
    pKSO                    % Projection of the Kohn-Sham orbitals
    X                       % For dipole signal
    Y
    Z
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % EFG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function val = get.externalFieldGauge(obj),     val     =   obj.EFG;    end
    function set.externalFieldGauge(obj,val),       obj.EFG =   val;        end
    % potentialVector ~~~~~~~~~~~~~~~~~~~~~~
    function val = get.potentialVector(obj)
        if ~isempty(obj.FA),                        val     =   obj.FA;     %#ok<ALIGN> 
        elseif ~isempty(obj.EF),                    val     =   obj.EF.A;
        else,                                       val     =   [];         end
    end
    function set.potentialVector(obj,val),          obj.FA  =   val;        end
    % electricField ~~~~~~~~~~~~~~~~~~~~~~~~
    function val = get.electricField(obj)
        if ~isempty(obj.FE),                        val     =   obj.FE;     %#ok<ALIGN> 
        elseif ~isempty(obj.EF),                    val     =   obj.EF.E;
        else,                                       val     =   [];         end
    end
    function set.electricField(obj,val),            obj.FE  =   val;        end
    % electricFieldDerivative ~~~~~~~~~~~~~~
    function val = get.electricFieldDerivative(obj)
        if ~isempty(obj.FDE),                       val     =   obj.FDE;    %#ok<ALIGN> 
        elseif ~isempty(obj.EF),                    val     =   obj.EF.DE;
        else,                                       val     =   [];         end
    end
    function set.electricFieldDerivative(obj,val),  obj.FDE =   val;        end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_TDDFT(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.diffDT),     obj.diffDT  =   1e-5;                           end
end
end
methods (Access=protected)
    initializeChildren(obj,isRst)
    setOutputExternalField(obj,sN,opt)
    setOutputOrbitalDensity(obj,opt)
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    % Parent-class components
    [~,PropNames]       =   QMol_TDDFT.propertyNames;

    ClassName           =   'QMol_TDDFT_sympSplitOp';
    PropNames           =   [PropNames,{'splitMotif','EFG','diffDT','externalFieldGauge', ...
                                'potentialVector','electricField','electricFieldDerivative'}];    
end
end
%% Propagate TDDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=protected)
    % Propagator components
    setExpT(obj,d)
    setExpV(obj,c)
    % Output results
    addOutputExternalField(obj,sN,~,~)
    saveOutputOrbitalDensity(obj,k,t)
end
methods (Access=protected)
function E = getFE(obj,t)
%getFE compute and return the electric field
    
    if obj.uE,      E   =   obj.EF.electricField(t);
    else,           E   =  -(obj.EF.potentialVector(t+.5*obj.diffDT) - obj.EF.potentialVector(t-.5*obj.diffDT)) / obj.diffDT;
    end
end
function DE = getFDE(obj,t)
%getDFE compute and return the derivative of the electric field
    
    if obj.uDE,     DE  =   obj.EF.electricFieldDerivative(t);
    elseif obj.uE,  DE  =   (obj.EF.electricField(t+.5*obj.diffDT)-obj.EF.electricField(t-.5*obj.diffDT)) / obj.diffDT; 
    else,           DE  =   (2*obj.EF.potentialVector(t)-obj.EF.potentialVector(t+obj.diffDT)-obj.EF.potentialVector(t-obj.diffDT))/obj.diffDT^2;
    end
end
end

methods (Access=protected)
function saveRestartChildren(obj,~,~)
%saveOutputChildren
    
    % Copy expT, expV, ...
    obj.restart.c_      =   obj.c_;
    obj.restart.expT    =   obj.expT;
    obj.restart.expV    =   obj.expV;
    obj.restart.expVup  =   obj.expVup;
    obj.restart.expVdw  =   obj.expVdw;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

