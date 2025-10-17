classdef QMol_TDCI_sympSplitOp < QMol_TDCI
%QMol_TDCI_sympSplitOp symplectic split-operator TDCI propagator

%   Version     Date        Author
%   01.23.000   06/04/2025  F. Mauger
%       Creation (from QMol_TDSE_sympSplitOp)
%   01.23.001   06/07/2025  F. Mauger
%       Fix saving 'all' wave functions and run-time documentation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.001','06/07/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDCI})
function showInfo
    fprintf('  * QMol_TDCI_sympSplitOp:\n');
    fprintf('      > TDCI propagator\n'); 
    fprintf('      > Symplectic split-operator interface\n'); 
    QMol_TDCI_sympSplitOp.version;
end
end
methods (Access=protected)
    ref = showDoc(obj)
    mem = getMemoryProfileWaveFunction(obj,opt)
    mem = getMemoryProfilePropagator(obj,opt)
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (GetAccess=public,SetAccess=?QMol_suite)
    % Split motif
    splitMotif              % (default) HDH or DHD
    % External field
    diffDT              =   1e-5
end
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % External field
    EFG
end
properties (Transient,GetAccess=public,SetAccess=?QMol_suite)
    % External fields (to be repurposed in EF)
    externalFieldGauge      % EFG ('none'/'off' or 'length')
    potentialVector         %
    electricField           %
    electricFieldDerivative %
end
properties (Access=protected)
    % Time propagation
    isHDH                   % whether to use the HDH scheme
    t_                      % internal time variable
    dt_                     % current time step in propagator
    % External fields
    FG                      % Field gauge (0 = ignore any field, 1 = length)
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
    H0                      % CI matrix, including the CAP
    D                       % damping mask matrix
    % Miscellaneous
    nWfcn                   % Number of wave functions
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
    function set.potentialVector(obj,val),          obj.FA  =   val;        
                                                    obj.EF  =   [];         end
    % electricField ~~~~~~~~~~~~~~~~~~~~~~~~
    function val = get.electricField(obj)
        if ~isempty(obj.FE),                        val     =   obj.FE;     %#ok<ALIGN> 
        elseif ~isempty(obj.EF),                    val     =   obj.EF.E;
        else,                                       val     =   [];         end
    end
    function set.electricField(obj,val),            obj.FE  =   val;        
                                                    obj.EF  =   [];         end
    % electricFieldDerivative ~~~~~~~~~~~~~~
    function val = get.electricFieldDerivative(obj)
        if ~isempty(obj.FDE),                       val     =   obj.FDE;    %#ok<ALIGN> 
        elseif ~isempty(obj.EF),                    val     =   obj.EF.DE;
        else,                                       val     =   [];         end
    end
    function set.electricFieldDerivative(obj,val),  obj.FDE =   val;        
                                                    obj.EF  =   [];         end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_TDCI(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.diffDT),     obj.diffDT  =   1e-5;                           end
end
end
methods (Access=protected)
    initializeChildren(obj,isRst)
    setOutputExternalField(obj,sN,opt)
    setOutputWaveFunction(obj,opt)
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    % Parent-class components
    [~,PropNames]       =   QMol_TDCI.propertyNames;

    ClassName           =   'QMol_TDCI_sympSplitOp';
    PropNames           =   [PropNames,{'splitMotif','EFG','diffDT','externalFieldGauge', ...
                                'potentialVector','electricField','electricFieldDerivative'}];    
end
end
%% Propagate TDDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=protected)
    % Output results
    addOutputExternalField(obj,sN,~,~)
    saveOutputWaveFunction(obj,k,t)
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
    
    % Copy H0, D, etc.
    obj.restart.H0      =   obj.H0;
    obj.restart.D       =   obj.D;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

