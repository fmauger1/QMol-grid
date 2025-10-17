classdef QMol_TDDFT_extSymp < QMol_TDDFT
%QMol_TDDFT_extSymp ext symplectic TDDFT propagator 

%   Version     Date        Author
%   01.23.000   07/21/2025  F. Mauger
%       Creation (from QMol_TDDFT_sympSplitOp)
%               07/22/2025  F. Mauger
%       Continue implementation
%               07/23/2025  F. Mauger
%       Fix restart
%       Add HRH, TVR, and TRV split motifs
%   01.23.001   07/30/2025  F. Mauger
%       Add HHR hybrid and HRH hybrid split motifs
%   01.23.002   08/06/20205 F. Mauger
%       Add option to save distance between extended phase-space copies 

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.002','08/06/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDDFT})
function showInfo
    fprintf('  * QMol_TDDFT_extSymp:\n');
    fprintf('      > TDDFT propagator\n'); 
    fprintf('      > Extended symplectic interface\n'); 
    QMol_TDDFT_extSymp.version;
end
end
methods (Access=protected)
    ref = showDoc(obj)
    mem = getMemoryProfileOrbitalDensity(obj,opt)
    mem = getMemoryProfilePropagator(obj,opt)
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (GetAccess=public,SetAccess=?QMol_suite)
    splitMotif          =   'TRV'    % split motif [ 'HHR' | 'HRH' | 'TVR' | 'TRV' (default) | 'HHR hybrid' | 'HRH hybrid' ]
    % External field
    diffDT              =   1e-5
end
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % External field
    EFG
    % Save distance between extended-phase space copies
    sDist               =   false
    sDistT
    oDist
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    saveDistanceCopy        % (sDist) save the distance between the extended phase-space copies [ true | false (default) ]
    saveDistanceCopyTime    % (sDistT) time sampling for saving the distance 
    outDistanceCopy         % (oDist) output structure for distance between the extended phase-space copies [ structure | [] (default) ]
end
properties (Transient,GetAccess=public,SetAccess=?QMol_suite)
    omega                   % (w) restrain parameter [ scalar (default 10) ]
    splitPotential          % (spltV) whether to split the potential operator [ true (default) | false ]
    externalFieldGauge      % EFG ('none'/'off', 'length', or 'velocity')
    potentialVector         %
    electricField           %
    electricFieldDerivative %
end
properties (Access=protected)
    % Time propagation
    algo                    % split motif to be used (1 = 'HHR', 2 = 'HRH', 3 = 'TVR', 4 = 'TRV', 5 = 'HHR hybrid')
    w                   =   10      % restrain
    spltV               =   true    % whether to split the potential
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
    expT                    % exponentiated kinetic-term propagator
    expV                    % exponentiated potential-term propagator
    expVup
    expVdw
    % Extended phase space
    pDFT                    % KSO attached to the DFT object
    p1                      % 1st extended phase-space KSOs
    p2                      % 2nd extended phase-space KSOs
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
    % w ~~~~~~~~~~~~~~~~
    function set.omega(obj,val),                    obj.w       =   val;    end
    function val = get.omega(obj),                  val         =   obj.w;  end
    % spltV ~~~~~~~~~~~~
    function set.splitPotential(obj,val),           obj.spltV   =   val;    end
    function val = get.splitPotential(obj),         val         =   obj.spltV; end

    % sDist ~~~~~~~~~~~~                          ==== Save distance =============
    function set.saveDistanceCopy(obj,val),         obj.sDist   =   val;        end
    function val = get.saveDistanceCopy(obj),       val         =   obj.sDist;  end
    % sDistT ~~~~~~~~~~~
    function set.saveDistanceCopyTime(obj,val),     obj.sDistT  =   val;        end
    function val = get.saveDistanceCopyTime(obj),   val         =   obj.sDistT; end
    % oDist ~~~~~~~~~~~~
    function set.outDistanceCopy(obj,val),          obj.oDist   =   val;        end
    function val = get.outDistanceCopy(obj),        val         =   obj.oDist;  end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_TDDFT(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.diffDT),     obj.diffDT  =   1e-5;                       end
    if isempty(obj.w),          obj.w       =   10;                         end
    if isempty(obj.splitMotif), obj.splitMotif= 'TRV';                      end
    if isempty(obj.sDist),      obj.sDist   =   false;                      end
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

    ClassName           =   'QMol_TDDFT_extSymp';
    PropNames           =   [PropNames,{'splitMotif','EFG','diffDT','externalFieldGauge', ...
                                'potentialVector','electricField','electricFieldDerivative', ...
                                'omega','w','splitPotential','spltV', ...
                                'sDist','sDistT','saveDistanceCopy','saveDistanceCopyTime', ...
                                'oDist','outDistanceCopy'}];    
end
end
%% Propagate TDDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=protected)
    % Propagator components
    setExpT(obj,d)
    setExpV(obj,p,c)
    % Output results
    addOutputExternalField(obj,sN,~,~)
    saveOutputOrbitalDensity(obj,k,t)
function ref = docOutputChildren(obj)           
%docOutputChildren

    % Output results
    if obj.sDist,   fprintf('  * Distance between extended phase-space copies\n');  end
    
    % Display warning for dipole velocity and acceleration
    if obj.sVel   ||   obj.sAcc
        fprintf('\n   WARNING: Calculations of the dipole velocity and acceleration only use\n');
        fprintf('     !!!!   the explicit part of the exchange-correlation functional and\n');
        fprintf('     !!!!   ignore contributions for the implicit one. This will lead to\n');
        fprintf('     !!!!   incorrect results for\n');
        fprintf('     !!!!   > orbital-resolved contributions with exact exchange and\n');
        fprintf('     !!!!   > meta GGA functional.\n');
    end
    ref = {};                   
end
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
    obj.restart.expT    =   obj.expT;
    obj.restart.expV    =   obj.expV;
    obj.restart.expVup  =   obj.expVup;
    obj.restart.expVdw  =   obj.expVdw;

    obj.DFT.KSO         =   obj.pDFT;                                       % Make sure DFT.KSO is properly aligned with pDFT
    obj.restart.pDFT    =   obj.pDFT;
    obj.restart.p1      =   obj.p1;
    obj.restart.p2      =   obj.p2;

end
function setOutputOrbital(obj)
%setOutputOrbital
    
    % Average orbitals
    if obj.DFT.isSpinPol,   obj.pDFT.KSOup  =   .5*(obj.p1.KSOup+obj.p2.KSOup); %#ok<ALIGN>
                            obj.pDFT.KSOdw  =   .5*(obj.p1.KSOdw+obj.p2.KSOdw);
    else,                   obj.pDFT.KSO    =   .5*(obj.p1.KSO+obj.p2.KSO);     end

    obj.DFT.KSO         =   obj.pDFT;
end
function setOutputChildren(obj,opt)
%setOutputChildren set oDist output parameters

    switch lower(opt)
        case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%
        
            % Clean up any old data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.oDist           =   [];
            
            % DFT energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if obj.sDist
                % When to save the distance
                obj.oDist.ind   =   [obj.getOutputIndex(obj.sDistT),NaN];
                obj.oDist.time  =   obj.tspan(obj.oDist.ind(1:end-1));
                obj.oDist.n     =   1;
        
                if obj.sEF,         obj.setOutputExternalField('oDist','init');     end
        
                % Distance data allocation
                obj.oDist.distance = NaN(1,numel(obj.oDist.time));
            else
                obj.oDist.ind   =   NaN;
                obj.oDist.n     =   1;
            end
        
        case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Main components
            if obj.sDist,   obj.oDist   =   rmfield(obj.oDist,{'ind','n'});         if obj.sEF %#ok<ALIGN> 
                            obj.setOutputExternalField('oDist','clean');            end
            else,           obj.oDist   =   [];                                     end
        
        otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Unexpected option
            error('QMol:TDDFT:setOutputChildren',['Unknown option ' opt]);
    end
end
end
%% Overwrite save and finalize methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%saveOutputResults already limits call of the various methods to when
%saving is required
methods (Access=protected)
function saveOutputDFT(obj,t),          obj.setOutputOrbital;           saveOutputDFT@QMol_TDDFT(obj,t);        end
function saveOutputEnergy(obj,k,t),     obj.setOutputOrbital;           saveOutputEnergy@QMol_TDDFT(obj,k,t);   end
function saveOutputFunction(obj,k,t),   obj.setOutputOrbital;           saveOutputFunction@QMol_TDDFT(obj,k,t); end
% At the end of the propagation, only return the average output orbitals
function finalize(obj),                 obj.setOutputOrbital;           finalize@QMol_TDDFT(obj);               end
function saveOutputChildren(obj,k,t)
%saveOutputChildren save oDist

    % Extended phase-space distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oDist.ind(obj.oDist.n)
        % Distance
        if obj.DFT.isSpinPol
            obj.oDist.distance(obj.oDist.n)     =   sqrt(sum(abs(obj.p1.KSOup-obj.p2.KSOup).^2,'all')*obj.dv + ...
                                                         sum(abs(obj.p1.KSOdw-obj.p2.KSOdw).^2,'all')*obj.dv);
        else
            obj.oDist.distance(obj.oDist.n)     =   sqrt(sum(abs(obj.p1.KSO  -obj.p2.KSO  ).^2,'all')*obj.dv);

        end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oDist',k,t);        end

        % Update counter
        obj.oDist.n     =   obj.oDist.n + 1;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

