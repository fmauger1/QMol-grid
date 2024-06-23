classdef (Abstract) QMol_TDSE < QMol_suite
%QMol_TDSE time propagation interface for Schrodinger-equation models. It 
%   defines:
%   > (abstract) interface class
%   > overall time propagation workflow
%     (subclasses need only implement elementary-step for propagation and
%     one-the-fly analysis/computations)
%   > common components for TDSE simulations

% QUESTIONS:
%   > do I want to add a a unique identifyer for each specific propagator

% NOTES:
%   > general template for definitions of what to save is
%     - saveWhat        =   activates saving of 'what' (e.g., wave function);
%                           generally a boolean (true to save)
%     - saveWhatTime    =   specifies the times at which 'what' should be
%                           saved. This should support:
%           vector      =   user-defined times at which to save
%           []          =   use the global time vector
%           'all'       =   save at every time step
%     - saveWhatIndex   =   if 'what' is a wave-function-resolved quantity,
%                           specify the indexes of the wave functions to 
%                           include in the saved quantity. This should 
%                           support:
%           vector      =   user-defined list of indexes
%           []          =   save result from all wave functions
%           'all'       =   save results from all wave functions
%     Then, the associated saved result is contained in the 'outWhat'
%     structure, which includes both the result and the time-sampling
%     vector (and the indexes of selected wave functionss for wave-
%     function-resolved outputs).

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDSE})
function showInfo
    fprintf('  * QMol_TDSE:\n      > TDSE interface\n'); 
    QMol_TDSE.version;
end
end
methods (Access=public)
    mem = getMemoryProfile(obj,opt)
    showDocumentation(obj)
end
methods (Abstract,Access=protected)
    ref = showDoc(obj)
end
methods (Access=protected)
function mem = getMemoryProfilePropagator(~,opt)

    if nargin < 2,  opt     =   false;  end
    if opt, QMol_SE_profiler.showMemoryFootprint('TDSE propagator has an undefined memory footprint', 0,1); end
    mem                 =   0; 

end
function mem = getMemoryProfileWaveFunction(obj,opt)

    if nargin < 2,  opt     =   false;  end
    
    if (obj.sWfcn || obj.sWfcnP)   &&   opt
        QMol_SE_profiler.showMemoryFootprint('Wave functions', 0,1);
        if obj.sWfcn,   QMol_DSE_profiler.showMemoryFootprint('Wave functions ouput has an undefined memory footprint', 0,2); end
        if obj.sWfcnP,  QMol_SE_profiler.showMemoryFootprint('Projection of the wave functions ouput has an undefined memory footprint', 0,2); end
    end

    mem                 =   0;

end
function mem = getMemoryProfileChildren(~,~),   mem = 0;                    end
function ref = docOutputChildren(~),            ref = {};                   end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Time propagation
    disp                =   true
    T
    dt                  =   1e-2
    % Gobbler
    ABC
    % External field
    EF
    % Save SE object into individual files
    sSE                 =   false
    sSEF                =   'QMolGrid--TDSE--SE'
    sSET
    % Save dipole/velocity/acceleration
    sDip                =   false
    sDipI
    sDipT
    sVel                =   false
    sVelI
    sVelT               =   'dipole'
    sAcc                =   false
    sAccI
    sAccT               =   'dipole'
    % Save energies
    sESE                =   false
    sESET
    sEWfcn              =   false
    sEWfcnT
    % Save external field
    sEF                 =   false
    % Save ionization
    sIon                =   false
    sIWfcnI
    sIonT
    % Save wave functions
    sWfcn               =   false
    sWfcnI
    sWfcnT
    sWfcnP              =   false
    sWfcnB
    sWfcnPI
    sWfcnPT
    % Save output functions
    sF
    sFT
    % Save restart data
    sRest               =   false
    sRestF              =   'QMolGrid--TDSE--Restart.mat'
    sRestT
    % Output results
    oSE
    oDip
    oVel
    oAcc
    oESE
    oEWfcn
    oIon
    oWfcn
    oWfcnP
    oF
    oRest
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    % Time propagation ==================
    time                                % T
    timeStep                            % dt
    display                             % disp
    % Gobbler ===========================
    absorbingBoundary                   % ABC
    % External field ====================
    externalField                       % EF
    % Save SE object into individual files
    saveSE                              % sSE
    saveSEFileName                      % sSEF
    saveSETime                          % sSET
    % Save dipole/velocity/acceleration =
    saveDipole                          % sDip
    saveDipoleWaveFunctionIndex         % sDipI
    saveDipoleTime                      % sDipT
    saveDipoleVelocity                  % sVel
    saveDipoleVelocityWaveFunctionIndex % sVelI
    saveDipoleVelocityTime              % sVelT
    saveDipoleAcceleration              % sAcc
    saveDipoleAccelerationWaveFunctionIndex  % sAccI
    saveDipoleAccelerationTime          % sAccT
    % Save energies =====================
    saveEnergySE                        % sESE
    saveEnergySETime                    % sESET
    saveEnergyWaveFunction              % sEWfcn
    saveEnergyWaveFunctionTime          % sEWfcnT
    % Save external field ===============
    saveExternalField                   % sEF
    % Save ionization ===================
    saveIonization                      % sIon
    saveIonizationWaveFunctionIndex     % sIWfcnI
    saveIonizationTime                  % sIonT
    % Save wafe function ================
    saveWaveFunction                    % sWfcn
    saveWaveFunctionIndex               % sWfcnI
    saveWaveFunctionTime                % sWfcnT
    saveWaveFunctionProjection          % sWfcnP
    saveWaveFunctionProjectionBasis     % sWfcnB
    saveWaveFunctionProjectionIndex     % sWfcnPI
    saveWaveFunctionProjectionTime      % sWfcnPT
    % Save output functions =============
    saveOutputFunction                  % sF
    saveOutputFunctionTime              % sFT
    % Save restart data =================
    saveRestart                         % sRest
    saveRestartFileName                 % sRestF
    saveRestartTime                     % sRestT
    % Output results ====================
    outSE                               % oSE
    outDipole                           % oDip
    outDipoleVelocity                   % oVel
    outDipoleAcceleration               % oAcc
    outEnergySE                         % oESE
    outEnergyWaveFunction               % oEWfcn
    outIonization                       % oIon
    outWaveFunction                     % oWfcn
    outWaveFunctionProjection           % oWfcnP
    outOutputFunction                   % oF
    outRestart                          % oRest
end
properties (Access=protected)
    % Time propagation
    SE
    tspan
    iref
    % Restart
    restart
end
properties (Transient,Access=private)
    isRun               =   false
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % T ~~~~~~~~~~~~~~~~                          ==== Time =======================
    function set.time(obj,val),                     obj.T       =   val;        end
    function val = get.time(obj),                   val         =   obj.T;      end
    % dt ~~~~~~~~~~~~~~~
    function set.timeStep(obj,val),                 obj.dt      =   val;        end
    function val = get.timeStep(obj),               val         =   obj.dt;     end
    % disp ~~~~~~~~~~~~~
    function set.display(obj,val),                  obj.disp    =   val;        end
    function val = get.display(obj),                val         =   obj.disp;   end
    % ABC ~~~~~~~~~~~~~~                          ==== Gobbler ====================
    function set.absorbingBoundary(obj,val),        obj.ABC     =   val;        end
    function val = get.absorbingBoundary(obj),      val         =   obj.ABC;    end
    % EF ~~~~~~~~~~~~~~~                          ==== External field =============
    function set.externalField(obj,val),            obj.EF      =   val;        end
    function val = get.externalField(obj),          val         =   obj.EF;     end
    % sSE ~~~~~~~~~~~~~~                          ==== Save SE to files ===========
    function set.saveSE(obj,val),                   obj.sSE     =   val;        end
    function val = get.saveSE(obj),                 val         =   obj.sSE;    end
    % sSEF ~~~~~~~~~~~~~
    function set.saveSEFileName(obj,val),           obj.sSEF    =   val;        end
    function val = get.saveSEFileName(obj),         val         =   obj.sSEF;   end
    % sSET ~~~~~~~~~~~~~
    function set.saveSETime(obj,val),               obj.sSET    =   val;        end
    function val = get.saveSETime(obj),             val         =   obj.sSET;   end
    % sDip ~~~~~~~~~~~~~                          ==== Save dipole/vel/acc ========
    function set.saveDipole(obj,val),               obj.sDip    =   val;        end
    function val = get.saveDipole(obj),             val         =   obj.sDip;   end
    % sDipI ~~~~~~~~~~~~
    function set.saveDipoleWaveFunctionIndex(obj,val),   obj.sDipI   =   val;        end
    function val = get.saveDipoleWaveFunctionIndex(obj), val         =   obj.sDipI;  end
    % sDipT ~~~~~~~~~~~~
    function set.saveDipoleTime(obj,val),           obj.sDipT   =   val;        end
    function val = get.saveDipoleTime(obj),         val         =   obj.sDipT;  end
    % sVel ~~~~~~~~~~~~~
    function set.saveDipoleVelocity(obj,val),       obj.sVel    =   val;        end
    function val = get.saveDipoleVelocity(obj),     val         =   obj.sVel;   end
    % sVelI ~~~~~~~~~~~~
    function set.saveDipoleVelocityWaveFunctionIndex(obj,val), obj.sVelI=val;        end
    function val = get.saveDipoleVelocityWaveFunctionIndex(obj), val =   obj.sVelI;  end
    % sVelT ~~~~~~~~~~~~
    function set.saveDipoleVelocityTime(obj,val),   obj.sVelT   =   val;        end
    function val = get.saveDipoleVelocityTime(obj), val         =   obj.sVelT;  end
    % sAcc ~~~~~~~~~~~~~
    function set.saveDipoleAcceleration(obj,val),   obj.sAcc    =   val;        end
    function val = get.saveDipoleAcceleration(obj), val         =   obj.sAcc;   end
    % sAccI ~~~~~~~~~~~~
    function set.saveDipoleAccelerationWaveFunctionIndex(obj,val), obj.sAccI=val;    end
    function val = get.saveDipoleAccelerationWaveFunctionIndex(obj), val=obj.sAccI;  end
    % sVelT ~~~~~~~~~~~~
    function set.saveDipoleAccelerationTime(obj,val), obj.sAccT =   val;        end
    function val = get.saveDipoleAccelerationTime(obj), val     =   obj.sAccT;  end
    % sESE ~~~~~~~~~~~~~                          ==== Save energies ==============
    function set.saveEnergySE(obj,val),             obj.sESE    =   val;        end
    function val = get.saveEnergySE(obj),           val         =   obj.sESE;   end
    % sESET ~~~~~~~~~~~~
    function set.saveEnergySETime(obj,val),         obj.sESET   =   val;        end
    function val = get.saveEnergySETime(obj),       val         =   obj.sESET;  end
    % sEWfcn ~~~~~~~~~~~
    function set.saveEnergyWaveFunction(obj,val),   obj.sEWfcn  =   val;        end
    function val = get.saveEnergyWaveFunction(obj), val         =   obj.sEWfcn; end
    % sEWfcnT ~~~~~~~~~~
    function set.saveEnergyWaveFunctionTime(obj,val),    obj.sEWfcnT  =   val;        end
    function val = get.saveEnergyWaveFunctionTime(obj),  val          =   obj.sEWfcnT; end
    % sEF ~~~~~~~~~~~~~~                          ==== Save external field ========
    function set.saveExternalField(obj,val),        obj.sEF     =   val;        end
    function val = get.saveExternalField(obj),      val         =   obj.sEF;    end
    % sIon ~~~~~~~~~~~~~                          ==== Save ionization ============
    function set.saveIonization(obj,val),           obj.sIon    =   val;        end
    function val = get.saveIonization(obj),         val         =   obj.sIon;   end
    % sIWfcnI ~~~~~~~~~~~
    function set.saveIonizationWaveFunctionIndex(obj,val), obj.sIWfcnI=   val;        end
    function val = get.saveIonizationWaveFunctionIndex(obj), val     =   obj.sIWfcnI; end
    % sIonT ~~~~~~~~~~~~
    function set.saveIonizationTime(obj,val),       obj.sIonT   =   val;        end
    function val = get.saveIonizationTime(obj),     val         =   obj.sIonT;  end
    % sWfcn ~~~~~~~~~~~~
    function set.saveWaveFunction(obj,val),         obj.sWfcn   =   val;        end
    function val = get.saveWaveFunction(obj),       val         =   obj.sWfcn;  end
    % sWfcnI ~~~~~~~~~~~
    function set.saveWaveFunctionIndex(obj,val),    obj.sWfcnI  =   val;        end
    function val = get.saveWaveFunctionIndex(obj),  val         =   obj.sWfcnI; end
    % sWfcnT ~~~~~~~~~~~
    function set.saveWaveFunctionTime(obj,val),     obj.sWfcnT   =   val;        end
    function val = get.saveWaveFunctionTime(obj),   val          =   obj.sWfcnT; end
    % sWfcnP ~~~~~~~~~~~
    function set.saveWaveFunctionProjection(obj,val),    obj.sWfcnP   =   val;        end
    function val = get.saveWaveFunctionProjection(obj),  val         =   obj.sWfcnP;  end
    % sWfcnB ~~~~~~~~~~~
    function set.saveWaveFunctionProjectionBasis(obj,val), obj.sWfcnB=   val;        end
    function val = get.saveWaveFunctionProjectionBasis(obj), val     =   obj.sWfcnB; end
    % sWfcnPI ~~~~~~~~~~
    function set.saveWaveFunctionProjectionIndex(obj,val), obj.sWfcnPI=   val;        end
    function val = get.saveWaveFunctionProjectionIndex(obj), val     =   obj.sWfcnPI; end
    % sWfcnPT ~~~~~~~~~~
    function set.saveWaveFunctionProjectionTime(obj,val), obj.sWfcnPT =   val;        end
    function val = get.saveWaveFunctionProjectionTime(obj), val      =   obj.sWfcnPT; end
    % sFWfcn ~~~~~~~~~~~
    function set.saveOutputFunction(obj,val),       obj.sF      =   val;        end
    function val = get.saveOutputFunction(obj),     val         =   obj.sF;     end
    % sFWfcnT ~~~~~~~~~~~
    function set.saveOutputFunctionTime(obj,val), obj.sFT=val;        end
    function val = get.saveOutputFunctionTime(obj), val  =   obj.sFT; end
    % sRest ~~~~~~~~~~~~                          ==== Save restart data ==========
    function set.saveRestart(obj,val),              obj.sRest   =   val;        end
    function val = get.saveRestart(obj),            val         =   obj.sRest;  end
    % sRestF ~~~~~~~~~~~
    function set.saveRestartFileName(obj,val),      obj.sRestF  =   val;        end
    function val = get.saveRestartFileName(obj),    val         =   obj.sRestF; end
    % sRestT ~~~~~~~~~~~
    function set.saveRestartTime(obj,val),          obj.sRestT  =   val;        end
    function val = get.saveRestartTime(obj),        val         =   obj.sRestT; end
    % oSE ~~~~~~~~~~~~~~                         ==== Output results ==============
    function set.outSE(obj,val),                    obj.oSE     =   val;        end
    function val = get.outSE(obj),                  val         =   obj.oSE;    end
    % oDip ~~~~~~~~~~~~~
    function set.outDipole(obj,val),                obj.oDip    =   val;        end
    function val = get.outDipole(obj),              val         =   obj.oDip;   end
    % oVel ~~~~~~~~~~~~~
    function set.outDipoleVelocity(obj,val),        obj.oVel    =   val;        end
    function val = get.outDipoleVelocity(obj),      val         =   obj.oVel;   end
    % oAcc ~~~~~~~~~~~~~
    function set.outDipoleAcceleration(obj,val),    obj.oAcc    =   val;        end
    function val = get.outDipoleAcceleration(obj),  val         =   obj.oAcc;   end
    % oESE ~~~~~~~~~~~~~
    function set.outEnergySE(obj,val),              obj.oESE    =   val;        end
    function val = get.outEnergySE(obj),            val         =   obj.oESE;   end
    % oEWfcn ~~~~~~~~~~~
    function set.outEnergyWaveFunction(obj,val),    obj.oEWfcn  =   val;        end
    function val = get.outEnergyWaveFunction(obj),  val         =   obj.oEWfcn; end
    % oIon ~~~~~~~~~~~~~
    function set.outIonization(obj,val),            obj.oIon    =   val;        end
    function val = get.outIonization(obj),          val         =   obj.oIon;   end
    % oWfcn ~~~~~~~~~~~~
    function set.outWaveFunction(obj,val),          obj.oWfcn   =   val;        end
    function val = get.outWaveFunction(obj),        val         =   obj.oWfcn;  end
    % oWfcnP ~~~~~~~~~~~
    function set.outWaveFunctionProjection(obj,val),     obj.oWfcnP   =   val;        end
    function val = get.outWaveFunctionProjection(obj),   val         =   obj.oWfcnP;  end
    % oF ~~~~~~~~~~~~~~~
    function set.outOutputFunction(obj,val), obj.oF   =   val;        end
    function val = get.outOutputFunction(obj),val        =   obj.oF;  end
    % oFRest ~~~~~~~~~~~
    function set.outRestart(obj,val),               obj.oRest   =   val;        end
    function val = get.outRestart(obj),             val         =   obj.oRest;  end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.disp),       obj.disp    =   true;                           end
    if isempty(obj.dt),         obj.dt      =   1e-2;                           end
    if isempty(obj.sSE),        obj.sSE     =   false;                          end
    if isempty(obj.sSEF),       obj.sSEF    =   'QMolGrid--TDSE--SE';           end
    if isempty(obj.sDip),       obj.sDip    =   false;                          end
    if isempty(obj.sVel),       obj.sVel    =   false;                          end
    if isempty(obj.sVelT),      obj.sVelT   =   'dipole';                       end
    if isempty(obj.sAcc),       obj.sAcc    =   false;                          end
    if isempty(obj.sAccT),      obj.sAccT   =   'dipole';                       end
    if isempty(obj.sESE),       obj.sESE    =   false;                          end
    if isempty(obj.sEWfcn),     obj.sEWfcn  =   false;                          end
    if isempty(obj.sEF),        obj.sEF     =   false;                          end
    if isempty(obj.sIon),       obj.sIon    =   false;                          end
    if isempty(obj.sWfcn),      obj.sWfcn   =   false;                          end
    if isempty(obj.sWfcnP),     obj.sWfcnP  =   false;                          end
    if isempty(obj.sRestF),     obj.sRestF  =   'QMolGrid--TDSE--Restart.mat';  end
    if isempty(obj.sRest),      obj.sRest   =   false;                          end
end
% Methods in separate files
    initialize(obj,SE,opt)
end
methods (Access=protected)
function ind = getOutputIndex(obj,t)
%getOutputIndex convert sampling time(s) to time indexes
    
    % Sempling indexes
    if isempty(t),              ind     =   obj.iref;                       % Default time sampling
    elseif ischar(t)
        if strcmpi(t,'all'),    ind     =   1:numel(obj.tspan);             % All time steps
        else,                   ind     =   obj.iref;
            warning('QMol:TDSE:outputIndex',['Unknown sampling option ' t'; default sampling used instead']);
        end
    elseif numel(t) > 1,        ind     =   obj.findOutputIndex(t);         % User defined time steps
    elseif t > 0,               ind     =   obj.T(1):sign(obj.dt)*t:obj.T(end)+.5*obj.dt;
                                ind     =   obj.findOutputIndex(ind);       % User defined time sampling
    else,                       ind     =   1:round(-t):numel(obj.tspan);   % User defined time step sampling
    end
end
function initializeChildren(~,~);       end                                 % Initialize children
function setOutputChildren(~,~);        end                                 % Children class implement additional output
function setOutputExternalField(~,~,~); end                                 % Enable saving of output external field
% Methods in separate files
    finalize(obj)
    setOutputSE(obj,opt)
    setOutputDipole(obj,opt)
    setOutputEnergy(obj,opt)
    setOutputFunction(obj,opt)
    setOutputIonization(obj,opt)
    setOutputRestart(obj,opt)
end
methods (Abstract,Access=protected)
    setOutputWaveFunction(obj,opt)
end
methods (Access=private)
function ind = findOutputIndex(obj,t)
%findOutputIndex finds the propagation indexes matching the input times
    
    ind                 =   uniquetol(interp1(obj.tspan,1:numel(obj.tspan),t,'nearest','extrap'),1e-10,'DataScale',1);
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_TDSE';

    PropNames           =  {... Time ======================================
                            'disp','display',                       'T','time',                                         'dt','timeStep', ...
                            ... Gobbler ===================================
                            'ABC','absorbingBoundary',                      ...
                            ... External field ============================
                            'EF','externalField',                           ...
                            ... Save SE object into individual files =====
                            'sSE','saveSE',                             'sSEF','saveSEFileName',                            'sSET','saveSETime', ...
                            ... Save dipole/velocity/acceleration =========
                            'sDip','saveDipole',                        'sDipI','saveDipoleWaveFunctionIndex',              'sDipT','saveDipoleTime', ...
                            'sVel','saveDipoleVelocity',                'sVelI','saveDipoleVelocityWaveFunctionIndex',      'sVelT','saveDipoleVelocityTime', ...
                            'sAcc','saveDipoleAcceleration',            'sAccI','saveDipoleAccelerationWaveFunctionIndex',  'sAccT','saveDipoleAccelerationTime', ...
                            ... Save energies =============================
                            'sESE','saveEnergySE',                      'sESET','saveEnergySETime', ...
                            'sEWfcn','saveEnergyWaveFunction',          'sEWfcnT','saveEnergyWaveFunctionTime', ...
                            ... Save external field =======================
                            'sEF','saveExternalField', ...
                            ... Save ionization ===========================
                            'sIon','saveIonization',                    'sIWfcnI','saveIonizationWaveFunctionIndex',        'sIonT','saveIonizationTime', ...
                            ... Save wave function ========================
                            'sWfcn','saveWaveFunction',                 'sWfcnT','saveWaveFunctionTime',                    'sWfcnI','saveWaveFunctionIndex', ...
                            'sWfcnP','saveWaveFunctionProjection',      'sWfcnB','saveWaveFunctionProjectionBasis', ...
                            'sWfcnPI','saveWaveFunctionProjectionIndex','sWfcnPT','saveWaveFunctionProjectionTime',...
                            ... Save output functions =====================
                            'sF','saveOutputFunction',                  'sFT','saveOutputFunctionTime', ...
                            ... Save restart data =========================
                            'sRest','saveRestart',                      'sRestF','saveRestartTime',                         'sRestT','saveRestartFileName', ...
                            ... Output results ============================
                            'oSE','outSE', ...
                            'oDip','outDipole',                         'oVel','outDipoleVelocity',                         'oAcc','outDipoleAcceleration', ...
                            'oESE','outEnergySE',                       'oEWfcn','outEnergyWaveFunction', ...
                            'oIon','outIonization', ...
                            'oWfcn','outWaveFunction',                  'oWfcnP','outWaveFunctionProjection', ...
                            'oF','outOutputFunction', ...
                            'oRest','outRestart'};
    
end
end
%% Propagate TDDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
    propagate(obj,opt)
end
methods (Access=protected)
    saveOutputSE(obj,t)
    saveOutputEnergy(obj,k,t)
    saveOutputFcn(obj,k,t)                  % Resolve conflicting name with output structure
function saveOutputChildren(~,~,~);         end                             % Children class implement additional output
function addOutputExternalField(~,~,~,~);   end                             % Add the external field info to an output structure
function saveRestartChildren(~,~,~);        end                             % Children class save needed info in restart structure
end
methods (Abstract,Access=protected)
    setTimeStep(obj,dt,t)
    applyTimeStep(obj,t)
    saveOutputDipole(obj,k,t)
    [E,DE] = getExternalFieldEnergy(obj,t)
    saveOutputIonization(obj,k,t)
    saveOutputWaveFunction(obj,k,t);
end
methods (Access=private)
    saveOutputResults(obj,k,t)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

