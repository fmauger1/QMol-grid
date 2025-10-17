classdef (Abstract) QMol_TDCI < QMol_suite
%QMol_TDCI interface for time-dependent configuration-interaction (TDCI)
%   This class provides the common backbone for the various types of TDCI
%   propagators available in the QMol-grid package.
%
%   Available TDCI propagators are
%     * QMol_TDCI_SSO       symplectic split-operator.
%
%   For large CI models, with a few hundreds or few thousands of
%   configuration states in the CI basis, using sparse matrices can
%   drastically speed up TDCI calculations. Use the sparsify method of the
%   CI object to be propagated to convert the CI and dipole-coupling
%   matrices to sparse ones. Alternatively, TDCI propagators offer the
%   option to perform the propagation in a subspace of the system,
%   specified as a basis of eigen states of the CI matrix. Note that the
%   sub-space propagation systematically lead to dense matrices, which can
%   slow down the propagation.
%
%   TDCI propagators support on-the-fly calculation and saving of a number 
%   of observables and quantities during the propagation. These are listed
%   below, in the "save" editable properties. The general name convention 
%   and template for supported options are (see each propagator 
%   documentation for specific options and support)
%     * save<What>      Activates the calculation and saving of <what> and
%                       is generally a boolean. If not activated, the
%                       associated time and index properties are irrelevant
%                       and ignored.
%     * save<What>Time  Specifies the times at which <what> should be
%                       saved. Supported types of entries generally are
%           vector  ->  user defined times at which to save
%           scalar  ->  a positive scalar specifies the sampling time step
%                       between successive saving times. For forward time 
%                       propagation, it is equivalent to 
%                           time(1):save<What>Time:time(end)
%                       and for backward propagation, it is equivalent to 
%                           time(1):-save<What>Time:time(end).
%           'all'   ->  save at every time step
%           []      ->  use the global time vector 
%     * save<What>Index Specifies the indexes of the wave functions for
%                       which to save the results. In general, supported
%                       entries are
%           vector  ->  user defined list of indexes
%           'all'   ->  saves the results for all wave functions
%           []      ->  depending on the type of output, either saves the
%                       results for all wave function, i.e., equivalent to
%                       'all', or deactivates saving altogether, i.e.,
%                       equivalent to save<What>=false. The latter case is
%                       a guardrail against output that may lead to very
%                       large memory requirements.
%   The output of each observable or quantity saved in then stored in an
%   output structure in the propagator object, with name out<What>. In this
%   structure,
%     * out<What>.time  contains the exact times at which the quantities
%                       have been saved.
%     * if the electric field was recorded (saveExternalField=true), 
%                       out<What>.electricField, out<What>.vectorPotential,
%                       and out<What>.electricFieldDerivative contain the
%                       corresponding values for the electric field, vector
%                       potential, and derivative of the electric field,
%                       respectively.
%                       Note: dependig on the choice of gauge and input
%                       driving field, only a subset of these may be saved.  
%
%   Editable properties (actual propagators may define additional ones):
%     * Time propagation: display, time, timeStep, absorbingBoundary, 
%           externalField
%     * Save the configuration-interaction object into individual files: 
%           saveCI, saveCIFileName, saveCITime
%     * Save the dipole, dipole velocity, and dipole acceleration signals:
%           saveDipole, saveDipoleWaveFunctionIndex, saveDipoleTime,
%           saveDipoleVelocity, saveDipoleVelocityWaveFunctionIndex,
%           saveDipoleVelocityTime, saveDipoleAcceleration,
%           saveDipoleAccelerationWaveFunctionIndex, 
%           saveDipoleAccelerationTime
%     * Save the configuration-interaction and wave function energies: 
%           saveEnergyCI, saveEnergyCITime, saveEnergyWaveFunction,
%           saveEnergyWaveFunctionTime
%     * Save the external field information: saveExternalField
%     * Save the ionization signal: saveIonization, 
%           saveIonizationWaveFunctionIndex, saveIonizationTime
%     * Save the output of a functions of the wave function: 
%           saveOutputFunction, saveOutputFunctionTime
%     * Save restart data file: saveRestart, saveRestartFileName, 
%           saveRestartTime
%
%   Output results (read-only properties): outCI, outDipole, 
%       outDipoleVelocity, outDipoleAcceleration, outEnergyCI, 
%       outEnergyWaveFunction, outIonization, outWaveFunction,
%       outWaveFunctionProjection, outOutputFunction, outRestart
%
%   Methods:
%     * Changing class properties: set, reset, clear
%     * Initialization: initialize
%     * Run-time documentation: showDocumentation, getMemoryProfile
%     * Time propagation: propagate
%   
%   See also QMol_TDCI_SSO

%   Version     Date        Author
%   01.23.000   05/25/2025  F. Mauger
%       Creation (from QMol_TDSE)

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.000','05/25/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDCI})
function showInfo
    fprintf('  * QMol_TDCI:\n      > TDCI interface\n'); 
    QMol_TDCI.version;
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
    if opt, QMol_CI_profiler.showMemoryFootprint('TDCI propagator has an undefined memory footprint', 0,1); end
    mem                 =   0; 

end
function mem = getMemoryProfileWaveFunction(obj,opt)

    if nargin < 2,  opt     =   false;  end
    
    if (obj.sWfcn || obj.sWfcnP)   &&   opt
        QMol_CI_profiler.showMemoryFootprint('Wave functions', 0,1);
        if obj.sWfcn,   QMol_CI_profiler.showMemoryFootprint('Wave functions ouput has an undefined memory footprint', 0,2); end
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
    pV
    % Gobbler
    ABC
    % External field
    EF
    % Save CI object into individual files
    sCI                 =   false
    sCIF                =   'QMolGrid--TDCI--CI'
    sCIT
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
    sECI                =   false
    sECIT
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
    % Save output functions
    sF
    sFT
    % Save restart data
    sRest               =   false
    sRestF              =   'QMolGrid--TDCI--Restart.mat'
    sRestT
    % Output results
    oCI
    oDip
    oVel
    oAcc
    oECI
    oEWfcn
    oIon
    oWfcn
    oF
    oRest
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    % Time propagation ==================
    time                                % T
    timeStep                            % dt
    propagationSpace                    % pV
    display                             % disp
    % Gobbler ===========================
    absorbingBoundary                   % ABC
    % External field ====================
    externalField                       % EF
    % Save CI object into individual files
    saveCI                              % sCI
    saveCIFileName                      % sCIF
    saveCITime                          % sCIT
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
    saveEnergyCI                        % sECI
    saveEnergyCITime                    % sECIT
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
    % Save output functions =============
    saveOutputFunction                  % sF
    saveOutputFunctionTime              % sFT
    % Save restart data =================
    saveRestart                         % sRest
    saveRestartFileName                 % sRestF
    saveRestartTime                     % sRestT
    % Output results ====================
    outCI                               % oCI
    outDipole                           % oDip
    outDipoleVelocity                   % oVel
    outDipoleAcceleration               % oAcc
    outEnergyCI                         % oECI
    outEnergyWaveFunction               % oEWfcn
    outIonization                       % oIon
    outWaveFunction                     % oWfcn
    outOutputFunction                   % oF
    outRestart                          % oRest
end
properties (Access=protected)
    % Time propagation
    CI
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
    % pV ~~~~~~~~~~~~~~~
    function set.propagationSpace(obj,val),         obj.pV      =   val;        end
    function val = get.propagationSpace(obj),       val         =   obj.pV;     end
    % disp ~~~~~~~~~~~~~
    function set.display(obj,val),                  obj.disp    =   val;        end
    function val = get.display(obj),                val         =   obj.disp;   end
    % ABC ~~~~~~~~~~~~~~                          ==== Gobbler ====================
    function set.absorbingBoundary(obj,val),        obj.ABC     =   val;        end
    function val = get.absorbingBoundary(obj),      val         =   obj.ABC;    end
    % EF ~~~~~~~~~~~~~~~                          ==== External field =============
    function set.externalField(obj,val),            obj.EF      =   val;        end
    function val = get.externalField(obj),          val         =   obj.EF;     end
    % sCI ~~~~~~~~~~~~~~                          ==== Save CI to files ===========
    function set.saveCI(obj,val),                   obj.sCI     =   val;        end
    function val = get.saveCI(obj),                 val         =   obj.sCI;    end
    % sCIF ~~~~~~~~~~~~~
    function set.saveCIFileName(obj,val),           obj.sCIF    =   val;        end
    function val = get.saveCIFileName(obj),         val         =   obj.sCIF;   end
    % sCIT ~~~~~~~~~~~~~
    function set.saveCITime(obj,val),               obj.sCIT    =   val;        end
    function val = get.saveCITime(obj),             val         =   obj.sCIT;   end
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
    % sECI ~~~~~~~~~~~~~                          ==== Save energies ==============
    function set.saveEnergyCI(obj,val),             obj.sECI    =   val;        end
    function val = get.saveEnergyCI(obj),           val         =   obj.sECI;   end
    % sECIT ~~~~~~~~~~~~
    function set.saveEnergyCITime(obj,val),         obj.sECIT   =   val;        end
    function val = get.saveEnergyCITime(obj),       val         =   obj.sECIT;  end
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
    % oCI ~~~~~~~~~~~~~~                         ==== Output results ==============
    function set.outCI(obj,val),                    obj.oCI     =   val;        end
    function val = get.outCI(obj),                  val         =   obj.oCI;    end
    % oDip ~~~~~~~~~~~~~
    function set.outDipole(obj,val),                obj.oDip    =   val;        end
    function val = get.outDipole(obj),              val         =   obj.oDip;   end
    % oVel ~~~~~~~~~~~~~
    function set.outDipoleVelocity(obj,val),        obj.oVel    =   val;        end
    function val = get.outDipoleVelocity(obj),      val         =   obj.oVel;   end
    % oAcc ~~~~~~~~~~~~~
    function set.outDipoleAcceleration(obj,val),    obj.oAcc    =   val;        end
    function val = get.outDipoleAcceleration(obj),  val         =   obj.oAcc;   end
    % oECI ~~~~~~~~~~~~~
    function set.outEnergyCI(obj,val),              obj.oECI    =   val;        end
    function val = get.outEnergyCI(obj),            val         =   obj.oECI;   end
    % oEWfcn ~~~~~~~~~~~
    function set.outEnergyWaveFunction(obj,val),    obj.oEWfcn  =   val;        end
    function val = get.outEnergyWaveFunction(obj),  val         =   obj.oEWfcn; end
    % oIon ~~~~~~~~~~~~~
    function set.outIonization(obj,val),            obj.oIon    =   val;        end
    function val = get.outIonization(obj),          val         =   obj.oIon;   end
    % oWfcn ~~~~~~~~~~~~
    function set.outWaveFunction(obj,val),          obj.oWfcn   =   val;        end
    function val = get.outWaveFunction(obj),        val         =   obj.oWfcn;  end
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
    if isempty(obj.sCI),        obj.sCI     =   false;                          end
    if isempty(obj.sCIF),       obj.sCIF    =   'QMolGrid--TDCI--CI';           end
    if isempty(obj.sDip),       obj.sDip    =   false;                          end
    if isempty(obj.sVel),       obj.sVel    =   false;                          end
    if isempty(obj.sVelT),      obj.sVelT   =   'dipole';                       end
    if isempty(obj.sAcc),       obj.sAcc    =   false;                          end
    if isempty(obj.sAccT),      obj.sAccT   =   'dipole';                       end
    if isempty(obj.sECI),       obj.sECI    =   false;                          end
    if isempty(obj.sEWfcn),     obj.sEWfcn  =   false;                          end
    if isempty(obj.sEF),        obj.sEF     =   false;                          end
    if isempty(obj.sIon),       obj.sIon    =   false;                          end
    if isempty(obj.sWfcn),      obj.sWfcn   =   false;                          end
    if isempty(obj.sRestF),     obj.sRestF  =   'QMolGrid--TDCI--Restart.mat';  end
    if isempty(obj.sRest),      obj.sRest   =   false;                          end
end
% Methods in separate files
    initialize(obj,CI,opt)
end
methods (Access=protected)
function ind = getOutputIndex(obj,t)
%getOutputIndex convert sampling time(s) to time indexes
    
    % Sempling indexes
    if isempty(t),              ind     =   obj.iref;                       % Default time sampling
    elseif ischar(t)
        if strcmpi(t,'all'),    ind     =   1:numel(obj.tspan);             % All time steps
        else,                   ind     =   obj.iref;
            warning('QMol:TDCI:outputIndex',['Unknown sampling option ' t '; default sampling used instead']);
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
    setOutputCI(obj,opt)
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
    
    ClassName           =   'QMol_TDCI';

    PropNames           =  {... Time ======================================
                            'disp','display',                       'T','time',                                         'dt','timeStep', ...
                            ... Gobbler ===================================
                            'ABC','absorbingBoundary',                      ...
                            ... External field ============================
                            'EF','externalField',                           ...
                            ... Save CI object into individual files =====
                            'sCI','saveCI',                             'sCIF','saveCIFileName',                            'sCIT','saveCITime', ...
                            ... Save dipole/velocity/acceleration =========
                            'sDip','saveDipole',                        'sDipI','saveDipoleWaveFunctionIndex',              'sDipT','saveDipoleTime', ...
                            'sVel','saveDipoleVelocity',                'sVelI','saveDipoleVelocityWaveFunctionIndex',      'sVelT','saveDipoleVelocityTime', ...
                            'sAcc','saveDipoleAcceleration',            'sAccI','saveDipoleAccelerationWaveFunctionIndex',  'sAccT','saveDipoleAccelerationTime', ...
                            ... Save energies =============================
                            'sECI','saveEnergyCI',                      'sECIT','saveEnergyCITime', ...
                            'sEWfcn','saveEnergyWaveFunction',          'sEWfcnT','saveEnergyWaveFunctionTime', ...
                            ... Save external field =======================
                            'sEF','saveExternalField', ...
                            ... Save ionization ===========================
                            'sIon','saveIonization',                    'sIWfcnI','saveIonizationWaveFunctionIndex',        'sIonT','saveIonizationTime', ...
                            ... Save wave function ========================
                            'sWfcn','saveWaveFunction',                 'sWfcnT','saveWaveFunctionTime',                    'sWfcnI','saveWaveFunctionIndex', ...
                            ... Save output functions =====================
                            'sF','saveOutputFunction',                  'sFT','saveOutputFunctionTime', ...
                            ... Save restart data =========================
                            'sRest','saveRestart',                      'sRestF','saveRestartTime',                         'sRestT','saveRestartFileName', ...
                            ... Output results ============================
                            'oCI','outCI', ...
                            'oDip','outDipole',                         'oVel','outDipoleVelocity',                         'oAcc','outDipoleAcceleration', ...
                            'oECI','outEnergyCI',                       'oEWfcn','outEnergyWaveFunction', ...
                            'oIon','outIonization', ...
                            'oWfcn','outWaveFunction', ...
                            'oF','outOutputFunction', ...
                            'oRest','outRestart'};
    
end
end
%% Propagate TDDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
    propagate(obj,opt)
end
methods (Access=protected)
    saveOutputCI(obj,t)
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

