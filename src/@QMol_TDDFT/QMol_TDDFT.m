classdef (Abstract) QMol_TDDFT < QMol_suite
%QMol_TDDFT time propagation interface for DFT and HF models. It defines:
%   > (abstract) interface class
%   > overall time propagation workflow
%     (subclasses need only implement elementary-step for propagation and
%     one-the-fly analysis/computations)
%   > common components for TDDFT simulations

% QUESTIONS:
%   > do I want to add a a unique identifyer for each specific propagator

% NOTES:
%   > general template for definitions of what to save is
%     - saveWhat        =   activates saving of 'what' (e.g., orbital);
%                           generally a boolean (true to save)
%     - saveWhatTime    =   specifies the times at which 'what' should be
%                           saved. This should support:
%           vector      =   user-defined times at which to save
%           []          =   use the global time vector
%           'all'       =   save at every time step
%     - saveWhatIndex   =   if 'what' is an orbital-resolved quantity,
%                           specify the indexes of the orbitals to include
%                           in the saved quantity. This should support
%           vector      =   user-defined list of indexes
%           []          =   save result from all orbitals
%           'all'       =   save results from all orbitals
%     Then, the associated saved result is contained in the 'outWhat'
%     structure, which includes both the result and the time-sampling
%     vector (and the indexes of selected orbitals for orbital-resolved
%     outputs).

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
    fprintf('  * QMol_TDDFT:\n      > TDDFT interface\n'); 
    QMol_TDDFT.version;
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
    if opt, QMol_DFT_profiler.showMemoryFootprint('TDDFT propagator has an undefined memory footprint', 0,1); end
    mem                 =   0; 

end
function mem = getMemoryProfileOrbitalDensity(obj,opt)

    if nargin < 2,  opt     =   false;  end
    
    if (obj.sRho || obj.sKSO || obj.sKSOP)   &&   opt
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham orbitals and one-body density', 0,1);
        if obj.sRho,    QMol_DFT_profiler.showMemoryFootprint('One-body density ouput has an undefined memory footprint', 0,2); end
        if obj.sKSO,    QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham orbitals ouput has an undefined memory footprint', 0,2); end
        if obj.sKSOP,   QMol_DFT_profiler.showMemoryFootprint('Orbitals projection ouput has an undefined memory footprint', 0,2); end
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
    % Save DFT object into individual files
    sDFT                =   false
    sDFTF               =   'QMolGrid--TDDFT--DFT'
    sDFTT
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
    sEDFT               =   false
    sEDFTT
    sEKSO               =   false
    sEKSOT
    % Save external field
    sEF                 =   false
    % Save ionization
    sIon                =   false
    sIKSOI
    sIonT
    % Save orbitals/density
    sRho                =   false
    sRhoT
    sKSO                =   false
    sKSOI
    sKSOT
    sKSOP               =   false
    sKSOPB
    sKSOPI
    sKSOPT
    % Save output functions
    sFRho
    sFRhoT
    sFKSO
    sFKSOT
    % Save restart data
    sRest               =   false
    sRestF              =   'QMolGrid--TDDFT--Restart.mat'
    sRestT
    % Output results
    oDFT
    oDip
    oVel
    oAcc
    oEDFT
    oEKSO
    oIon
    oRho
    oKSO
    oKSOP
    oFRho
    oFKSO
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
    % Save DFT object into individual files
    saveDFT                             % sDFT
    saveDFTFileName                     % sDFTF
    saveDFTTime                         % sDFTT
    % Save dipole/velocity/acceleration =
    saveDipole                          % sDip
    saveDipoleOrbitalIndex              % sDipI
    saveDipoleTime                      % sDipT
    saveDipoleVelocity                  % sVel
    saveDipoleVelocityOrbitalIndex      % sVelI
    saveDipoleVelocityTime              % sVelT
    saveDipoleAcceleration              % sAcc
    saveDipoleAccelerationOrbitalIndex  % sAccI
    saveDipoleAccelerationTime          % sAccT
    % Save energies =====================
    saveEnergyDFT                       % sEDFT
    saveEnergyDFTTime                   % sEDFTT
    saveEnergyOrbital                   % sEKSO
    saveEnergyOrbitalTime               % sEKSOT
    % Save external field ===============
    saveExternalField                   % sEF
    % Save ionization ===================
    saveIonization                      % sIon
    saveIonizationOrbitalIndex          % sIKSOI
    saveIonizationTime                  % sIonT
    % Save orbital/density ==============
    saveDensity                         % sRho
    saveDensityTime                     % sRhoT
    saveOrbital                         % sKSO
    saveOrbitalIndex                    % sKSOI
    saveOrbitalTime                     % sKSOT
    saveOrbitalProjection               % sKSOP
    saveOrbitalProjectionBasis          % sKSOPB
    saveOrbitalProjectionIndex          % sKSOPI
    saveOrbitalProjectionTime           % sKSOPT
    % Save output functions =============
    saveOutputFunctionDensity           % sFRho
    saveOutputFunctionDensityTime       % sFRhoT
    saveOutputFunctionOrbital           % sFKSO
    saveOutputFunctionOrbitalTime       % sFKSOT
    % Save restart data =================
    saveRestart                         % sRest
    saveRestartFileName                 % sRestF
    saveRestartTime                     % sRestT
    % Output results ====================
    outDFT                              % oDFT
    outDipole                           % oDip
    outDipoleVelocity                   % oVel
    outDipoleAcceleration               % oAcc
    outEnergyDFT                        % oEDFT
    outEnergyOrbital                    % oEKSO
    outIonization                       % oIon
    outDensity                          % oRho
    outOrbital                          % oKSO
    outOrbitalProjection                % oKSOP
    outOutputFunctionDensity            % oFRho
    outOutputFunctionOrbital            % oFKSO
    outRestart                          % oRest
end
properties (Access=protected)
    % Time propagation
    DFT
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
    % sDFT ~~~~~~~~~~~~~                          ==== Save DFT to files ==========
    function set.saveDFT(obj,val),                  obj.sDFT    =   val;        end
    function val = get.saveDFT(obj),                val         =   obj.sDFT;   end
    % sDFTF ~~~~~~~~~~~~
    function set.saveDFTFileName(obj,val),          obj.sDFTF   =   val;        end
    function val = get.saveDFTFileName(obj),        val         =   obj.sDFTF;  end
    % sDFTT ~~~~~~~~~~~~
    function set.saveDFTTime(obj,val),              obj.sDFTT   =   val;        end
    function val = get.saveDFTTime(obj),            val         =   obj.sDFTT;  end
    % sDip ~~~~~~~~~~~~~                          ==== Save dipole/vel/acc ========
    function set.saveDipole(obj,val),               obj.sDip    =   val;        end
    function val = get.saveDipole(obj),             val         =   obj.sDip;   end
    % sDipI ~~~~~~~~~~~~
    function set.saveDipoleOrbitalIndex(obj,val),   obj.sDipI   =   val;        end
    function val = get.saveDipoleOrbitalIndex(obj), val         =   obj.sDipI;  end
    % sDipT ~~~~~~~~~~~~
    function set.saveDipoleTime(obj,val),           obj.sDipT   =   val;        end
    function val = get.saveDipoleTime(obj),         val         =   obj.sDipT;  end
    % sVel ~~~~~~~~~~~~~
    function set.saveDipoleVelocity(obj,val),       obj.sVel    =   val;        end
    function val = get.saveDipoleVelocity(obj),     val         =   obj.sVel;   end
    % sVelI ~~~~~~~~~~~~
    function set.saveDipoleVelocityOrbitalIndex(obj,val), obj.sVelI=val;        end
    function val = get.saveDipoleVelocityOrbitalIndex(obj), val =   obj.sVelI;  end
    % sVelT ~~~~~~~~~~~~
    function set.saveDipoleVelocityTime(obj,val),   obj.sVelT   =   val;        end
    function val = get.saveDipoleVelocityTime(obj), val         =   obj.sVelT;  end
    % sAcc ~~~~~~~~~~~~~
    function set.saveDipoleAcceleration(obj,val),   obj.sAcc    =   val;        end
    function val = get.saveDipoleAcceleration(obj), val         =   obj.sAcc;   end
    % sAccI ~~~~~~~~~~~~
    function set.saveDipoleAccelerationOrbitalIndex(obj,val), obj.sAccI=val;    end
    function val = get.saveDipoleAccelerationOrbitalIndex(obj), val=obj.sAccI;  end
    % sVelT ~~~~~~~~~~~~
    function set.saveDipoleAccelerationTime(obj,val), obj.sAccT =   val;        end
    function val = get.saveDipoleAccelerationTime(obj), val     =   obj.sAccT;  end
    % sEDFT ~~~~~~~~~~~~                          ==== Save energies ==============
    function set.saveEnergyDFT(obj,val),            obj.sEDFT   =   val;        end
    function val = get.saveEnergyDFT(obj),          val         =   obj.sEDFT;  end
    % sEDFTT ~~~~~~~~~~~
    function set.saveEnergyDFTTime(obj,val),        obj.sEDFTT  =   val;        end
    function val = get.saveEnergyDFTTime(obj),      val         =   obj.sEDFTT; end
    % sEKSO ~~~~~~~~~~~~
    function set.saveEnergyOrbital(obj,val),        obj.sEKSO   =   val;        end
    function val = get.saveEnergyOrbital(obj),      val         =   obj.sEKSO;  end
    % sEDFTT ~~~~~~~~~~~
    function set.saveEnergyOrbitalTime(obj,val),    obj.sEKSOT  =   val;        end
    function val = get.saveEnergyOrbitalTime(obj),  val         =   obj.sEKSOT; end
    % sEF ~~~~~~~~~~~~~~                          ==== Save external field ========
    function set.saveExternalField(obj,val),        obj.sEF     =   val;        end
    function val = get.saveExternalField(obj),      val         =   obj.sEF;    end
    % sIon ~~~~~~~~~~~~~                          ==== Save ionization ============
    function set.saveIonization(obj,val),           obj.sIon    =   val;        end
    function val = get.saveIonization(obj),         val         =   obj.sIon;   end
    % sIKSOI ~~~~~~~~~~~
    function set.saveIonizationOrbitalIndex(obj,val), obj.sIKSOI=   val;        end
    function val = get.saveIonizationOrbitalIndex(obj), val     =   obj.sIKSOI; end
    % sIonT ~~~~~~~~~~~~
    function set.saveIonizationTime(obj,val),       obj.sIonT   =   val;        end
    function val = get.saveIonizationTime(obj),     val         =   obj.sIonT;  end
    % sRho ~~~~~~~~~~~~~                          ==== Save orbital/density =======
    function set.saveDensity(obj,val),              obj.sRho    =   val;        end
    function val = get.saveDensity(obj),            val         =   obj.sRho;   end
    % sRhoT ~~~~~~~~~~~~
    function set.saveDensityTime(obj,val),          obj.sRhoT   =   val;        end
    function val = get.saveDensityTime(obj),        val         =   obj.sRhoT;  end
    % sKSO ~~~~~~~~~~~~~
    function set.saveOrbital(obj,val),              obj.sKSO    =   val;        end
    function val = get.saveOrbital(obj),            val         =   obj.sKSO;   end
    % sKSOI ~~~~~~~~~~~~
    function set.saveOrbitalIndex(obj,val),         obj.sKSOI   =   val;        end
    function val = get.saveOrbitalIndex(obj),       val         =   obj.sKSOI;  end
    % sKSOT ~~~~~~~~~~~~
    function set.saveOrbitalTime(obj,val),          obj.sKSOT   =   val;        end
    function val = get.saveOrbitalTime(obj),        val         =   obj.sKSOT;  end
    % sKSOP ~~~~~~~~~~~~
    function set.saveOrbitalProjection(obj,val),    obj.sKSOP   =   val;        end
    function val = get.saveOrbitalProjection(obj),  val         =   obj.sKSOP;  end
    % sKSOPB ~~~~~~~~~~~
    function set.saveOrbitalProjectionBasis(obj,val), obj.sKSOPB=   val;        end
    function val = get.saveOrbitalProjectionBasis(obj), val     =   obj.sKSOPB; end
    % sKSOPI ~~~~~~~~~~~
    function set.saveOrbitalProjectionIndex(obj,val), obj.sKSOPI=   val;        end
    function val = get.saveOrbitalProjectionIndex(obj), val     =   obj.sKSOPI; end
    % sKSOPT ~~~~~~~~~~~
    function set.saveOrbitalProjectionTime(obj,val), obj.sKSOPT =   val;        end
    function val = get.saveOrbitalProjectionTime(obj), val      =   obj.sKSOPT; end
    % sFRho ~~~~~~~~~~~~                          ==== Save output function =======
    function set.saveOutputFunctionDensity(obj,val), obj.sFRho  =   val;        end
    function val = get.saveOutputFunctionDensity(obj), val      =   obj.sFRho;  end
    % sFRhoT ~~~~~~~~~~~
    function set.saveOutputFunctionDensityTime(obj,val), obj.sFRhoT=val;        end
    function val = get.saveOutputFunctionDensityTime(obj), val  =   obj.sFRhoT; end
    % sFKSO ~~~~~~~~~~~~
    function set.saveOutputFunctionOrbital(obj,val), obj.sFKSO  =   val;        end
    function val = get.saveOutputFunctionOrbital(obj), val      =   obj.sFKSO;  end
    % sFRhoT ~~~~~~~~~~~
    function set.saveOutputFunctionOrbitalTime(obj,val), obj.sFKSOT=val;        end
    function val = get.saveOutputFunctionOrbitalTime(obj), val  =   obj.sFKSOT; end
    % sRest ~~~~~~~~~~~~                          ==== Save restart data ==========
    function set.saveRestart(obj,val),              obj.sRest   =   val;        end
    function val = get.saveRestart(obj),            val         =   obj.sRest;  end
    % sRestF ~~~~~~~~~~~
    function set.saveRestartFileName(obj,val),      obj.sRestF  =   val;        end
    function val = get.saveRestartFileName(obj),    val         =   obj.sRestF; end
    % sRestT ~~~~~~~~~~~
    function set.saveRestartTime(obj,val),          obj.sRestT  =   val;        end
    function val = get.saveRestartTime(obj),        val         =   obj.sRestT; end
    % oDFT ~~~~~~~~~~~~~                         ==== Output results ==============
    function set.outDFT(obj,val),                   obj.oDFT    =   val;        end
    function val = get.outDFT(obj),                 val         =   obj.oDFT;   end
    % oDip ~~~~~~~~~~~~~
    function set.outDipole(obj,val),                obj.oDip    =   val;        end
    function val = get.outDipole(obj),              val         =   obj.oDip;   end
    % oVel ~~~~~~~~~~~~~
    function set.outDipoleVelocity(obj,val),        obj.oVel    =   val;        end
    function val = get.outDipoleVelocity(obj),      val         =   obj.oVel;   end
    % oAcc ~~~~~~~~~~~~~
    function set.outDipoleAcceleration(obj,val),    obj.oAcc    =   val;        end
    function val = get.outDipoleAcceleration(obj),  val         =   obj.oAcc;   end
    % oEDFT ~~~~~~~~~~~~
    function set.outEnergyDFT(obj,val),             obj.oEDFT   =   val;        end
    function val = get.outEnergyDFT(obj),           val         =   obj.oEDFT;  end
    % oEKSO ~~~~~~~~~~~~
    function set.outEnergyOrbital(obj,val),         obj.oEKSO   =   val;        end
    function val = get.outEnergyOrbital(obj),       val         =   obj.oEKSO;  end
    % oIon ~~~~~~~~~~~~~
    function set.outIonization(obj,val),            obj.oIon    =   val;        end
    function val = get.outIonization(obj),          val         =   obj.oIon;   end
    % oRho ~~~~~~~~~~~~~
    function set.outDensity(obj,val),               obj.oRho    =   val;        end
    function val = get.outDensity(obj),             val         =   obj.oRho;   end
    % oKSO ~~~~~~~~~~~~~
    function set.outOrbital(obj,val),               obj.oKSO    =   val;        end
    function val = get.outOrbital(obj),             val         =   obj.oKSO;   end
    % oKSOP ~~~~~~~~~~~~
    function set.outOrbitalProjection(obj,val),     obj.oKSOP   =   val;        end
    function val = get.outOrbitalProjection(obj),   val         =   obj.oKSOP;  end
    % oFRho ~~~~~~~~~~~~
    function set.outOutputFunctionDensity(obj,val), obj.oFRho   =   val;        end
    function val = get.outOutputFunctionDensity(obj),val        =   obj.oFRho;  end
    % oFRho ~~~~~~~~~~~~
    function set.outOutputFunctionOrbital(obj,val), obj.oFKSO   =   val;        end
    function val = get.outOutputFunctionOrbital(obj),val        =   obj.oFKSO;  end
    % oFRho ~~~~~~~~~~~~
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
    if isempty(obj.sDFT),       obj.sDFT    =   false;                          end
    if isempty(obj.sDFTF),      obj.sDFTF   =   'QMolGrid--TDDFT--DFT';         end
    if isempty(obj.sDip),       obj.sDip    =   false;                          end
    if isempty(obj.sVel),       obj.sVel    =   false;                          end
    if isempty(obj.sVelT),      obj.sVelT   =   'dipole';                       end
    if isempty(obj.sAcc),       obj.sAcc    =   false;                          end
    if isempty(obj.sAccT),      obj.sAccT   =   'dipole';                       end
    if isempty(obj.sEDFT),      obj.sEDFT   =   false;                          end
    if isempty(obj.sEKSO),      obj.sEKSO   =   false;                          end
    if isempty(obj.sEF),        obj.sEF     =   false;                          end
    if isempty(obj.sIon),       obj.sIon    =   false;                          end
    if isempty(obj.sRho),       obj.sRho    =   false;                          end
    if isempty(obj.sKSO),       obj.sKSO    =   false;                          end
    if isempty(obj.sKSOP),      obj.sKSOP   =   false;                          end
    if isempty(obj.sRestF),     obj.sRestF  =   'QMolGrid--TDDFT--Restart.mat'; end
    if isempty(obj.sRest),      obj.sRest   =   false;                          end
end
% Methods in separate files
    initialize(obj,DFT,opt)
end
methods (Access=protected)
function ind = getOutputIndex(obj,t)
%getOutputIndex convert sampling time(s) to time indexes
    
    % Sempling indexes
    if isempty(t),              ind     =   obj.iref;                       % Default time sampling
    elseif ischar(t)
        if strcmpi(t,'all'),    ind     =   1:numel(obj.tspan);             % All time steps
        else,                   ind     =   obj.iref;
            warning('QMol:TDDFT:outputIndex',['Unknown sampling option ' t'; default sampling used instead']);
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
    setOutputDFT(obj,opt)
    setOutputDipole(obj,opt)
    setOutputEnergy(obj,opt)
    setOutputFunction(obj,opt)
    setOutputIonization(obj,opt)
    setOutputRestart(obj,opt)
end
methods (Abstract,Access=protected)
    setOutputOrbitalDensity(obj,opt)
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
    
    ClassName           =   'QMol_TDDFT';
    PropNames           =  {... Time ======================================
                            'disp','display',                       'T','time',                                     'dt','timeStep', ...
                            ... Gobbler ===================================
                            'ABC','absorbingBoundary',                      ...
                            ... External field ============================
                            'EF','externalField',                           ...
                            ... Save DFT object into individual files =====
                            'sDFT','saveDFT',                       'sDFTF','saveDFTFileName',                      'sDFTT','saveDFTTime',                          ...
                            ... Save dipole/velocity/acceleration =========
                            'sDip','saveDipole',                    'sDipI','saveDipoleOrbitalIndex',               'sDipT','saveDipoleTime', ...
                            'sVel','saveDipoleVelocity',            'sVelI','saveDipoleVelocityOrbitalIndex',       'sVelT','saveDipoleVelocityTime', ...
                            'sAcc','saveDipoleAcceleration',        'sAccI','saveDipoleAccelerationOrbitalIndex',   'sAccT','saveDipoleAccelerationTime',...
                            ... Save energies =============================
                            'sEDFT','saveEnergyDFT',                'sEDFTT','saveEnergyDFTTime', ...
                            'sEKSO','saveEnergyOrbital',            'sEKSOT','saveEnergyOrbitalTime'...
                            ... Save external field =======================
                            'sEF','saveExternalField', ...
                            ... Save ionization ===========================
                            'sIon','saveIonization',                'sIKSOI','saveIonizationOrbitalIndex',          'sIonT','saveIonizationTime', ...
                            ... Save orbitals/density =====================
                            'sRho','saveDensity',                   'sRhoT','saveDensityTime',                      'sKSO','saveOrbital', ...
                            'sKSOT','saveOrbitalTime',              'sKSOI','saveOrbitalIndex', ...
                            'sKSOP','saveOrbitalProjection',        'sKSOPB','saveOrbitalProjectionBasis', ...
                            'sKSOPI','saveOrbitalProjectionIndex',  'sKSOPT','saveOrbitalProjectionTime',...
                            ... Save output functions =====================
                            'sFRho','saveOutputFunctionDensity',    'sFRhoT','saveOutputFunctionDensityTime', ...
                            'sFKSO','saveOutputFunctionOrbital',    'sFKSOT','saveOutputFunctionOrbitalTime', ...
                            ... Save restart data =========================
                            'sRest','saveRestart',                  'sRestF','saveRestartTime',                     'sRestT','saveRestartFileName', ...
                            ... Output results ============================
                            'oDFT','outDFT', ...
                            'oDip','outDipole',                     'oVel','outDipoleVelocity',                     'oAcc','outDipoleAcceleration', ...
                            'oEDFT','outEnergyDFT',                 'oEKSO','outEnergyOrbital', ...
                            'oIon','outIonization', ...
                            'oRho','outDensity',                    'oKSO','outOrbital',                            'oKSOP','outOrbitalProjection', ...
                            'oFRho','outOutputFunctionDensity'      'oFKSO','outOutputFunctionOrbital', ...
                            'oRest','outRestart'};
    
end
end
%% Propagate TDDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
    propagate(obj,opt)
end
methods (Access=protected)
    saveOutputDFT(obj,t)
    saveOutputEnergy(obj,k,t)
    saveOutputFunction(obj,k,t)
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
    saveOutputOrbitalDensity(obj,k,t);
end
methods (Access=protected)
    saveOutputResults(obj,k,t)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

