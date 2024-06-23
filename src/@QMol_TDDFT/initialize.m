function initialize(obj,DFT,isRst)
%initialize initializes the TDDFT object and all its components
    
    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nargin < 3,      isRst   =   false;
    if nargin < 2,      DFT     =   [];         end, end

    if ~isempty(DFT),   obj.DFT =   DFT;        end

    % Set output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.isRun
        % External field 
        if ~isempty(obj.EF),    obj.EF.initialize(obj.DFT.disc);            end

        % Gobbler
        if ~isempty(obj.ABC),   obj.ABC.initialize(obj.DFT,obj.dt>0);       end

        % Restart mode
        if isRst,               obj.initializeChildren(true);               else

        % Initialization/house keeping
        obj.restart     =   [];

        % Time sampling
        obj.tspan       =   uniquetol([obj.T(1):obj.dt:obj.T(end), obj.T(end)],1e-10,'DataScale',1);
        if obj.dt < 0,      obj.tspan   =   flip(obj.tspan);                end
        obj.iref        =   obj.getOutputIndex(obj.T);

        % Implementation-specific initialization
        obj.initializeChildren(false);

        % Output data (will perform necessary cleanup)
        obj.setOutputDFT('init');                                           % DFT object in separate files
        obj.setOutputDipole('init');                                        % dipole/velocity/acceleration
        obj.setOutputEnergy('init');                                        % DFT and orbital energies
        obj.setOutputIonization('init');                                    % Ionization
        obj.setOutputOrbitalDensity('init');                                % Orbital/density
        obj.setOutputFunction('init');                                      % Output function
        obj.setOutputRestart('init');                                       % Restart file

        obj.setOutputChildren('init');                                      % Implementation-specific output
        end
    end

    % House keeping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.isInit       =   true;                                              % Not sure how useful this is ...
end