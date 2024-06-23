classdef (Abstract) QMol_SEq < QMol_suite
%QMol_QMol_SEq Schrodinger-equation (SE) level of theory. It defines
%   > (abstract) interface class
%   > common components for SE-based simulations
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_SEq})
function showInfo
    fprintf('  * QMol_SEq:\n      > SE interface\n'); 
    QMol_SEq.version;
end
end
methods (Access=public)
function showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Header
    QMol_doc.showHeader;

    % Theory level
    QMol_doc.showSection('Theory');
    ref                 =   obj.showDoc;

    fprintf('  * The variational SE model matches the canonical Hamiltonian formalism\n');
    fprintf('    describing the time evolution of the system [Mauger 2024].\n');
    obj.version;
    
    ref                 =   [ref, {'Mauger 2024'}];

    fprintf('\n');

    % Discretization
    QMol_doc.showSection('Discretization');

    if isempty(obj.disc)
        fprintf('  * No discretization specified\n')
    else
        ref             =   [ref, obj.disc.showDocumentation];
    end
    fprintf('\n')

    % Atomic/molecular geometry
    QMol_doc.showSection('Potential');
    
    if isempty(obj.V)
        fprintf('  * No potential specified\n')
    else
        ref             =   [ref, obj.V.showDocumentation];
    end

    fprintf('\n')

    % Electronic DOF
    QMol_doc.showSection('System model');

    fprintf('  * Electronic structure                                    wave functions\n');
    fprintf('    %i wave function(s)\n',obj.N);
    fprintf('\n')

    % Bibliography
    QMol_doc.showBibliography(ref);

    % Funding
    QMol_doc.showFunding;

    % Footer
    QMol_doc.showFooter;
end
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint of the SE object with all its components initialized and
%   used. The methods will try to return the various objects in the same
%   configuration they were passed, but it's not 100% guarantied. 
%   getMemoryProfile should be used in isolation and not, e.g., preceding 
%   ground-state or time propagation simulations.
    
    % Initialization
    if nargin < 2,  opt     =   false;  end

    % Domain discretization
    disc                =   obj.getDiscCopy;                                %#ok<PROPLC> 
    mem                 =   disc.getMemoryProfile(opt);                     %#ok<PROPLC> 

    inDisc              =   obj.disc;
    obj.disc            =   disc;                                           %#ok<PROPLC> 

    % Wave function(s)
    wfcn                =   disc.SE_allocateWaveFunction(0);                %#ok<PROPLC> 
    mem                 =   mem + wfcn.getMemoryProfile(opt,obj.N);         %#ok<PROPLC> 

    % Potential
    if ~isempty(obj.V)
        if ~obj.V.isInit,       obj.V.SE        =   obj;                    end
        mem             =   mem + obj.V.getMemoryProfile(opt);
        if ~obj.V.isInit,       obj.V.reset();                              end
    else
        obj.V           =   QMol_SE_V;          obj.V.SE        =   obj;
        mem             =   mem + obj.V.getMemoryProfile(opt);
        obj.V           =   [];
    end

    % Clean up
    obj.disc            =   inDisc;
    clear disc 
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Abstract,Constant,Access=public)
    dim
end
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Discretization
    disc
    % Electronic structure
    wfcn
    % Potential
    V
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    discretization          % disc
    waveFunction            % wfcn
    potential               % V
end
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    N                       % Number of wave functions (by default 1)
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    numberWaveFunction      % N
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % disc ~~~~~~~~~~~~~
    function set.discretization(obj,val),                   obj.disc    =   val;        end
    function val = get.discretization(obj),                 val         =   obj.disc;   end
    % wfcn ~~~~~~~~~~~~~
    function set.waveFunction(obj,val),                     obj.wfcn    =   val;        end
    function val = get.waveFunction(obj),                   val         =   obj.wfcn;   end
    % N ~~~~~~~~~~~~~~~~
    function val = get.N(obj),  if obj.isInit,              val         =   size(obj.wfcn.wfcn,obj.dim+1);      %#ok<ALIGN> % obj.N should be empty after the object initialization
                                elseif isempty(obj.N),      val         =   1;
                                else,                       val         =   obj.N;      end,    end
    function set.numberWaveFunction(obj,val),               obj.N       =   val;        end
    function val = get.numberWaveFunction(obj),             val         =   obj.N;      end
    % V ~~~~~~~~~~~~~~~~
    function set.potential(obj,val),                        obj.V       =   val;        end
    function val = get.potential(obj),                      val         =   obj.V;      end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object

    % Reset all components 
    % (cross-dependencies might be broken);
    if ~isempty(obj.disc),  obj.disc.reset;     end
    if ~isempty(obj.wfcn),  obj.wfcn.reset;     end
    if ~isempty(obj.V),     obj.V.reset;        end
    
    % Initialization status
    obj.isInit          =   false;
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

end
function initialize(obj,showBuild)
%initialize initializes the DFT object and all its components
    
    % Initialization
    if nargin == 1,         showBuild   =   [];     end
    if isempty(showBuild),  showBuild   =   0;      end

    function buildSection(txt)
        if showBuild > 0,   fprintf('  * %-66s',txt);   end
    end
    function buildStatus(t)
        if showBuild > 0,   if t, fprintf('  OK\n'); else, fprintf('FAIL\n'); end, end
    end

    if showBuild == 1
        QMol_doc.showHeader;
        QMol_doc.showFunding;
        QMol_doc.showSection('Building:');
    end

    obj.reset;
    
    % Discretization
    buildSection('Discretization');
    if isempty(obj.disc)
        warning('QMol:SE:noDiscretization','Missing the discretization component; initialization skipped on it.')
        buildStatus(false);
    else
        % Check if discretization has changed
        if ~obj.disc.isInit,        obj.reset;      end
        
        % Initialize the discretization component
        obj.disc.initialize(obj);
        buildStatus(obj.disc.isInit);
    end

    % Wave function(s)
    buildSection('Wave function(s)');
    if isempty(obj.disc)   ||   ~obj.disc.isInit
        warning('QMol:SE:wfcn',['Unable to initialized the wave function(s) (missing or uninitialized discretization component).\n'...
                                'No initialization performed on the wave function(s).'])
        buildStatus(false);
    else
        % Allocate wave function
        obj.wfcn        =   obj.disc.SE_allocateWaveFunction(obj.N,obj.wfcn);

        % Initialization status
        buildStatus(obj.wfcn.isInit);
    end

    % Potential
    buildSection('Potential');
    if isempty(obj.V)
        warning('QMol:SE:noPotential','Missing the potential component; initialization skipped on it.')
        buildStatus(false)
    else
        obj.V.initialize(obj);
        buildStatus(obj.V.isInit);
    end
    
    % House keeping
    obj.isInit      =   true;

    if showBuild > 0,   fprintf('\n'); 
    if showBuild == 1,  QMol_doc.showFooter;      end, end
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_SEq';
    PropNames           =  {'disc','wfcn','V','N',                          ...
                            'discretization','waveFunction','potential','numberWaveFunction'};
end
end
%% Interface methods (to be overloaded) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Abstract,Access=?QMol_SEq)
    % Display
    showDoc                 % Show implementation-specific documentation
end
%% Compute SE components (get) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function [E1,E2,E3] = getEnergy(obj,I1,I2)
%getEnergy computes the Schrodinger-equation or wave function (wfcn) 
%   energy components.
%
%   [E,err] = getEnergy('wave function') computes the energies E and errors
%       err for the member-property wfcns. Optionally, the potential object 
%       to be used for the computation of the energies and error can be 
%       included.
%
%   [Etot,Ekin,Epot] = getEnergy('Schrodinger equation') computes the 
%       total, kinetic, and potential energy components, respectively, for 
%       the member-property wfcns. 

    % Initialization 
    if nargin == 1,     opt     =   'se';   IC  =   [];                     %#ok<ALIGN> 
    elseif nargin == 2
        if ischar(I1),  opt     =   I1;     IC  =   [];
        else,           opt     =   'se';   IC  =   I1;                     end
    else
        if ischar(I1),  opt     =   I1;     IC  =   I2;
        else,           opt     =   I2;     IC  =   I1;                     end, end

    % Compute energy components
    switch lower(opt)
        % DFT energy -----------------------
        case {'se','schrodinger','schrodinger equation','schrodinger_equation'}

            % Get energy components
            E2          =   obj.disc.SE_energyKinetic(obj.wfcn);

            if ~isempty(obj.V),     E3  =   obj.V.getEnergy(obj.wfcn);
            else,                   E3  =   0;                              end
            
            % Total energy
            E1 = E2 + E3;
        % Orbital energies -----------------
        case {'wfcn','wave function','wave_function'}
            % Initialization
            if isempty(IC)
                Vse     =   obj.V;
            elseif isa(IC,obj.disc.SE_classPotential)
                Vse     =   IC;
            else
                Vse     =   obj.V;
                warning('QMol:SE:getEnergy', ...
                    'Expecting a potential object; use member potential instead.')
            end

            if nargout == 1
                E1      =   obj.disc.SE_energyWaveFunction(Vse,obj.wfcn);
            else
                [E1,E2] =   obj.disc.SE_energyWaveFunction(Vse,obj.wfcn);
            end
        % Unexpected case ------------------
        otherwise
            E1 = NaN;   E2 = NaN;   E3 = NaN;
            warning('QMol:SE:getEnergy', ...
                ['Unknown energy-computation case ' opt '; no energy computed.'])

    end
end
function showEnergy(obj,opt)
%showEnergy displays the Schrodinger-equation or wave function (wfcn) 
%   energy components.
%
%   showEnergy('wave function') displays the energies and errors for the 
%       member-property wfcn. 
%
%   showEnergy('Schrodinger equation') displays the total, kinetic, and 
%       potential energy components for the member-property wfcn. 
    
    % Defaul option
    if nargin < 2,      opt     =   [];         end
    if isempty(opt),    opt     =   'wfcn';     end

    % Show energy
    switch lower(opt)
        % DFT energy -----------------------
        case {'se','schrodinger','schrodinger equation','schrodinger_equation'}
            % Initialization
            [Etot,Ekin,Epot]    =   obj.getEnergy('SE');

            
            fprintf('  Component      Energy (a.u.)      Energy (eV)\n')
            fprintf('  -----------    -------------     -------------\n');
            
            % Show energy component
            fprintf('  Kinetic        %#10.3f         %#10.3f\n',Ekin,convertUnit.au2ev(Ekin));
            fprintf('  Potential      %#10.3f         %#10.3f\n',Epot,convertUnit.au2ev(Epot));
            fprintf('  -----------    -------------     -------------\n');
            fprintf('  Total          %#10.3f         %#10.3f\n',Etot,convertUnit.au2ev(Etot));

            % Finalization
            fprintf('  ----------------------------------------------\n');

        % Orbital energies -----------------
        case {'wfcn','wave function','wave_function'}
            % Initialization
            [E,err]     =   obj.getEnergy('wfcn');

            fprintf('  Wave fcn      Energy (-eV)         Error(a.u.)\n')
            fprintf('  --------     ------------          -----------\n');

            % Display results
            fprintf('    %3i        %9.3f             %#10.3e\n', ...
                [1:obj.N; convertUnit.au2ev(-E.'); err.']);

            % Finalization
            fprintf('  ----------------------------------------------\n');

        % Unexpected case ------------------
        otherwise
            warning('QMol:SE:showEnergy', ...
                ['Unknown energy-display case ' opt '; no energy shown.'])

    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

