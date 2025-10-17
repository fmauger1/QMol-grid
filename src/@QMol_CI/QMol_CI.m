classdef (Abstract) QMol_CI < QMol_suite
%QMol_CI common interface for configuration interaction calculations
%   QMol_CI is an abstract class that define the common interface and
%   features for configuration interaction (CI) type calculations. End
%   users use one of the overloadping classes for CI calculations:
%     > QMol_CI_conv, which uses an explicit convolution to calculate the
%       two-electron integrals in the CI matrix (for 1D systems only)
%
%
%   See also QMol_CI_conv

%   Note: this is the QMol-grid implementation of an algorithm developed
%     by Giorgio Visentin (for Psi4)
%
%   Version     Date        Author
%   01.22.000   05/20/2025  F. Mauger
%       Prepare release version
%   01.23.001   06/05/2025  F. Mauger
%       Fix computeGroundState when the CI matrix is sparse
%   01.23.002   06/10/2025  F. Mauger
%       Fix implicit definition of Vext as a pseudopotential, cell of 
%       pseudopotentials, function handle, or vector

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.002','06/10/2025','G. Visentin & F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_CI})
function showInfo
    fprintf('  * QMol_CI:\n      > Configuration interaction (CI)\n'); 
    QMol_CI.version;
end
end
methods (Access=public)
function showDocumentation(obj)
%showDocumentation displays the run-time documentation
%   The run-time documentation prints out the specific configuration of the
%   CI object, as defined by its editable properties. The CI object must be
%   initialized before using showDocumentation.
%
%   obj.showDocumentation;
    
    % Header
    QMol_doc.showHeader;

    % Theory level
    QMol_doc.showSection('Theory');
    
    fprintf('  * Configuration interaction (CI) [Visentin 2025]\n');
    if obj.tol > 0,     fprintf('    Tolerance = %7.3g\n',obj.tol);             end
    if obj.isDip,       fprintf('    Including the dipole-coupling matrix\n');  end
    obj.version;

    ref                 =   obj.showDoc;                                    %#ok<PROP>
    fprintf('\n')
    
    ref                 =   [ref,{'Visentin 2025'}];                        %#ok<PROP>

    % Discretization
    QMol_doc.showSection('Discretization');

    if isempty(obj.disc)
        fprintf('  * No discretization specified\n')
    else
        ref             =   [ref, obj.disc.showDocumentation];              %#ok<PROP>
    end
    fprintf('\n')

    % Electronic structure
    QMol_doc.showSection('Electronic structure');

    fprintf('  * %i electrons\n',obj.N);

    if isempty(obj.SOB)
        fprintf('  * No spatial orbital specified\n')
    else
        fprintf('  * %i spatial orbitals\n',size(obj.SOB,2));
    end

    if obj.isSetConfigBasis
        obj.showConfigBasisSet
    elseif isempty(obj.CSB)
        fprintf('  * No configuration state basis specified\n')
    else
        fprintf('  * %i configuration-state basis elements\n',size(obj.CSB,1));
    end

    fprintf('\n')

    % Potentials
    QMol_doc.showSection('Potentials');

    if isempty(obj.Vext)
        fprintf('  * No external (atomic or molecular) potential specified\n')
    else
        ref             =   [ref, obj.Vext.showDocumentation];              %#ok<PROP>
    end
    fprintf('\n')

    fprintf('  * Electron-electron interaction potential           explicit convolution\n');
    fprintf('    Interaction pot. = %s (elec.-elec.)\n',func2str(obj.Vee));
    fprintf('\n')

    % Bibliography
    QMol_doc.showBibliography(ref);                                         %#ok<PROP>

    % Funding
    QMol_doc.showFunding;

    % Footer
    QMol_doc.showFooter;
end
% Methods in separate files
    mem = getMemoryProfile(obj,opt)
end
methods (Abstract,Access=?QMol_CI)
    ref = showDoc(obj)                  % implementation specific part of the run-time documentation
end
methods (Access=private)
function showConfigBasisSet(obj)
%showConfigBasisSet displays the configuration basis setup

    % Initialization
    fprintf('  * Configuration-state basis\n');

    [algo,~,locRef,locAct,locFrz] = obj.getConfigurationBasisParameters(false);
    nbRef               =   size(locRef,1);
    totAct              =   locAct{1};                                      for k = 2:nbRef
    totAct              =  [totAct,locAct{k}];                              end %#ok<AGROW>
    totAct              =   unique(totAct);

    % Display the configuration-state basis
    fprintf('    Type         = '); if nbRef > 1, fprintf('MR-'); end, fprintf('%s\n',obj.type)
    fprintf('    Total spin   = %g\n', 0.5*sum(sign(locRef)));                              if algo > 0
    fprintf('    Reference(s) =');   fprintf(' %i',locRef(1,:));  fprintf('\n');    for k = 2:nbRef
    fprintf('                  ');   fprintf(' %i',locRef(k,:));  fprintf('\n');    end,    end,    if ~isempty(locFrz)
    fprintf('    Frozen       =');   fprintf(' %i',locFrz(locFrz > 0));    fprintf('\n');
    fprintf('                  ');   fprintf(' %i',locFrz(locFrz < 0));    fprintf('\n');           end
    fprintf('    Active       =');   fprintf(' %i',totAct(totAct > 0));    fprintf('\n');
    fprintf('                  ');   fprintf(' %i',totAct(totAct < 0));    fprintf('\n');
    
    % Constraints
    if ~isempty(obj.noDbl) || ~isempty(obj.noEmpt),             fprintf('    Constraints on spatial orbitals:');    end
    if ~isempty(obj.noDbl),     fprintf('    > Not full  =');   fprintf(' %i',obj.noDbl);   fprintf('\n');          end
    if ~isempty(obj.noEmpt),    fprintf('    > Not empty =');   fprintf(' %i',obj.noEmpt);  fprintf('\n');          end

    % Size of configuration state basis
    fprintf('    %i configuration-state basis elements\n',size(obj.CSB,1));

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Constant,Access=public)
    isSpinPol           =   true    % By construction, all CI models are spin polarized
    dim                 =   1       % Model dimensionn (of the spin orbitals)
end
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Discretization
    disc                    % discretization object
    % Electronic structure
    N                       % number of electrons
    SOB                     % spatial orbital basis
    CSB                     % configuration state basis
    wfcn                    % wave functions
    % Potential
    Vext                    % external (atomic or molecular) potential
    Vee                 =   @(x) 1./sqrt(x.^2+2)
    % CI calculation
    CI                      % CI matrix
    DX                      % dipole-couping mmatrix
    isDip               =   true    % whether to calculate the dipole-coupling matrices
    disp                =   true
    tol                 =   1e-10
    ref                     % reference for CIS/D/...
    act                     % active spin orbitals
    frz                     % frozen spin orbitals
    noDbl                   % forbidden doubly occupied spatial orbitals
    noEmpt                  % forbidden empty orbitals
end
properties (GetAccess=public,SetAccess=?QMol_suite)
    type                    % type of CI calculation [ 'CIS' (default) | 'CISD' | 'RAS' ]
    isSetConfigBasis    =   false   % whether the configuration basis has been built from the setConfigurationBasis method [ true | false (default) ]
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    discretization          % (disc) domain discretization object [ QMol_disc (default []) ]
    numberElectron          % (N) number of electrons [ integer (default []) ]
    orbitalBasis            % (SOB) spatial orbital basis [ matrix (defaut []) ]
    configurationBasis      % (CSB) configuration basis [ matrix (default []) ]
    waveFunction            % (wfcn) wave function [ [] (default) | column vector | matrix ]
    externalPotential       % (Vext) external (atomic or molecular) potential [ QMol_DFT_Vext | pseudopotential | cell of pseudopotentials | function handle | vector (default []) ]
    interactionPotential    % (Vee) electron-electron interaction potential [ function handle (@(x) 1./sqrt(x.^2+2)) ]
    display                 % (disp) whether to display the progress of CI calculations [ true (default) | false ]
    tolerance               % (tol) threshold for setting CI matrix elements to zero [ scalar (default 1e-10) ]
    CImatrix                % (CI) CI matrix [ matrix (default []) ]
    DXmatrix                % (DX) dipole-coupling matrix [ matrix (default []) ]
    isDipole                % (isDip) whether to calculate the dipole-coupling matrix together with the CI one [ true (default) | false ]
    reference               % (ref) configuration state reference [ vector (default []) ]
    active                  % (act) list of spin orbitals that may contain an electron [ vector (default []) ]
    frozen                  % (frz) list of spin orbitals that must contain an electron [ vector (default []) ]
    noDouble                % (noDbl) list of spatial orbitals that must not contain two electons [ vector (default) ]
    noEmpty                 % (noEmp) list of spatial orbital that must not be empty [ vector (default) ]
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=?QMol_CI)
    isReal                  % whether the spatial orbitals are real valued
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % disc ~~~~~~~~~~~~~
    function set.discretization(obj,val),                   obj.disc    =   val;        end
    function val = get.discretization(obj),                 val         =   obj.disc;   end
    % N ~~~~~~~~~~~~~~~~
    function set.numberElectron(obj,val),                   obj.N       =   val;        end
    function val = get.numberElectron(obj),                 val         =   obj.N;      end
    % SOB ~~~~~~~~~~~~~~
    function set.orbitalBasis(obj,val),                     obj.SOB     =   val;        end
    function val = get.orbitalBasis(obj),                   val         =   obj.SOB;    end
    % CSB ~~~~~~~~~~~~~~
    function set.configurationBasis(obj,val),               obj.CSB     =   val;        end
    function val = get.configurationBasis(obj),             val         =   obj.CSB;    end
    % wfcn ~~~~~~~~~~~~~
    function set.waveFunction(obj,val),                     obj.wfcn    =   val;        end
    function val = get.waveFunction(obj),                   val         =   obj.wfcn;   end
    % Vext ~~~~~~~~~~~~~
    function set.externalPotential(obj,val),                obj.Vext    =   val;        end
    function val = get.externalPotential(obj),              val         =   obj.Vext;   end
    % Vee ~~~~~~~~~~~~~~
    function set.interactionPotential(obj,val),             obj.Vee     =   val;        end
    function val = get.interactionPotential(obj),           val         =   obj.Vee;    end
    % disp ~~~~~~~~~~~~~
    function set.display(obj,val),                          obj.disp    =   val;        end
    function val = get.display(obj),                        val         =   obj.disp;   end
    % tol ~~~~~~~~~~~~~~
    function set.tolerance(obj,val),                        obj.tol     =   val;        end
    function val = get.tolerance(obj),                      val         =   obj.tol;    end
    % CI ~~~~~~~~~~~~~~~
    function set.CImatrix(obj,val),                         obj.CI      =   val;        end
    function val = get.CImatrix(obj),                       val         =   obj.CI;     end
    % DX ~~~~~~~~~~~~~~~
    function set.DXmatrix(obj,val),                         obj.DX      =   val;        end
    function val = get.DXmatrix(obj),                       val         =   obj.DX;     end
    % isDip ~~~~~~~~~~~~
    function set.isDipole(obj,val),                         obj.isDip   =   val;        end
    function val = get.isDipole(obj),                       val         =   obj.isDip;  end
    % ref ~~~~~~~~~~~~~~
    function set.reference(obj,val),                        obj.ref     =   val;        end
    function val = get.reference(obj),                      val         =   obj.ref;    end
    % act ~~~~~~~~~~~~~~
    function set.active(obj,val),                           obj.act     =   val;        end
    function val = get.active(obj),                         val         =   obj.act;    end
    % frz ~~~~~~~~~~~~~~
    function set.frozen(obj,val),                           obj.frz     =   val;        end
    function val = get.frozen(obj),                         val         =   obj.frz;    end
    % noDbl ~~~~~~~~~~~~
    function set.noDouble(obj,val),                         obj.noDbl   =   val;        end
    function val = get.noDouble(obj),                       val         =   obj.noDbl;  end
    % noEmpt ~~~~~~~~~~~
    function set.noEmpty(obj,val),                          obj.noEmpt  =   val;        end
    function val = get.noEmpty(obj),                        val         =   obj.noEmpt; end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function obj = QMol_CI(varargin)
%QMol_CI constructor for configuration interaction models
%   Specify name properties with their corresponding values. Several name-
%   value pairs can be specified consecutively, and names are case 
%   insensitive
%
%   CI = QMol_CI('name1','value1',___) creates a CI object with the name
%     properties set to the specified values.
%
%   CI = QMol_CI(DFT,'name1','value1',___), DFT is a spin-restricted DFT 
%     object, reads off the domain, number of electrons, orbital basis,
%     atomic or molecular potential, and electron-electron interaction
%     potential off DFT. Additional or overwriting parameters can be
%     specified as a list of name-value pair arguments.
    
    % If any property input, pass them to the set member method
    obj.set(varargin{:});
end
function set(obj,varargin)
%set updates the name properties to the specified values
%   Several name-value pairs can be specified consecutively and names are
%   case insensitive
%       obj.set('name',value)                   
%       obj.set(name1,value1,name2,value2,___)
%   Read the domain, number of electrons, orbital basis, atomic or 
%   molecular potential, and electron-electron interaction potential off a 
%   spin-restricted DFT object by passing it as the first argument
%       obj.set(DFT,name1,value1,___)
%   The following name-value pairs specify additional or overwriting
%   parameters.

    % Anything to set
    if nargin == 1,     obj.reset;                                          else %#ok<ALIGN>

    if isa(varargin{1},'QMol_DFT_spinRes')
        % Initialization
        DFT             =   varargin{1};

        % Copy domain
        obj.disc        =   DFT.getDiscCopy;
        obj.disc.QM     =   obj;

        % Copy number of electrons
        obj.N           =   round(sum(DFT.occ));                            if abs(obj.N-sum(DFT.occ)) > 1e-10
        warning('QMol:CI:numberElectron', ...
            ['Fractional number of electrons in the DFT object\n' ...
             'Number of electrons rounded to ' num2str(obj.N)]);            end
        
        % Copy orbital basis
        if DFT.disc.isBasis
            % Reconstruct the orbitals
            KSO         =   DFT.disc.DFT_reconstructOrbital(DFT.KSO);
            obj.SOB     =   KSO.orbital;
        else
            % Copy the orbitals
            obj.SOB     =   DFT.orbital.orbital;
        end

        % Copy the atomic or molecular potential
        obj.Vext        =   DFT.Vext.copyPotential;

        % Copy electron-electron interaction potential (from Hartree)
        obj.Vee         =   DFT.Vh.Vee;

        % Call parent set method on name-value pairs
        set@QMol_suite(obj,varargin{2:end});

    else
        % Call parent set method on all parameters
        set@QMol_suite(obj,varargin{:});

    end,                                                                    end
end
function reset(obj)
%reset resets all temporary (transient) properties of the object

    % Also reset all components 
    % (cross-dependencies might be broken);
    if ~isempty(obj.disc),  obj.disc.reset;     end
    if ~isempty(obj.Vext),  obj.Vext.reset;     end
    
    % Initialization status
    obj.isInit          =   false;
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.Vee),        obj.Vee         =   @(x) 1./sqrt(x.^2+2);   end
    if isempty(obj.disp),       obj.disp        =   true;                   end
    if isempty(obj.tol),        obj.tol         =   1e-10;                  end
    if isempty(obj.isDip),      obj.isDip       =   true;                   end
    if isempty(obj.type),       obj.type        =   'CIS';                  end
end
% Methods in separate files
    initialize(obj,showBuild)
    setConfigurationBasis(obj)
    cleanConfigurationBasis(obj)
end
methods (Abstract,Access=?QMol_CI)
    isOK = initializeVee(obj)
end
methods (Access=private)
    [algo,locN,locRef,locAct,locFrz] = getConfigurationBasisParameters(obj,isBld)
end
methods (Hidden,Access = ?QMol_CI)
function disc = getDiscCopy(obj)
%getDiscCopy returns a streamlined copy of the discretization (object)
    
    % Create discretization object
    if isempty(obj.disc)
        % Get copy from x (implicit definition)
        disc            =   QMol_disc('x',obj.x);
    elseif obj.disc.isBasis
        % Basis-set model
        disc            =   QMol_disc_basis('x',obj.disc.x,'v',obj.disc.v);
    else
        % Regular grid discretization
        disc            =   QMol_disc('x',obj.disc.x);
    end

    % Set QM component
    disc.QM             =   obj;

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_CI_conv';
    PropNames           =  {'disc','N','SOB','CSB','wfcn','Vext','Vee','disp','tol','CI','DX','isDip',...
                            'ref','act','frz','noDbl','noEmp',...
                            'discretization','numberElectron','orbitalBasis','waveFunction', ...
                            'configurationBasis','externalPotential','interactionPotential', ...
                            'display','tolerance','CImatrix','DXmatrix','isDipole', ...
                            'type','reference','active','frozen','noDouble','noEmpty'};
end
end
%% CI calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function sparsify(obj)
%sparsify converts the CI and dipole-coupling into sparse matrices
%
%   obj.sparsify

    if ~isempty(obj.CI) && ~issparse(obj.CI),   obj.CI  =   sparse(obj.CI); end
    if ~isempty(obj.DX) && ~issparse(obj.DX),   obj.DX  =   sparse(obj.DX); end
end
    % Defined in separate file
    CI = computeCImatrix(obj)
    computeGroundState(obj,n)
    [E,DE] = getEnergy(obj)
    showEnergy(obj)
    analyzeCIspectrum(obj,n)
end
%% Core and two-electron integrals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Abstract,Access=?QMol_CI)
    h = getCoreIntegral(obj,k,l)
    v = getTwoElectronIntegral(obj,k,l,m,n)
    d = getDipoleIntegral(obj,k,l)
end
methods (Access=?QMol_CI)
function [a,r,s] = getExcitationIndexes(obj,CS1,CS2)
%getExcitationIndexes returns the excitation indexes and sign
%   The method identifies the indexes a in |CS1> and r in |CS2> such that
%   the configuration state |CS2> = s*|CS1_a^r>, with the scalar s (+/-1)
%   corresponding to the sign of the permutation. If the two configuration
%   states are identical or differ by more than two excitations, the output
%   are all set to empty vectors.

    % Initialization
    a                   =   [];
    r                   =   [];
    s                   =   [];

    % Excitation level
    [C,ic]              =   intersect(CS1,CS2);
    if numel(C) == obj.N || numel(C) < obj.N-2,     return,                 end

    % Signature of common components
    s                   =   1;
    for k = ic'
        n               =   find(CS2 == CS1(k));                            if n~=k % where is the matching index in CS2
        CS2([n,k])      =   CS2([k,n]);             s       =   -1*s;       end     % align indexes and update signature
    end

    if ~all(CS1(ic) == CS2(ic))
        warning('QMol:CI:getExcitationIndexes','Excitation indexes are not properly identified. Contact developers\n')
    end

    % Excitation indexes
    ie                  =   setdiff(1:obj.N,ic);
    a                   =   CS1(ie);
    r                   =   CS2(ie);

end
end
%% One-body density and one-particle density matrices %%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
    [rUp,rDw] = getDensityMatrix(obj,psi)
end
methods (Abstract,Access=public)
    rho = getDensity(obj,varargin)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

