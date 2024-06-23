classdef QMol_DFT_SCF_Anderson < QMol_suite
%QMol_DFT_SCF_Anderson Anderson's mixing scheme for density-functional
%   theory (DFT) self-consistent-field (SCF) calculations. It defines
%   > tol           Convergence tolerance
%                   AKA: tolerance
%   > maxit         Maximum number of iterations
%                   AKA: maxIterations
%   > mix           Mixing coefficient
%                   AKA: mixing
%   > mode          Mixing mode ('density' or 'potential')
%                   AKA: mixingMode
%   > test          Convergence test
%                   AKA: convergenceTest ('density', 'eigenvalues')
%   > disp          Display option
%                   AKA: display
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT_SCF_Anderson})
function showInfo
    fprintf( '  * QMol_DFT_SCF_Anderson:\n');
    fprintf(['      > DFT self-consistent field (SCF) solver\n' ...
             '      > Anderson''s mixing scheme\n']); 
    QMol_DFT_SCF_Anderson.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Header
    if ~obj.isRun
        QMol_doc.showHeader;
    end
    
    QMol_doc.showSection('Self-consistent-field (SCF) methods');

    % Eigensolver
    if ~isempty(obj.eigs),  ref =   obj.eigs.showDocumentation;
    else,                   ref =   {};                             end

    % Anderson's mixing parameters
    fprintf('  * Self-consistent-field (SCF)                          Anderson''s mixing\n');
    fprintf('    DFT-SCF solver using an Anderson''s mixing scheme [Anderson 1965], as\n');
    fprintf('    described in [Johnson 1988].\n');

    fprintf('    Tolerance   = %-6.2g\n', obj.tol);
    fprintf('    Max. iter.  = %i\n', obj.maxit);
    fprintf('    Mix. coeff. = %-6.2g\n', obj.mix);
    if ~isempty(obj.mixMode)
        if      obj.mixMode == 1,   fprintf('    Mix. mode   = density\n');     %#ok<ALIGN> 
        elseif  obj.mixMode == 2,   fprintf('    Mix. mode   = potential\n');
        else,                       fprintf('    Unexpected mixing mode\n');    end
    else
        fprintf('    Mix. mode   = %s\n', obj.mode);
    end
    if ~isempty(obj.testMode)
        if      obj.testMode == 1,  fprintf('    Conv. test  = density\n');     %#ok<ALIGN> 
        elseif  obj.testMode == 2,  fprintf('    Conv. test  = eigenvalues\n');
        else,                       fprintf('    Unexpected convergence test\n');    end
    else
        fprintf('    Conv. test  = %s\n', obj.test);
    end
    
    obj.version; fprintf('\n')

    ref                 =   [ref {'Anderson1965','Johnson1988'}];

    % Bibliography, funding, and footer
    QMol_doc.showBibliography(ref);
    QMol_doc.showFunding;
    if ~obj.isRun,          QMol_doc.showFooter;    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    tol                 =   1e-10
    maxit               =   100
    mix                 =   0.5
    mode                =   'density'
    test                =   'density'
    disp                =   true
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    tolerance               % tol
    maxIterations           % maxit
    mixing                  % mix
    mixingMode              % mode
    convergenceTest         % test
    display                 % disp
end
properties (Transient,Access=private)
    % Eigensolver
    eigs
    % Computation modes
    mixMode
    testMode
    % Currently running an SCF
    isRun               =   false
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % tol ~~~~~~~~~~~~~~
    function set.tolerance(obj,val),            obj.tol     =   val;        end
    function val = get.tolerance(obj),          val         =   obj.tol;    end
    % maxit ~~~~~~~~~~~~
    function set.maxIterations(obj,val),        obj.maxit   =   val;        end
    function val = get.maxIterations(obj),      val         =   obj.maxit;  end
    % mix ~~~~~~~~~~~~~~
    function set.mixing(obj,val),               obj.mix     =   val;        end
    function val = get.mixing(obj),             val         =   obj.mix;    end
    % mode ~~~~~~~~~~~~~
    function set.mixingMode(obj,val),           obj.mode    =   val;        end
    function val = get.mixingMode(obj),         val         =   obj.mode;   end
    % test ~~~~~~~~~~~~~
    function set.convergenceTest(obj,val),      obj.test    =   val;        end
    function val = get.convergenceTest(obj),    val         =   obj.test;   end
    % disp ~~~~~~~~~~~~~
    function set.display(obj,val),              obj.disp    =   val;        end
    function val = get.display(obj),            val         =   obj.disp;   end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset clears all temporary (transient) properties of the object
    
    % Run-time variables
    obj.eigs            =   [];
    obj.mixMode         =   [];
    obj.testMode        =   [];
    obj.isRun           =   false;
    
    % Initialization status
    obj.isInit          =   false;
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.tol),    obj.tol     =   1e-10;      end
    if isempty(obj.maxit),  obj.maxit   =   100;        end
    if isempty(obj.mix),    obj.mix     =   0.5;        end
    if isempty(obj.mode),   obj.mode    =   'density';  end
    if isempty(obj.test),   obj.test    =   'density';  end
    if isempty(obj.disp),   obj.disp    =   true;       end
end
end
methods (Access=public)
function initialize(obj,ES)
%initialize initializes all the atom (center) components. For compatibility
%   with user-defined molecular potential, it supports having an empty atom
%   list
    
    % Set links
    if nargin > 1,  obj.eigs    =   ES;         end

    % Set computation mode
    switch lower(obj.mode)
        case {'density','den','d'},     obj.mixMode     =   1;
        case {'potential','pot','p'},   obj.mixMode     =   2;
        otherwise,                      obj.mixMode     =   1;
            warning('QMolGrid:DFT_SCF_Anderson:mixingMode', ...
                ['Unknown mixing mode ' obj.mode '; using (default) density mixing.'])
    end
    switch lower(obj.test)
        case {'density','den','d'},     obj.testMode    =   1;
        case {'eigenvalues','eigenvalue','eigs','eig'}
                                        obj.testMode    =   2;
        otherwise,                      obj.testMode    =   1;
            warning('QMolGrid:DFT_SCF_Anderson:convergenceTest', ...
                ['Unknown convergence test ' obj.mode '; using (default) density.'])
    end

    % House keeping
    obj.isInit       =   true;
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_DFT_SCF_Anderson';
    PropNames           =  {'tol','maxit','mix','mode','test','disp',...
                            'tolerance','maxIterations','mixing', ...
                            'mixingMode','convergenceTest','display'};
end
end
%% Solve SCF problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function solveSCF(obj,DFT,ES,IC)
%solveSCF solves the SCF problem for the input DFT model using the ES
%   eigensolver. Optionally, the starting density or Kohn-Sham potential
%   objects can be specified as an 4th input
    
    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nargin < 4,  IC  =   [];
    if nargin < 3,  ES  =   [];                 end, end
    
    obj.reset;
    obj.isRun   =   true;
    
    if obj.disp,    QMol_doc.showHeader;        end
    
    % (Re)link everything ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.disp
        QMol_doc.showSection('Build density-functional-theory (DFT) model');
        DFT.initialize(2);
    else
        DFT.initialize(0);
    end

    % Eigen solver
    if isempty(ES), if DFT.disc.isBasis && DFT.disc.nV < 101                %#ok<ALIGN> 
            ES          =   QMol_DFT_eig_basis;                 % full matrix diagonalization
    else,   ES          =   QMol_DFT_eigs;                  end % subset of eigenvalues
    else,   ES.reset;                                       end % input eigensolver


    ES.initialize(DFT);
    obj.initialize(ES);     TM  =   obj.testMode;   MM  =   obj.mixMode;
    
    if obj.disp,            obj.showDocumentation;          end

    % Initialize display ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.disp
        % Header
        QMol_doc.showSection('SCF iterations');
        fprintf('  Iter.  Orbital energies (-eV)                                     Error\n');
        fprintf('  -----  -----------------------------------------------------    --------\n');
    end

    % Allocate necessary objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if MM == 1          % Mixing densities (need 4 different density objects + 1 KS potential)
        RHOi            =   {DFT.disc.DFT_allocateDensity, DFT.disc.DFT_allocateDensity};
        RHOo            =   {DFT.disc.DFT_allocateDensity, DFT.disc.DFT_allocateDensity};
        
        VKSi            =   DFT.disc.DFT_allocatePotential;   VKSi    =   {VKSi,VKSi};
        VKSo            =   VKSi;
    elseif MM == 2      % Mixing potential (need 4 different potential objects)
        VKSi            =   {DFT.disc.DFT_allocatePotential, DFT.disc.DFT_allocatePotential};
        VKSo            =   {DFT.disc.DFT_allocatePotential, DFT.disc.DFT_allocatePotential};

        if TM == 1      % Test convergence from density (need 2 different density objects)
            RHOi        =   DFT.disc.DFT_allocateDensity;   RHOi    =   {RHOi,RHOi};
            RHOo        =   DFT.disc.DFT_allocateDensity;   RHOo    =   {RHOo,RHOo};
        elseif TM == 2  % Test convergence from eigenvalues (need a single density object)
            RHOi        =   DFT.disc.DFT_allocateDensity;   RHOi    =   {RHOi,RHOi};
            RHOo        =   RHOi;
        else,   error('QMolGrid:DFT_SCF_Anderson:testMode', ...
                    'Unexpected testMode (shouldn''t have happened)')
        end
    else
        error('QMolGrid:DFT_SCF_Anderson:mixingMode', ...
            'Unexpected mixMode (shouldn''t have happened)')
    end
    
    % Initial condition ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DFT.setPotentialKernel;
    if isempty(IC)
        % Completely start from scratch (zero density)
        DFT.disc.DFT_allocateDensity(RHOi{1},true);
        
        % Get initial density with proper total charge
        ES.computeEigenstates(DFT.getPotential(RHOi{1},VKSi{1}));
        DFT.getDensity(RHOi{1});

        % Need to update the potential?
        if MM == 2,     DFT.getPotential(RHOi{1},VKSi{1});  end

    elseif ischar(IC) && any(strcmpi(IC,{'self','DFT'}))
        % Start from density computed from the DFT object KSO
        DFT.getDensity(RHOi{1});

        % Make sure we have the proper charge
        N               =   RHOi{1}.getCharge;

        if any(abs(N-DFT.Ntot) < .01*DFT.Ntot)      % charge mismatch > 1%
            % Get initial density with proper total charge
            DFT.getPotential(RHOi{1},VKSi{1});
            VKSi{1}.initialize(DFT.disc);
            ES.computeEigenstates(VKSi{1});
            DFT.getDensity(RHOi{1});
        else
            % Rescale density
            DFT.disc.DFT_mix(RHOi{1},DFT.Ntot./N,RHOi{1});
        end

        % Need to update the potential?
        if MM == 2,     DFT.getPotential(DFT.getDensity(RHOi{1}),VKSi{1});  end

    elseif isa(IC,DFT.disc.DFT_classDensity)
        % Start from input one-body density (charge adjusted)
        N               =   IC.getCharge;

        if any(abs(N-DFT.Ntot) < .01*DFT.Ntot)      % charge mismatch > 1%
            % Get initial density with proper total charge
            DFT.getPotential(IC,VKSi{1});
            VKSi{1}.initialize(DFT.disc);
            ES.computeEigenstates(VKSi{1});
            DFT.getDensity(RHOi{1});
        else
            % Rescale density
            DFT.disc.DFT_mix(RHOi{1},DFT.Ntot./N,IC);
        end

        % Need to update the potential?
        if MM == 2,     DFT.getPotential(RHOi{1},VKSi{1});                  end

    elseif isa(IC,DFT.disc.DFT_classPotential)
        % Start from input Kohn-Sham potential (copy it)
        if DFT.isSpinPol,   DFT.disc.DFT_mix(VKSi{1},[1 1],IC);
        else,               DFT.disc.DFT_mix(VKSi{1},1,IC);                 end

        % Need to get the one-body density?
        if MM == 1,     VKSi{1}.initialize(DFT.disc);                       %#ok<ALIGN> 
                        ES.computeEigenstates(VKSi{1});
                        DFT.getDensity(RHOi{1});                            end
    else
        % Unexpected initial condition
        warning('QMolGrid:DFT_SCF_Anderson:solveSCF', ...
            'Unexpected initial condition for SCF computation; No SCF performed');
        return
    end
    
    % Initial condition is converged? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if MM == 1          % Mixing densities
        VKSi{1}.initialize(DFT.disc);
        Ei              =   ES.computeEigenstates(DFT.getPotential(RHOi{1},VKSi{1}));
        DFT.getDensity(RHOo{1});

        if TM == 2,     VKSo{1}.initialize(DFT.disc);
                        Eo  =   ES.computeEigenstates(DFT.getPotential(RHOo{1},VKSo{1})); end

    elseif MM == 2      % Mixing potentials
        VKSi{1}.initialize(DFT.disc);
        Ei              =   ES.computeEigenstates(VKSi{1});
        DFT.getPotential(DFT.getDensity(RHOi{1}),VKSo{1});

        VKSo{1}.initialize(DFT.disc);
        Eo              =   ES.computeEigenstates(VKSo{1});
        if TM == 1,     DFT.getDensity(RHOo{1});                            end

    else
        % Unexpected mixing mode
        error('QMolGrid:DFT_SCF_Anderson:mixingMode', ...
            'Unexpected mixMode (shouldn''t have happened)')
    end

    % Test for convergence
    DFT.setPotentialKernel;
    if TM == 1          % Change in one-body density
        cv              =   all(DFT.disc.DFT_dist(RHOi{1},RHOo{1}) < obj.tol);

    elseif TM == 2      % Change in sum of eigenvalues
        if DFT.isSpinPol,   cv  =   abs(sum(Ei{1}-Eo{1})) < obj.tol && ...
                                    abs(sum(Ei{2}-Eo{2})) < obj.tol;        %#ok<ALIGN> 
        else,               cv  =   abs(sum(Ei   -Eo   )) < obj.tol;        end

    else                % Unexpected convergence test mode
        error('QMolGrid:DFT_SCF_Anderson:testMode', ...
            'Unexpected convergence test mode (shouldn''t have happened)')
    end
    
    % First iteration (linear mixing) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ~cv
    % Initialization -----------------------
    RHOi    =   {RHOi{2}, RHOi{1}};     RHOo    =   {RHOo{2},RHOo{1}};
    VKSi    =   {VKSi{2}, VKSi{1}};     VKSo    =   {VKSo{2},VKSo{1}};

    if DFT.isSpinPol,   M   =   obj.mix*[1 1];  else,   M = obj.mix; end
    
    % Mixing densities ---------------------
    if MM == 1,         DFT.disc.DFT_mix(RHOi{1},1-M,RHOi{2},M,RHOo{2});    %#ok<ALIGN> 
        % Update results
        DFT.getPotential(RHOi{1},VKSi{1});
        VKSi{1}.initialize(DFT.disc);
        Ei              =   ES.computeEigenstates(VKSi{1});
        DFT.getDensity(RHOo{1});

        if TM == 2,     DFT.getPotential(RHOo{1},VKSo{1});                  %#ok<ALIGN> 
                        VKSo{1}.initialize(DFT.disc);
                        Eo  =   ES.computeEigenstates(VKSo{1});             end

    % Mixing potentials --------------------
    elseif MM == 2,     DFT.disc.DFT_mix(VKSi{1},1-M,VKSi{2},M,VKSo{2});
        % Update results
        VKSi{1}.initialize(DFT.disc);
        Ei              =   ES.computeEigenstates(VKSi{1});
        DFT.getPotential(DFT.getDensity(RHOi{1}),VKSo{1});

        VKSo{1}.initialize(DFT.disc);
        Eo              =   ES.computeEigenstates(VKSo{1});
        if TM == 1,     DFT.getDensity(RHOo{1});                            end

    % Unexpected mixing case -----------
    else,   error('QMolGrid:DFT_SCF_Anderson:mixMode', ...
                'Unexpected convergence test mode (shouldn''t have happened)');     end

    % Compute error ------------------------
    DFT.setPotentialKernel;
    if TM == 1,             err =   DFT.disc.DFT_dist(RHOi{1},RHOo{1});     %#ok<ALIGN> 
    elseif TM == 2
        if DFT.isSpinPol,   err =  [abs(sum(Ei{1}-Eo{1})), abs(sum(Ei{2}-Eo{2}))];
        else,               err =   abs(sum(Ei   -Eo   ));                  end
    else,   error('QMolGrid:DFT_SCF_Anderson:testMode', ...
                'Unexpected convergence test mode (shouldn''t have happened)'); end
    
    if obj.disp,    obj.showIterationResult(0,Ei,err);  end

    % Test for convergence -----------------
    cv                  =   all(err < obj.tol*obj.mix);
end
    % Anderson iterations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ~cv
    for k = 1:obj.maxit
        % Mixing Densities -----------------
        if MM == 1                                                          %#ok<ALIGN> 
            % Mix density and stuff
            b           =   DFT.disc.DFT_SCF_AndersonMixCoeff(RHOi{1},RHOi{2},RHOo{1},RHOo{2});
            DFT.disc.DFT_mix(RHOi{2},(1-obj.mix)*b,RHOi{2},(1-obj.mix)*(1-b),RHOi{1},...
                                        obj.mix *b,RHOo{2},   obj.mix *(1-b),RHOo{1});

            RHOi        =   {RHOi{2}, RHOi{1}};     RHOo    =   {RHOo{2},RHOo{1}};
            VKSi        =   {VKSi{2}, VKSi{1}};     VKSo    =   {VKSo{2},VKSo{1}};

            % Update results
            DFT.getPotential(RHOi{1},VKSi{1});
            VKSi{1}.initialize(DFT.disc);
            Ei          =   ES.computeEigenstates(VKSi{1});
            DFT.getDensity(RHOo{1});
    
            if TM == 2,     DFT.getPotential(RHOo{1},VKSo{1})               %#ok<ALIGN> 
                            VKSo{1}.initialize(DFT.disc);
                            Eo  =   ES.computeEigenstates(VKSo{1});         end

        % Mixing potentials ----------------
        elseif MM == 2
            % Mix potential and stuff
            b           =   DFT.disc.DFT_SCF_AndersonMixCoeff(VKSi{1},VKSi{2},VKSo{1},VKSo{2});
            DFT.disc.DFT_mix(VKSi{2},(1-obj.mix)*b,VKSi{2},(1-obj.mix)*(1-b),VKSi{1},...
                                        obj.mix *b,VKSo{2},   obj.mix *(1-b),VKSo{1});

            RHOi        =   {RHOi{2}, RHOi{1}};     RHOo    =   {RHOo{2},RHOo{1}};
            VKSi        =   {VKSi{2}, VKSi{1}};     VKSo    =   {VKSo{2},VKSo{1}};

            % Update results
            VKSi{1}.initialize(DFT.disc);
            Ei          =   ES.computeEigenstates(VKSi{1});
            DFT.getPotential(DFT.getDensity(RHOi{1}),VKSo{1});
    
            VKSo{1}.initialize(DFT.disc);
            Eo          =   ES.computeEigenstates(VKSo{1});
            if TM == 1,     DFT.getDensity(RHOo{1});                        end
        % Unexpected mixing case -----------
        else,   error('QMolGrid:DFT_SCF_Anderson:mixMode', ...
                    'Unexpected convergence test mode (shouldn''t have happened)');     end

        % Compute error ------------------------
        DFT.setPotentialKernel;
        if TM == 1,             err =   DFT.disc.DFT_dist(RHOi{1},RHOo{1}); %#ok<ALIGN> 
        elseif TM == 2
            if DFT.isSpinPol,   err =  [abs(sum(Ei{1}-Eo{1})), abs(sum(Ei{2}-Eo{2}))];
            else,               err =   abs(sum(Ei   -Eo   ));              end
        else,   error('QMolGrid:DFT_SCF_Anderson:testMode', ...
                    'Unexpected convergence test mode (shouldn''t have happened)'); end
        
        if obj.disp,    obj.showIterationResult(k,Ei,err);                  end

        % Test convergence
        cv              =   all(err < obj.tol*obj.mix);     if cv,  break;  end
    end
end
    % Finalization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.disp,fprintf('  ------------------------------------------------------------------------\n');
        if cv,  fprintf('  Kohn-Sham orbitals converged to tolerance\n\n'); %#ok<ALIGN> 
        else,   fprintf('  WARNING: The Kohn-Sham orbitals did not converge to specified tolerance.\n');
                fprintf('    Try increasing the number of iterations and/or decreasing the mixing\n    parameter.\n\n');  end
        
        % Show orbital and DFT energies
        QMol_doc.showSection('Orbital energies');
        DFT.showEnergy('KSO');      fprintf('\n');
        
        QMol_doc.showSection('DFT-component energies');
        DFT.showEnergy('DFT');      fprintf('\n');
        
        % Footer
        QMol_doc.showFooter;
    end
    
    % Clear member temporary variables
    obj.eigs            =   [];
    obj.isRun           =   false;
end
end
methods (Access=private)
function showIterationResult(~,nb,E,err)
%showIterationResult displays the result of an iteration of the SCF
    
    % Iteration
    fprintf('   %3i  ',nb);

    % Energies and errors
    if iscell(E)
        % Spin up component
        fprintf('%6.2f',convertUnit.au2ev(-E{1}(1)));

        if isscalar(E{1})       %#ok<ALIGN> nothing to do
        elseif numel(E{1}) < 8, fprintf(' |%6.2f',convertUnit.au2ev(-E{1}(2:end)));
        else,                   fprintf(' |%6.2f | (...)',convertUnit.au2ev(-E{1}(2)));
                                fprintf(' |%6.2f',convertUnit.au2ev(-E{1}(end-3:end)));    end
        
        fprintf(['%#' num2str(68-8*min(numel(E{1}),7)) '.3e\n'],err(1));

        % Spin down component
        fprintf('        %6.2f',convertUnit.au2ev(-E{2}(1)));

        if isscalar(E{2})       %#ok<ALIGN> nothing to do
        elseif numel(E{2}) < 8, fprintf(' |%6.2f',convertUnit.au2ev(-E{2}(2:end)));
        else,                   fprintf(' |%6.2f | (...)',convertUnit.au2ev(-E{2}(2)));
                                fprintf(' |%6.2f',convertUnit.au2ev(-E{2}(end-3:end)));    end

        fprintf(['%#' num2str(68-8*min(numel(E{2}),7)) '.3e\n'],err(2));
    else
        % First energy
        fprintf('%6.2f',convertUnit.au2ev(-E(1)));

        % Following energies
        if isscalar(E)          %#ok<ALIGN> nothing to do
        elseif numel(E) < 8,fprintf(' |%6.2f',convertUnit.au2ev(-E(2:end)));
        else,                   fprintf(' |%6.2f | (...)',convertUnit.au2ev(-E(2)));
                                fprintf(' |%6.2f',convertUnit.au2ev(-E(end-3:end)));    end
        % Error
        fprintf(['%#' num2str(68-8*min(numel(E),7)) '.3e\n'],err);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

