classdef (Abstract) QMol_DFT < QMol_suite
%QMol_DFT density-functional theory (DFT) level of theory. It defines
%   > (abstract) interface class
%   > common components for DFT-based simulations (spin-polarized models 
%     use cell arrays for the up- and down-spin occupations)
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT})
function showInfo
    fprintf('  * QMol_DFT:\n      > DFT interface\n'); 
    QMol_DFT.version;
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

    fprintf('  * The variational DFT model matches the canonical Hamiltonian formalism\n');
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
    QMol_doc.showSection('External potential');
    
    if isempty(obj.Vext)
        fprintf('  * No external potential specified\n')
    else
        ref             =   [ref, obj.Vext.showDocumentation];
    end

    fprintf('\n')

    % DFT potentials
    QMol_doc.showSection('DFT potential(s)');

    if isempty(obj.Vh)
        fprintf('  * No Hartree potential specified\n')
    else
        ref             =   [ref, obj.Vh.showDocumentation];
    end

    if isempty(obj.Vxc)
        fprintf('  * No exchange-correlation potential specified\n')
    elseif iscell(obj.Vxc)
        for k = 1:numel(obj.Vxc)
            ref         =  [ref, obj.Vxc{k}.showDocumentation];             %#ok<AGROW> 
        end
    else
        ref             =  [ref, obj.Vxc.showDocumentation];
    end

    fprintf('\n')

    % Electronic DOF
    QMol_doc.showSection('System model');

    fprintf('  * Electronic structure                                Kohn-Sham orbitals\n')
    if obj.isSpinPol
        % Spin polarized
        O               =   obj.occ{1};
        N               =   min(8,numel(O));
        fprintf(  '    Up-spin occ. = %4.2f',O(1)); if N > 1, fprintf(' | %4.2f',O(2:N)); end
        while N < numel(O)
            O           =   O(N+1:end);
            N           =   min(8,numel(O));
            fprintf(' |\n                   %4.2f',O(1)); if N > 1, fprintf(' | %4.2f',O(2:N)); end
        end
        
        O               =   obj.occ{2};
        N               =   min(8,numel(O));
        fprintf('\n    Down-spin    = %4.2f',O(1)); if N > 1, fprintf(' | %4.2f',O(2:N)); end
        while N < numel(O)
            O           =   O(N+1:end);
            N           =   min(8,numel(O));
            fprintf(' |\n                   %4.2f',O(1)); if N > 1, fprintf(' | %4.2f',O(2:N)); end
        end
        
        if obj.isInit
            fprintf('\n    Total charge = %5.2f (up) + %5.2f (down) = %5.2f (electrons)\n',obj.Ntot(1),obj.Ntot(2),obj.Ntot(1)+obj.Ntot(2));
        else
            fprintf('\n    DFT object has not beed initialized\n');
        end
        
    else
        % Spin restricted
        O               =   obj.occ;
        N               =   min(8,numel(O));
        fprintf(  '    Occupation   = %4.2f',O(1)); if N > 1, fprintf(' | %4.2f',O(2:N)); end
        while N < numel(O)
            O           =   O(N+1:end);
            N           =   min(8,numel(O));
            fprintf(' |\n                   %4.2f',O(1)); if N > 1, fprintf(' | %4.2f',O(2:N)); end
        end
        
        if obj.isInit
            fprintf('\n    Total charge = %5.2f (electrons)\n',obj.Ntot);
        else
            fprintf('\n    DFT object has not beed initialized\n');
        end
    end

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
%   footprint of the DFT object with all its components initialized and
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

    % Kohn-Sham orbitals
    if obj.isSpinPol,   KSO =   disc.DFT_allocateOrbital([0 0]);            %#ok<PROPLC> 
    else,               KSO =   disc.DFT_allocateOrbital(0);        end     %#ok<PROPLC>
    mem                 =   mem + KSO.getMemoryProfile(opt);                %#ok<PROPLC> 

    % One body density
    rho                 =   disc.DFT_allocateDensity([],false);             %#ok<PROPLC>
    mem                 =   mem + rho.getMemoryProfile(opt);                %#ok<PROPLC> 

    % Kohn-Sham potential
    Vks                 =   disc.DFT_allocatePotential();                   %#ok<PROPLC>
    Vks.disc            =   disc;                                           %#ok<PROPLC>
    mem                 =   mem + Vks.getMemoryProfile(opt);                %#ok<PROPLC>

    % Kohn-Sham gradient
    if ~strcmp(disc.DFT_classPotentialGradient,'N/A') %#ok<PROPLC> 
        DVks            =   disc.DFT_allocatePotentialGradient();           %#ok<PROPLC>
        DVks.disc       =   disc;                                           %#ok<PROPLC>
        mem             =   mem + DVks.getMemoryProfile(opt);               %#ok<PROPLC>
    end

    % External functional
    if ~isempty(obj.Vext)
        if ~obj.Vext.isInit,    obj.Vext.DFT    =   obj;                    end
        mem             =   mem + obj.Vext.getMemoryProfile(opt);
        if ~obj.Vext.isInit,    obj.Vext.reset();                           end
    end

    % Hartree functional
    if ~isempty(obj.Vh)
        if ~obj.Vh.isInit,      obj.Vh.DFT    =   obj;                      end
        mem             =   mem + obj.Vh.getMemoryProfile(opt);
        if ~obj.Vh.isInit,      obj.Vh.reset();                             end
    end

    % Exchange-correlation functionals
    if ~isempty(obj.Vh),        if iscell(obj.Vxc), for k = 1:numel(obj.Vxc) %#ok<ALIGN> 
        % Multiple exchange-correlation functionals
        if ~obj.Vxc{k}.isInit,  obj.Vxc{k}.DFT    =   obj;                  end
        mem             =   mem + obj.Vxc{k}.getMemoryProfile(opt);
        if ~obj.Vxc{k}.isInit,  obj.Vxc{k}.reset();                         end, end, else
        
        % Single exchange-correlation functional
        if ~obj.Vxc.isInit,     obj.Vxc.DFT    =   obj;                     end
        mem             =   mem + obj.Vxc.getMemoryProfile(opt);
        if ~obj.Vxc.isInit,     obj.Vxc.reset();                            end

    end, end

    % Clean up
    obj.disc            =   inDisc;
    clear disc KSO rho Vks
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Abstract,Constant,Access=public)
    isSpinPol                       % To be defined by the implementation
    dim
end
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Discretization
    disc
    % Electronic structure
    occ
    KSO
    % DFT potentials
    Vext
    Vh
    Vxc
    % Miscellaneous
    SIC                 =   'none'
end
properties (Hidden,Transient,GetAccess=public,SetAccess=?QMol_DFT)
    % Total charge
    Ntot
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    discretization                  % disc
    occupation                      % occ
    orbital                         % KSO
    externalPotential               % Vext
    HartreePotential                % Vh
    exchangeCorrelationPotential    % Vxc
    selfInteractionCorrection       % SIC
    totalCharge                     % Ntot
end
properties (Transient,Hidden,Access=?QMol_suite)
    % One-body density (holding memory)
    rho
    % Kohn-Sham potential (holding memory)
    Vks
    DVks
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % disc ~~~~~~~~~~~~~
    function set.discretization(obj,val),                   obj.disc    =   val;        end
    function val = get.discretization(obj),                 val         =   obj.disc;   end
    % occ ~~~~~~~~~~~~~~
    function set.occupation(obj,val),                       obj.occ     =   val;        end
    function val = get.occupation(obj),                     val         =   obj.occ;    end
    % KSO ~~~~~~~~~~~~~~
    function set.orbital(obj,val),                          obj.KSO     =   val;        end
    function val = get.orbital(obj),                        val         =   obj.KSO;    end
    % Vext ~~~~~~~~~~~~~
    function set.externalPotential(obj,val),                obj.Vext    =   val;        end
    function val = get.externalPotential(obj),              val         =   obj.Vext;   end
    % Vh ~~~~~~~~~~~~~~~
    function set.HartreePotential(obj,val),                 obj.Vh      =   val;        end
    function val = get.HartreePotential(obj),               val         =   obj.Vh;     end
    % Vxc ~~~~~~~~~~~~~~
    function set.exchangeCorrelationPotential(obj,val),     obj.Vxc     =   val;        end
    function val = get.exchangeCorrelationPotential(obj),   val         =   obj.Vxc;    end
    % SIC ~~~~~~~~~~~~~~
    function set.selfInteractionCorrection(obj,val),        obj.SIC     =   val;        end
    function val = get.selfInteractionCorrection(obj),      val         =   obj.SIC;    end
    % Ntot ~~~~~~~~~~~~~
    function val = get.totalCharge(obj),                    val         =   obj.Ntot;   end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object
    
    % Run-time variables
    obj.Ntot            =   [];
    if ~isempty(obj.rho),   delete(obj.rho);    obj.rho         =   [];     end
    if ~isempty(obj.Vks),   delete(obj.Vks);    obj.Vks         =   [];     end
    if ~isempty(obj.DVks),  delete(obj.DVks);   obj.DVks        =   [];     end

    % Also reset all components 
    % (cross-dependencies might be broken);
    if ~isempty(obj.disc),  obj.disc.reset;     end
    if ~isempty(obj.KSO),   obj.KSO.reset;      end
    if ~isempty(obj.Vext),  obj.Vext.reset;     end
    if ~isempty(obj.Vh),    obj.Vh.reset;       end
    if ~isempty(obj.Vxc)
        if iscell(obj.Vxc), for k=1:numel(obj.Vxc), obj.Vxc{k}.reset; end
        else,                                       obj.Vxc.reset;    end
    end
    
    % Initialization status
    obj.isInit          =   false;
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.SIC),    obj.SIC     =   'none'; end
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
        warning('QMol:DFT:noDiscretization','Missing the discretization component; initialization skipped on it.')
        buildStatus(false);
    else
        % Check if discretization has changed
        if ~obj.disc.isInit,        obj.reset;      end
        
        % Initialize the discretization component
        obj.disc.initialize(obj);
        buildStatus(obj.disc.isInit);
    end

    % Electronic structure
    if obj.isSpinPol
        % Spin polarized
        if iscell(obj.occ)
            % Expected format
            obj.occ     =  {obj.occ{1}(:).', obj.occ{2}(:).'};
        else
            % Implicit identical up- and down- spin occupations
            obj.occ     =  {obj.occ(:).', obj.occ(:).'};
        end

        obj.Ntot    =  [sum(obj.occ{1}), sum(obj.occ{2})];
    else
        % Spin restricted
        if iscell(obj.occ)
            % Correct occupation definition
            obj.occ     =   obj.occ{1}(:).';
        else
            % Expected format
            obj.occ     =   obj.occ(:).';
        end

        obj.Ntot        =   sum(obj.occ);
    end

    % Kohn-Sham orbitals
    buildSection('Kohn-Sham orbitals');
    if isempty(obj.disc)   ||   ~obj.disc.isInit
        warning('QMol:DFT:KSO',['Unable to initialized the Kohn-Sham orbitals (missing or uninitialized discretization component).\n'...
                                'No initialization performed on the Kohn-Sham orbitals.'])
        buildStatus(false);
    else
        % Number of orbitals
        if obj.isSpinPol,   N   =   [numel(obj.occ{1}) numel(obj.occ{2})];
        else,               N   =   numel(obj.occ);                         end

        % Allocate orbitals
        obj.KSO         =   obj.disc.DFT_allocateOrbital(N,obj.KSO);        % Allocation initializes orbitals

        buildStatus(obj.KSO.isInit);
    end

    % External potential
    buildSection('External potential');
    if isempty(obj.Vext)
        warning('QMol:DFT:noExternalPotential','Missing the external potential component; initialization skipped on it.')
        buildStatus(false)
    else
        obj.Vext.initialize(obj);
        buildStatus(obj.Vext.isInit);
    end
    
    % SIC options
    if iscell(obj.SIC)
        if numel(obj.SIC) == 1                                              %#ok<ISCL>
            % Same flavor of SIC for all DFT potentials
            SIC_Vh          =   obj.SIC;
            SIC_Vxc         =   obj.SIC;
        elseif numel(obj.SIC) == 2
            % Different SIC for Hartree and exchange-correlation potentials
            SIC_Vh          =   obj.SIC{1};
            SIC_Vxc         =   obj.SIC{2};
        elseif numel(obj.SIC) == numel(obj.Vxc) + 1
            % Different SIC for Hartree and exchange-correlation potentials
            SIC_Vh          =   obj.SIC{1};
            SIC_Vxc         =   obj.SIC(2:end);
        else
            % unexpected format, ignore SIC
            warning('QMol:DFT:SIC','Unexpected SIC specification; SIC imposed on none of the DFT functionals.')
            SIC_Vh          =   'none';
            SIC_Vxc         =   'none';
        end
    else
        % Same flavor of SIC for all DFT potentials
        SIC_Vh          =   obj.SIC;
        SIC_Vxc         =   obj.SIC;
    end

    if iscell(obj.Vxc)   &&   ~iscell(SIC_Vxc)
        SIC_Vxc         =   repmat({SIC_Vxc},1,numel(obj.Vxc));
    end

    % DFT potentials
    if isempty(obj.Vh)
        warning('QMol:DFT:noHartreePotential','Missing the Hartree potential component; initialization skipped on it.')
        buildSection('(no) Hartree potential'); buildStatus(true);
    else
        buildSection('Hartree potential');
        obj.Vh.initialize(obj,SIC_Vh);
        buildStatus(obj.Vh.isInit);
    end

    if isempty(obj.Vxc)
        buildSection('No exchange-correlation potential'); buildStatus(true);
    else
        buildSection('Exchange-correlation potential');

        R               =   true;
        if iscell(obj.Vxc)
            for k = 1:numel(obj.Vxc)
                obj.Vxc{k}.initialize(obj,SIC_Vxc{k});
                R       =   R   &&   obj.Vxc{k}.isInit;
            end
        else
            obj.Vxc.initialize(obj,SIC_Vxc);
            R           =   obj.Vxc.isInit;
        end

        buildStatus(R);
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
    
    ClassName           =   'QMol_DFT';
    PropNames           =  {'disc','occ','KSO','Vext','Vh','Vxc','SIC',     ...
                            'discretization','occupation','orbital',        ...
                            'externalPotential','HartreePotential',         ...
                            'exchangeCorrelationPotential',                 ...
                            'selfInteractionCorrection'};
end
end
%% Interface methods (to be overloaded) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Abstract,Access=?QMol_DFT)
    % Display
    showDoc                 % Show implementation-specific documentation
                            % Show energy
end
%% Compute DFT components (get) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function rho = getDensity(obj,rho)
%getDensity computes and returns the one-body density associated with the
%   occupation and KSO member properties. Optionally provide the one-body
%   density object as a second input for in-place computation. For
%   spin-polarized DFT systems, the up- and down-spin components are
%   combined in a cell array. In all cases, the DFT object is assumed to
%   have been (successfully) initialized.
    
    % Initialization
    if nargin == 1,         rho     =   [];                                 end
    if isempty(rho),        rho     =   obj.disc.DFT_allocateDensity;       end
    
    % Compute the one-body density (getDensity does the initialization)
    rho.getDensity(obj.occ,obj.KSO);

end
function setPotentialKernel(obj)
%setPotentialKernel parse through the exchange-correlation potential
%   components and (re)set their kernels. This step is typically required
%   for implicit functionals
    
    if ~isempty(obj.Vxc)                                                    % exchange-correlation potential
        if iscell(obj.Vxc), for k = 1:numel(obj.Vxc)                        %#ok<ALIGN> 
                            obj.Vxc{k}.setPotentialKernel;                  end
        else,               obj.Vxc.setPotentialKernel;                     end
    end

end
function V = getPotential(obj,rho,V)
%getPotential computes and returns the Kohn-Sham (KS) potential -- sum of
%   external, Hartree, and exchange-correlation components -- associated
%   with the DFT model. Optionally, as a second argument, provide the
%   one-body-density object to be used for the computation of the KS
%   potential. Empty or unspecified density object uses the one-body
%   density associated with the KS orbitals of the DFT object. Optionally,
%   as a third argument, provide the potential object in which the results
%   should be stored. If empty or unspecified, a new potential object is
%   allocated. In all cases, the DFT object is assumed to have been
%   (successfully) initialized.
    
    % Initialization
    if nargin < 3,  V   =   [];
    if nargin < 2,  rho =   [];     end, end
    
    if isempty(rho)
        obj.rho         =   obj.getDensity(obj.rho);                        % Store in one-body-density property
        rho             =   obj.rho;
    end

    % Compute KS potential
    V                   =   obj.Vext.getPotential(V,false);                 % external potential
    if ~isempty(obj.Vh),    obj.Vh.getPotential(rho,V,true);        end     % Hartree potential
    if ~isempty(obj.Vxc)                                                    % exchange-correlation potential
        if iscell(obj.Vxc), for k = 1:numel(obj.Vxc)                        %#ok<ALIGN> 
                            obj.Vxc{k}.getPotential(rho,V,true);    end
        else,               obj.Vxc.getPotential(rho,V,true);       end
    end

    % Initialize the potential
    V.initialize(obj.disc);

end
function DV = getPotentialGradient(obj,rho,DV,dim)
%getPotentialGradient computes and returns the Kohn-Sham (KS) potential 
%   gradient -- sum of external, Hartree, and exchange-correlation
%   components -- associated with the DFT model. Optionally, as a second
%   argument, provide the one-body-density object to be used for the
%   computation of the KS potential gradient. Empty or unspecified density
%   object uses the one-body density associated with the KS orbitals of the
%   DFT object. Optionally, as a third argument, provide the potential 
%   gradient object in which the results should be stored. If empty or
%   unspecified, a new potential gradient object is allocated. Optionally,
%   as a fourth argument, specify the dimensions along which the gradient
%   should be computed. If empty or unspecified, the gradient is computed
%   in all possible directions. In all cases, the DFT object is assumed to
%   have been (successfully) initialized.
    
    % Initialization
    if nargin < 4,  dim =   [];                                             %#ok<ALIGN> 
    if nargin < 3,  DV  =   [];
    if nargin < 2,  rho =   [];     end, end, end
    
    if isempty(dim),        dim     =   1:obj.dim;                          end
    if isempty(rho)
        obj.rho         =   obj.getDensity(obj.rho);                        % Store in one-body-density property
        rho             =   obj.rho;
    end

    DV                  =   obj.disc.DFT_allocatePotentialGradient(DV);

    % Compute KS potential
    for d = dim
        obj.Vext.getPotentialDerivative(d,DV,true);                                         % external potential
        if ~isempty(obj.Vh),    obj.Vh.getPotentialDerivative(d,rho,DV,true);       end     % Hartree potential
        if ~isempty(obj.Vxc), if iscell(obj.Vxc), for k = 1:numel(obj.Vxc)      %#ok<ALIGN> % exchange-correlation potential
                                obj.Vxc{k}.getPotentialDerivative(d,rho,DV,true);   end,    else
                                obj.Vxc.getPotentialDerivative(d,rho,DV,true);      end,    end
    end

    % Initialize the potential
    DV.initialize(obj.disc);

end
function [E1,E2,E3,E4,E5] = getEnergy(obj,I1,I2)
%getEnergy computes the DFT or Kohn-Sham-orbital (KSO) energy components.
%
%   [E,err] = getEnergy('orbital') computes the energies E and errors err
%       for the member-property KSOs. Optionally, the one-body-density or
%       Kohn-Sham potential object to be used for the computation of the
%       energies and error can be included (as a 1st or 2nd argument).
%
%   [Etot,Ekin,Eext,Eh,Exc] = getEnergy('DFT') computes the total, kinetic,
%       external, Hartree, and exchange-correlation energy components,
%       respectively, for the member-property KSOs. 

    % Initialization 
    if nargin == 1,     opt     =   'dft';  IC  =   [];                     %#ok<ALIGN> 
    elseif nargin == 2
        if ischar(I1),  opt     =   I1;     IC  =   [];
        else,           opt     =   'dft';  IC  =   I1;                     end
    else
        if ischar(I1),  opt     =   I1;     IC  =   I2;
        else,           opt     =   I2;     IC  =   I1;                     end, end
    
    obj.setPotentialKernel;

    % Compute energy components
    switch lower(opt)
        % DFT energy -----------------------
        case 'dft'
            % Initialization
            obj.getDensity(obj.rho);                        % Compute one-body density

            % Get energy components
            E2          =   obj.disc.DFT_energyKinetic(obj.occ,obj.KSO);

            if ~isempty(obj.Vext),  E3  =   obj.Vext.getEnergy(obj.rho);    %#ok<ALIGN> 
            elseif obj.isSpinPol,   E3  =   [0 0];
            else,                   E3  =   0;                              end

            if ~isempty(obj.Vh),    E4  =   obj.Vh.getEnergy(obj.rho);
            else,                   E4  =   0;                              end

            E5  =   0;
            if ~isempty(obj.Vxc)
                if iscell(obj.Vxc)                                          %#ok<ALIGN> 
                    for k=1:numel(obj.Vxc), E5  =   E5 + obj.Vxc{k}.getEnergy(obj.rho); end
                else,                       E5  =   obj.Vxc.getEnergy(obj.rho);         end
            end
            
            % Total energy
            E1 = sum(E2) + sum(E3) + sum(E4) + E5;
        % Orbital energies -----------------
        case {'kso','orbital','orb'}
            % Initialization
            if isempty(IC)
                V       =   obj.getPotential(obj.getDensity(obj.rho),obj.Vks);
            elseif isa(IC,obj.disc.DFT_classDensity)
                V       =   obj.getPotential(IC,obj.Vks);
            elseif isa(IC,obj.disc.DFT_classPotential)
                V       =   IC;
            else
                V       =   obj.getPotential(obj.getDensity(obj.rho),obj.Vks);
                warning('QMol:DFT:getEnergy', ...
                    'Expecting a potential or density object; use DFT-generated potential instead.')
            end

            if nargout == 1
                E1      =   obj.disc.DFT_energyOrbital(V,obj.KSO);
            else
                [E1,E2] =   obj.disc.DFT_energyOrbital(V,obj.KSO);
            end
        % Unexpected case ------------------
        otherwise
            E1 = NaN;   E2 = NaN;   E3 = NaN;   E4 = NaN;   E5 = NaN;
            warning('QMol:DFT:getEnergy', ...
                ['Unknown energy-computation case ' opt '; no energy computed.'])

    end
end
function showEnergy(obj,opt)
%showEnergy displays the DFT or Kohn-Sham-orbital (KSO) energy components.
%
%   showEnergy('orbital') displays the energies and errors for the member-
%       property KSOs. 
%
%   showEnergy('DFT') displays the total, kinetic, external, Hartree, and 
%       exchange-correlation energy components for the member-property 
%       KSOs. 
    
    % Defaul option
    if nargin < 2,      opt     =   [];         end
    if isempty(opt),    opt     =   'kso';      end

    % Show energy
    switch lower(opt)
        % DFT energy -----------------------
        case 'dft'
            % Initialization
            [Etot,Ekin,Eext,Eh,Exc]     =   obj.getEnergy('DFT');

            if obj.isSpinPol
                fprintf('  Component      Energy (a.u.)      Energy (eV)      [ spin up/down (eV) ]\n')
                fprintf('  -----------    -------------     -------------     ---------------------\n');
            else
                fprintf('  Component      Energy (a.u.)      Energy (eV)\n')
                fprintf('  -----------    -------------     -------------\n');
            end
            
            % Show energy component
            if obj.isSpinPol
                fprintf('  Kinetic        %#10.3f         %#10.3f       [%#8.2f/ %#8.2f ]\n',...
                        [sum(Ekin), convertUnit.au2ev(sum(Ekin)),convertUnit.au2ev(Ekin)]);
                fprintf('  External       %#10.3f         %#10.3f       [%#8.2f/ %#8.2f ]\n',...
                        [sum(Eext), convertUnit.au2ev(sum(Eext)),convertUnit.au2ev(Eext)]);
            else
                fprintf('  Kinetic        %#10.3f         %#10.3f\n',Ekin,convertUnit.au2ev(Ekin));
                fprintf('  External       %#10.3f         %#10.3f\n',Eext,convertUnit.au2ev(Eext));
            end
            fprintf('  Hartree        %#10.3f         %#10.3f\n',Eh  ,convertUnit.au2ev(Eh  ));
            fprintf('  Exch.-corr.    %#10.3f         %#10.3f\n',Exc ,convertUnit.au2ev(Exc ));
            fprintf('  -----------    -------------     -------------\n');
            fprintf('  Total          %#10.3f         %#10.3f\n',Etot,convertUnit.au2ev(Etot));

            % Finalization
            fprintf('  ----------------------------------------------\n');

        % Orbital energies -----------------
        case {'kso','orbital','orb'}
            % Initialization
            [E,err]     =   obj.getEnergy('KSO');

            fprintf('  Orbital      Occ. (elec.)         Energy (-eV)               Error(a.u.)\n')
            fprintf('  -------      ------------         ------------               -----------\n');

            % Display results
            if obj.isSpinPol
                fprintf('    %3i        %8.2f             %9.3f                  %#10.3e\n', ...
                    [1:numel(obj.occ{1}); obj.occ{1}; convertUnit.au2ev(-E{1}.'); err{1}.']);
                fprintf('  -------      ------------         ------------               -----------\n');
                fprintf('    %3i        %8.2f             %9.3f                  %#10.3e\n', ...
                    [1:numel(obj.occ{2}); obj.occ{2}; convertUnit.au2ev(-E{2}.'); err{2}.']);
            else
                fprintf('    %3i        %8.2f             %9.3f                  %#10.3e\n', ...
                    [1:numel(obj.occ); obj.occ; convertUnit.au2ev(-E.'); err.']);
            end

            % Finalization
            fprintf('  ------------------------------------------------------------------------\n');

        % Unexpected case ------------------
        otherwise
            warning('QMol:DFT:showEnergy', ...
                ['Unknown energy-display case ' opt '; no energy shown.'])

    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

