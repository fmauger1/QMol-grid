classdef QMol_DFT_Vx_LDA_soft < QMol_suite
%QMol_DFT_Vx_LDA_soft implementation of the local-density-approximation
%   (LDA) Slater-exchange potential for a soft-Coulomb electron-electron
%   interaction potential.

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release
%   01.21.001   07/10/2024  F. Mauger
%       Correct typo in run-time documentation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.001','07/10/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT_Vx_LDA_soft})
function showInfo
    fprintf( '  * QMol_DFT_Vx_LDA_soft:\n');
    fprintf(['      > LDA Slater exchange\n',...
             '      > Soft-Coulomb electron interactions\n']); 
    QMol_DFT_Vx_LDA_soft.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Header
    fprintf('  * Slater-exchange functional           local-density approximation (LDA)\n');
    fprintf('    for a one-dimensional (1D) soft-Coulomb electron-electron interaction\n');
    fprintf('    potential, of the form\n');
    fprintf('        Vee(x) = -%s / sqrt( x^2 + %s^2 )\n',num2str(obj.Z,'%9.3g'),num2str(obj.a,'%9.3g'));
    fprintf('    The LDA Slater exchange is approximated by the scaled energy-per-\n')
    fprintf('    particle, parameterized as [Mauger 2024]\n')
    fprintf('        eps_x(r_s) = -1/2 * (1+alpa*r_s) / (beta*r+2*alpa*m*r_s^2)\n')
    fprintf('                          * log(1+beta*r+gamma*r.^m),\n')
    fprintf('    where r_s = 1/(2*rho) is the 1D Wigner-Seitz radius, and the\n')
    fprintf('    parameters alpha, beta, gamma and m are obtained by fit against the\n')
    fprintf('    exact LDA exchange potential')
    ref                 =   {'Mauger 2024'};

    % SIC (if any)
    if obj.SIC == 0
        fprintf('.\n');
    elseif obj.SIC == 1
        fprintf(', with average-density self-interaction\n    correction (ADSIC) [Legrand 2002].\n');
        ref             =   [ref, {'Legrand2002'}];
    end

    % Version
    obj.version;

end
function mem = getMemoryProfile(~,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components
    mem                 =   0;

    if opt
        QMol_DFT_profiler.showMemoryFootprint('Exchange functional (LDA soft Coulomb)',mem,1);
    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Soft-Coulomb potential parameters
    Z                   =   1
    a                   =   sqrt(2)
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    charge                          % Z
    softeningParameter              % a
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_DFT_Vx_LDA_soft,?QMol_DFT})
    % Linked objects
    DFT                             % DFT.disc always defines the domain
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Run-time variables
    SIC                             % 0 = none, 1 = ADSIC
end
properties(Constant,Access=public)
    type                =   'LDA_X'
end
properties (Constant,Access=private)
    % Potential parametrization
    alpha               =   10.18001817
    beta                =   5.502143989
    gamma               =   14.64700068
    m                   =   2.301803657
    % Threshold density for zero potential
    tol                 =   1e-10
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % Z ~~~~~~~~~~~~~~~~
    function set.charge(obj,val),                   obj.Z       =   val;    end
    function val = get.charge(obj),                 val         =   obj.Z;  end
    % a ~~~~~~~~~~~~~~~~
    function set.softeningParameter(obj,val),       obj.a       =   val;    end
    function val = get.softeningParameter(obj),     val         =   obj.a;  end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.DFT             =   [];
    obj.SIC             =   [];
    
    % Initialization status
    obj.isInit          =   false;
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.V0),     obj.Z       =   1;                              end
    if isempty(obj.s),      obj.a       =   sqrt(2);                        end
end
function initialize(obj,DFT,SIC)
%initialize initializes the object
    
    % Initialization needed?
    obj.reset;
    
    % Set links
    obj.DFT             =   DFT;

    % Identify the flavor of SIC
    if nargin == 2,     SIC = [];       end
    if isempty(SIC),    SIC = 'none';   end

    switch lower(SIC)
        case {'none','off'}
            obj.SIC     =   0;
        case {'adsic','average density','averagedensity','average_density'}
            obj.SIC     =   1;
        otherwise
            warning('QMol:QMol_DFT_Vx_LDA_soft:SIC', ...
                    ['Unknown flavor of self-interaction correction (SIC) ' SIC '\n; SIC ignored for the LDA-exchange potential.'])
            obj.SIC     =   0;
    end
    
    % Miscellaneous
    obj.isInit          =   true;

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_DFT_Vx_LDA_soft';
    PropNames           =   {'Z','a','charge','softeningParameter'};
end
end
%% Get LDA exchange potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function V = getPotential(obj,rho,V,isAdd)
%getPotential returns the discretization of the LDA-exchange potential
%   associated with the member properties and input one-body density (rho).
%   Optionally, provide the potential object(s), where the discretization
%   should be stored, and if the exchange potential should be added to it.
    
    % Initialization
    if nargin < 4,  isAdd   =   false;                                      %#ok<ALIGN> 
    if nargin < 3,  V       =   [];
    if nargin < 2,  rho     =   [];             end, end, end
    
    if isempty(V)   ||   ~isAdd  
        V               =   obj.DFT.disc.DFT_allocatePotential(V);          % If any potential, reset it
    end
    
    if isempty(rho)
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);                % Use DFT one-body density storage
        rho             =   obj.DFT.rho;
    end

    scl                 =   obj.Z/obj.a;

    % Compute the LDA-exchange potential
    if rho.isSpinPol
        % Spin-polarized computation
        V.add(obj.scaledPotential( .25 ./ max(rho.rhoUp*obj.a,obj.tol) ) * scl, ...
              obj.scaledPotential( .25 ./ max(rho.rhoDw*obj.a,obj.tol) ) * scl);
    else
        % Spin-restricted computation
        V.add(obj.scaledPotential( .50 ./ max(rho.rho  *obj.a,obj.tol) ) * scl);
    end

    % Self-interaction correction (SIC == 0 for no correction)
    if obj.SIC == 1
        % ADSIC
        if rho.isSpinPol
            % Spin polarized
            V.add(-obj.scaledPotential(.25./max(rho.rhoUp*obj.a,obj.tol)*obj.DFT.Ntot(1))*scl, ...
                  -obj.scaledPotential(.25./max(rho.rhoDw*obj.a,obj.tol)*obj.DFT.Ntot(2))*scl);
        else
            % Spin restricted
            V.add(-obj.scaledPotential(.25./max(rho.rho  *obj.a,obj.tol)*obj.DFT.Ntot   )*scl);
        end
    end

end
function DV = getPotentialDerivative(obj,dim,rho,DV,isAdd)
%getPotentialDerivative returns the discretization of the derivative of the
%   LDA-exchange potential associated with the input one-body density
%   (rho) and its derivative (D_rho). Optionally, provide the potential
%   object(s), where the discretization should be stored, and if the
%   exchange potential should be added to it.
%   The potential derivative is computed using fast-Fourier transforms
%   (since the potential is local).

    % Initialization
    if nargin < 5,  isAdd   =   false;                                      %#ok<ALIGN> 
    if nargin < 4,  DV      =   [];
    if nargin < 3,  rho     =   [];     end, end, end
    
    if isempty(DV)   ||   ~isAdd  
        DV              =   obj.DFT.disc.DFT_allocatePotentialGradient(DV); % If any potential, reset it
    end

    scl                 =   obj.Z/obj.a;

    % Compute the LDA-exchange potential derivative
    if rho.isSpinPol
        % Spin-polarized computation
        if dim == 1
            DV.add(1,real(ifft(obj.DFT.disc.D .* fft(obj.scaledPotential( .25 ./ max(rho.rhoUp*obj.a,obj.tol) ) * scl))), ...
                     real(ifft(obj.DFT.disc.D .* fft(obj.scaledPotential( .25 ./ max(rho.rhoDw*obj.a,obj.tol) ) * scl))));
        else
            warning('QMol:QMol_Vx_LDA_exp:getPotentialDerivative',['Unexpected dimension (' num2str(dim) ') for potential-derivative computation.'])
        end
    else
        % Spin-restricted computation
        if dim == 1
            DV.add(1,real(ifft(obj.DFT.disc.D .* fft(obj.scaledPotential( .50 ./ max(rho.rho  *obj.a,obj.tol) ) * scl))));
        else
            warning('QMol:QMol_Vx_LDA_exp:getPotentialDerivative',['Unexpected dimension (' num2str(dim) ') for potential-derivative computation.'])
        end
    end

    % Self-interaction correction (SIC == 0 for no correction)
    if obj.SIC == 1
        % ADSIC
        if rho.isSpinPol
            % Spin polarized
            if dim == 1
                DV.add(1,real(ifft(obj.DFT.disc.D .* fft(-obj.scaledPotential(.25./max(rho.rhoUp*obj.a,obj.tol)*obj.DFT.Ntot(1))*scl))), ...
                         real(ifft(obj.DFT.disc.D .* fft(-obj.scaledPotential(.25./max(rho.rhoDw*obj.a,obj.tol)*obj.DFT.Ntot(2))*scl))));
            else
                warning('QMol:QMol_Vx_LDA_exp:getPotentialDerivative',['Unexpected dimension (' num2str(dim) ') for potential-derivative computation.'])
            end
        else
            % Spin restricted
            if dim == 1
                DV.add(1,real(ifft(obj.DFT.disc.D .* fft(-obj.scaledPotential(.25./max(rho.rho  *obj.a,obj.tol)*obj.DFT.Ntot   )*scl))));
            else
                warning('QMol:QMol_Vx_LDA_exp:getPotentialDerivative',['Unexpected dimension (' num2str(dim) ') for potential-derivative computation.'])
            end
        end
    end
end
function setPotentialKernel(~)
%setPotentialKernel no kernel to set
    
end
function E = getEnergy(obj,rho)
%getEnergy returns the exchange energy associated with the member 
%   properties and input one-body density (rho). Empty or missing density 
%   uses the one-body-density of the member DFT object (discouraged).
    
    % Initialization
    if nargin < 2,      rho     =   [];                                     end
    if isempty(rho),    rho     =   obj.DFT.getDensity(obj.DFT.rho);        end

    % Compute exchange energy
    if obj.SIC == 1
        % ADSIC
        if rho.isSpinPol
            % Spin polarized
            r_s         =   .25 ./ max(rho.rhoUp*obj.a,obj.tol);
            E           =       sum(rho.rhoUp.* (obj.scaledEnergyPerParticle(r_s) ...
                                                -obj.scaledEnergyPerParticle(r_s*obj.DFT.Ntot(1))));
            
            r_s         =   .25 ./ max(rho.rhoDw*obj.a,obj.tol);
            E           =   E + sum(rho.rhoDw.* (obj.scaledEnergyPerParticle(r_s) ...
                                                -obj.scaledEnergyPerParticle(r_s*obj.DFT.Ntot(2))));
        else
            % Spin restricted
            r_s         =   .50 ./ max(rho.rho  *obj.a,obj.tol);
            E           =       sum(rho.rho  .* (obj.scaledEnergyPerParticle(   r_s) ...
                                                -obj.scaledEnergyPerParticle(.5*r_s*obj.DFT.Ntot)));
        end
    else
        % No SIC
        if rho.isSpinPol
            % Spin polarized
            r_s         =   .25 ./ max(rho.rhoUp*obj.a,obj.tol);
            E           =       sum(rho.rhoUp.* obj.scaledEnergyPerParticle(r_s));
            
            r_s         =   .25 ./ max(rho.rhoDw*obj.a,obj.tol);
            E           =   E + sum(rho.rhoDw.* obj.scaledEnergyPerParticle(r_s));
        else
            % Spin restricted
            r_s         =   .50 ./ max(rho.rho  *obj.a,obj.tol);
            E           =       sum(rho.rho  .* obj.scaledEnergyPerParticle(r_s));
        end
    end

    E                   =   E * obj.DFT.disc.dx * obj.Z/obj.a;
end
end
%% LDA exchange functional components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=private)
function eps_x = scaledEnergyPerParticle(obj,r_s)
%scaledEnergyPerParticle scaled energy per particle
    
    eps_x               =   -.5 * (1+obj.alpha*r_s) ./ (obj.beta*r_s+2*obj.alpha*obj.m*r_s.^2) .* log(1 + obj.beta*r_s + obj.gamma*r_s.^obj.m);
end
function V_x = scaledPotential(obj,r_s)
%scaledPotential scaled potential

    V_x                 =   -.5 * (2*obj.beta*r_s+obj.alpha*obj.beta*r_s.^2+6*obj.alpha*obj.m*r_s.^2+4*obj.alpha^2*obj.m*r_s.^3) ...
                                ./(obj.beta*r_s+2*obj.alpha*obj.m*r_s.^2).^2 ...
                                .* log(1+obj.beta*r_s+obj.gamma*r_s.^obj.m) ...
                          	+.5 * (1+obj.alpha*r_s)./(obj.beta+2*obj.alpha*obj.m*r_s) ...
                                .* (obj.beta+obj.m*obj.gamma*r_s.^(obj.m-1)) ...
                                ./ (1+obj.beta*r_s+obj.gamma*r_s.^obj.m);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

