classdef QMol_DFT_Vx_XX_conv < QMol_suite
%QMol_DFT_Vx_XX_conv implementation of the exact-exchange (XX) functional.
%   The XX is an implicit functional (functional of the Kohn-Sham orbitals
%   rather than the one-body density) and ignores any self-interaction
%   correction option.

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release
%   01.23.000   07/23/2025  F. Mauger
%       Fix getPotentialDerivative

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.000','07/22/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT_Vx_XX_conv})
function showInfo
    fprintf( '  * QMol_DFT_Vx_XX_conv:\n');
    fprintf(['      > Exact-exchange\n',...
             '      > Explicit-convolution scheme\n']); 
    QMol_DFT_Vx_XX_conv.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Initialization
    ref                 =   {};

    % Header
    fprintf('  * Exact-exchange functional                         explicit convolution\n');
    fprintf('    Interaction pot. = %s (elec.-elec.)\n',func2str(obj.Vee));

    % Version
    obj.version;

end
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components (do not trust V, KSO, KSOup, or KSOdw, they might
    % have been initialized empty; fetch size from domain and parent DFT)
    mem_pot             =   QMol_DFT_profiler.getMemoryFootprint(2*numel(obj.DFT.disc.x)-1,'real');
    mem_KSO             =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.DFT.disc.x),'imag');

    if obj.DFT.isSpinPol,                                                                                   if opt
        % Spin polarized
        QMol_DFT_profiler.showMemoryFootprint('Exact-exchange functional (conv.)', 0,1);
        QMol_DFT_profiler.showMemoryFootprint('interaction potential',mem_pot,2);
        QMol_DFT_profiler.showMemoryFootprint('spin up kernel'       ,numel(obj.DFT.occ{1})*mem_KSO,2);
        QMol_DFT_profiler.showMemoryFootprint('spin down kernel'     ,numel(obj.DFT.occ{2})*mem_KSO,2);     end
        
        mem             =   (numel(obj.DFT.occ{1}) + numel(obj.DFT.occ{2})) * mem_KSO + mem_pot;
    else,                                                                                                   if opt
        QMol_DFT_profiler.showMemoryFootprint('Exact-exchange functional (conv.)', 0,1);
        QMol_DFT_profiler.showMemoryFootprint('interaction potential',mem_pot,2);
        QMol_DFT_profiler.showMemoryFootprint('kernel'               ,numel(obj.DFT.occ)*mem_KSO,2);        end
        
        mem             =   numel(obj.DFT.occ)*mem_KSO + mem_pot;
    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Electron-electron interaction potential
    Vee                 =   @(x) 1./sqrt(x.^2+2)
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    interactionPotential            % Vee
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_DFT_Vx_XX_conv,?QMol_DFT})
    % Linked objects
    DFT                             % DFT.disc always defines the domain
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Run-time variables
    V
end
properties(Constant,Access=public)
    type                =   'XX'
end
properties (Transient,Access=private)
    KSO
    KSOup
    KSOdw
end
properties (Constant,Access=private)
    tol                 =   1e-10   % Threshold occupation for empty Kohn-Sham orbital
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % Vee ~~~~~~~~~~~~~
    function set.interactionPotential(obj,val),     obj.Vee     =   val;        end
    function val = get.interactionPotential(obj),   val         =   obj.Vee;    end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.DFT             =   [];
    obj.V               =   [];

    obj.KSO             =   [];
    obj.KSOup           =   [];
    obj.KSOdw           =   [];
    
    % Initialization status
    obj.isInit          =   false;
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.Vee),    obj.Vee     =   @(x) 1./sqrt(x.^2+2);           end
end
function initialize(obj,DFT,~)  % support passing a SIC option (ignored)
%initialize initializes the object

    % Initialization
    obj.reset;
    
    % Set links
    obj.DFT             =   DFT;
    
    % Extended potential
    obj.V               =   obj.Vee(...
                                (0:2*numel(obj.DFT.disc.x)-2).' * obj.DFT.disc.dx ...
                               -(obj.DFT.disc.x(end)-obj.DFT.disc.x(1)) );
    
    % Miscellaneous
    obj.isInit          =   true;

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_DFT_Vx_XX_conv';
    PropNames           =   {'Vee','interactionPotential'};
end
end
%% Get exact-exchange functional %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function V = getPotential(obj,~,V,isAdd)
%getPotential returns exact-exchange (XX) potential operator associated
%   with the member properties. For consistency with explicit functionals,
%   a one-body-density can be specified as a second argument (but is
%   ignored). Optionally, provide the potential object(s), where the discretization
%   should be stored, and if the exchange potential should be added to it.
    
    % Initialization
    if nargin < 4,  isAdd   =   false;
    if nargin < 3,  V       =   [];                 end, end
    
    if isempty(V)   ||   ~isAdd  
        V               =   obj.DFT.disc.DFT_allocatePotential(V);          % If any potential, reset it
    end

    % Return handle to XX potential functional
    if obj.DFT.isSpinPol,   V.add(@(p,isUp) obj.applyPotential(p,isUp));
    else,                   V.add(@(p)      obj.applyPotential(p));         end

end
function DV = getPotentialDerivative(obj,~,~,DV,isAdd)
%getPotentialDerivative returns the exact-exchange (XX) potential-
%   derivative associated with the member properties. Optionally, provide 
%   the potential object(s), where the discretization should be stored, and
%   if the exchange potential should be added to it.

    % Initialization
    if nargin < 5,  isAdd   =   false;
    if nargin < 4,  DV       =   [];         end, end
    
    if isempty(DV)   ||   ~isAdd  
        DV              =   obj.DFT.disc.DFT_allocatePotential(DV);          % If any potential, reset it
    end

    % Return the handle to the XX potential derivative operator
    if obj.DFT.isSpinPol,   DV.add(1,@(opt,p,isUp) obj.applyPotentialDerivative(opt,p,isUp));   % first input is the dimension
    else,                   DV.add(1,@(opt,p)      obj.applyPotentialDerivative(opt,p));        end
end
function setPotentialKernel(obj)
%setPotentialKernel sets the kernel to be used for the computation of 
%   exact-exchange potentials

    if obj.DFT.disc.isBasis
        % Reconstruct KS orbitals
        disc            =   obj.DFT.disc;
        disc.KSO        =   disc.DFT_reconstructOrbital(obj.DFT.KSO,disc.KSO);

        % Save reconstructed wave function
        if obj.DFT.isSpinPol
            obj.KSOup   =   disc.KSO.KSOup;
            obj.KSOdw   =   disc.KSO.KSOdw;
        else
            obj.KSO     =   disc.KSO.KSO;
        end
    else % Save DFT-object KSO
        if obj.DFT.isSpinPol
            obj.KSOup   =   obj.DFT.KSO.KSOup;
            obj.KSOdw   =   obj.DFT.KSO.KSOdw;
        else
            obj.KSO     =   obj.DFT.KSO.KSO;
        end

    end
end
function Vp = applyPotential(obj,p,isUp)
%applyPotential computes the results of the exact-exchange potential
%   applied to the input one-electron wave function p
    
    % Initialization
    Vp                  =   zeros(numel(obj.DFT.disc.x),1);

    % Compute result of XX potential
    if obj.DFT.isSpinPol
        % Spin polarized
        if isUp
            % Local variable
            N           =   obj.DFT.occ{1};

            % Parse exchange
            for k = 1:numel(N), if N(k) >= obj.tol                          %#ok<ALIGN> 
                Vp      =   Vp - N(k)*obj.KSOup(:,k) .* conv(conj(obj.KSOup(:,k)).*p,obj.V,'same'); 
            end, end
        else
            % Local variable
            N           =   obj.DFT.occ{2};

            % Parse exchange
            for k = 1:numel(N), if N(k) >= obj.tol                          %#ok<ALIGN> 
                Vp      =   Vp - N(k)*obj.KSOdw(:,k) .* conv(conj(obj.KSOdw(:,k)).*p,obj.V,'same'); 
            end, end
        end
    else
        % Spin restricted
        N               =   0.5*obj.DFT.occ;        % include spin-channel resolved scailing

        % Parse exchange
        for k = 1:numel(N), if N(k) >= obj.tol                              %#ok<ALIGN> 
            Vp          =   Vp - N(k)*obj.KSO(:,k) .* conv(conj(obj.KSO(:,k)).*p,obj.V,'same'); 
        end, end

    end

    % Finalize results
    Vp                  =   Vp * obj.DFT.disc.dx;

end
function DVp = applyPotentialDerivative(obj,opt,p,isUp)
%applyPotentialDerivative computes the result of the exact-exchange
%   potential applied to the input one-electron wave function p for the
%   selected type of computation opt
    
    switch lower(opt),                                                      case 'dipacc',                                  if obj.DFT.isSpinPol
        % Dipole acceleration
        DVp             =   2*real(sum(obj.applyPotential(p,isUp) .* ifft(obj.DFT.disc.D.*fft(conj(p))))*obj.DFT.disc.dx);  else
        DVp             =   2*real(sum(obj.applyPotential(p     ) .* ifft(obj.DFT.disc.D.*fft(conj(p))))*obj.DFT.disc.dx);  end
        warning('The exact-exchange dipole acceleration is most likely wrong')
                                                                            otherwise
        % Unexpected case
        error('QMol:QMol_DFT_Vx_XX_conv:applyPotentialDerivative', ...
            ['Unknown application for the exact-exchange potential-derivative operator ' opt])
    end
end
function E = getEnergy(obj,~)
%getEnergy returns the exact-exchange (XX) energy associated with the
%   member properties (and linked DFT Kohn-Sham orbitals).
    
    % Initialization
    if obj.DFT.disc.isBasis
        % Reconstruct KS orbitals
        obj.DFT.disc.KSO=   obj.DFT.disc.DFT_reconstructOrbital(obj.DFT.KSO,obj.DFT.disc.KSO);
        KSO             =   obj.DFT.disc.KSO;                               %#ok<PROPLC> 
    else
        KSO             =   obj.DFT.KSO;                                    %#ok<PROPLC> 
    end
    E                   =   0;

    % Compute XX energy
    if obj.DFT.isSpinPol
        % Up-spin component
        N               =   obj.DFT.occ{1};

        for k = 1:numel(N),     if N(k) >= obj.tol                          %#ok<ALIGN> 
            % "Diagonal" term
            p           =   abs(KSO.KSOup(:,k)).^2;                         %#ok<PROPLC> 
            E           =   E - 0.5*N(k)^2 * sum( p .* conv(p,obj.V,'same') );

            % "Off-diagonal" terms
            for l=k+1:numel(N), if N(k) >= obj.tol                          %#ok<ALIGN> 
                p       =   conj(KSO.KSOup(:,k)) .* KSO.KSOup(:,l);         %#ok<PROPLC> 
                E       =   E - real( N(k)*N(l) * sum( p .* conv(conj(p),obj.V,'same') ) );
            end, end
        end, end

        % Down-spin component
        N               =   obj.DFT.occ{2};

        for k = 1:numel(N),     if N(k) >= obj.tol                          %#ok<ALIGN> 
            % "Diagonal" term
            p           =   abs(KSO.KSOdw(:,k)).^2;                         %#ok<PROPLC> 
            E           =   E - 0.5*N(k)^2 * sum( p .* conv(p,obj.V,'same') );

            % "Off-diagonal" terms
            for l=k+1:numel(N), if N(k) >= obj.tol                          %#ok<ALIGN> 
                p       =   conj(KSO.KSOdw(:,k)) .* KSO.KSOdw(:,l);         %#ok<PROPLC> 
                E       =   E - real( N(k)*N(l) * sum( p .* conv(conj(p),obj.V,'same') ) );
            end, end
        end, end
    else
        % Spin restricted
        N               =   0.5*obj.DFT.occ;

        for k = 1:numel(N),     if N(k) >= obj.tol                          %#ok<ALIGN> 
            % "Diagonal" term
            p           =   abs(KSO.KSO(:,k)).^2;                           %#ok<PROPLC> 
            E           =   E - 0.5*N(k)^2 * sum( p .* conv(p,obj.V,'same') );

            % "Off-diagonal" terms
            for l=k+1:numel(N), if N(k) >= obj.tol                          %#ok<ALIGN> 
                p       =   conj(KSO.KSO(:,k)) .* KSO.KSO(:,l);             %#ok<PROPLC> 
                E       =   E - real( N(k)*N(l) * sum( p .* conv(conj(p),obj.V,'same') ) );
            end, end
        end, end
        
        % Same energy in both up- and down-spin channels
        E               =   2*E;
    end
    
    % Finalize results
    E                   =   E * obj.DFT.disc.dx^2;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

