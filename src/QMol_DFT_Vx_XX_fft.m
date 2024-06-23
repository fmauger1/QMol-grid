classdef QMol_DFT_Vx_XX_fft < QMol_suite
%QMol_DFT_Vx_XX_fft implementation of the exact-exchange (XX) functional.
%   The XX is an implicit functional (functional of the Kohn-Sham orbitals
%   rather than the one-body density) and ignores any self-interaction
%   correction option. The convolutions involved in the computation of the
%   exact exchange are computed using fast-Fourier transforms over an
%   extended domain. It is almost always faster than the convolution-based
%   schenme but requires tuning of the extended domain.

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
    fprintf( '  * QMol_DFT_Vx_XX_fft:\n');
    fprintf(['      > Exact-exchange\n',...
             '      > Fast-Fourier-transform convolution scheme\n']); 
    QMol_DFT_Vx_XX_fft.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Initialization
    ref                 =   {};

    % Header
    fprintf('  * Exact-exchange functional           fast-Fourier transform convolution\n');

    % Electron-interaction potential
    fprintf('    Interaction pot. = %s (elec.-elec.)\n',func2str(obj.Vee));
    fprintf('    Ghost domain     = %s (%s)',num2str(obj.gL),obj.gU);
    if obj.aFL, fprintf(', plus\n'); else, fprintf(', including\n'); end
    fprintf('    Falloff          = %s (%s), %s shape\n',num2str(obj.fL),obj.fU,obj.fS);

    if obj.isInit
        nb              =   factor(obj.N);
        fprintf('    Ext. domain size = %u ',obj.N);
        if isscalar(nb),    fprintf('(prime) points\n');
        else,               fprintf('(%u',nb(1)); fprintf(' x %u',nb(2:end)); fprintf(') points\n'); end
    end

    % Version
    obj.version;

end
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

     % Evaluate extended potential length
    dx                  =   obj.DFT.disc.x(2)-obj.DFT.disc.x(1);
    Nb                  =   numel(obj.DFT.disc.x);

    switch lower(obj.gU)                                        % Ghost domain
        case {'point','points','grid','vertex','vertices'},     Nb      =   Nb + obj.gL;
        case {'a.u.','au','length','distance'},                 Nb      =   Nb + ceil(obj.gL/dx);
        otherwise,  warning('QMol:QMol_DFT_Vh_fft:gU', ['Unknown ghostUnit ' obj.gU '; Ghost domain ignored.'])
    end
    
    if obj.addFalloffLength, switch lower(obj.fU)      %#ok<ALIGN> % Add falloff to ghost domain
        case {'point','points','grid','vertex','vertices'},     Nb      =   Nb + 2*obj.fL;
        case {'a.u.','au','length','distance'},                 Nb      =   Nb + 2*ceil(obj.fL/dx);
        otherwise,  warning('QMol:QMol_DFT_Vh_fft:fU', ['Unknown falloffUnit ' obj.fU '; Ghost domain ignored.'])
    end, end

    % Display components (do not trust V, KSO, KSOup, or KSOdw, they might
    % have been initialized empty; fetch size from domain and parent DFT)
    mem_pot             =   QMol_DFT_profiler.getMemoryFootprint(Nb,'imag');
    mem_KSO             =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.DFT.disc.x),'imag');

    if obj.DFT.isSpinPol,                                                                                   if opt
        % Spin polarized
        QMol_DFT_profiler.showMemoryFootprint('Exact-exchange functional (fft)', 0,1);
        QMol_DFT_profiler.showMemoryFootprint('interaction potential',mem_pot,2);
        QMol_DFT_profiler.showMemoryFootprint('spin up kernel'       ,numel(obj.DFT.occ{1})*mem_KSO,2);
        QMol_DFT_profiler.showMemoryFootprint('spin down kernel'     ,numel(obj.DFT.occ{2})*mem_KSO,2);     end
        
        mem             =   (numel(obj.DFT.occ{1}) + numel(obj.DFT.occ{2})) * mem_KSO + mem_pot;
    else,                                                                                                   if opt
        QMol_DFT_profiler.showMemoryFootprint('Exact-exchange functional (fft)', 0,1);
        QMol_DFT_profiler.showMemoryFootprint('interaction potential',mem_pot,2);
        QMol_DFT_profiler.showMemoryFootprint('kernel'               ,numel(obj.DFT.occ)*mem_KSO,2);        end
        
        mem             =   numel(obj.DFT.occ)*mem_KSO + mem_pot;
    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (GetAccess=public,SetAccess=?QMol_suite)
    % Electron-electron interaction potential
    Vee                 =   @(x) 1./sqrt(x.^2+2)
    % Ghost domain
    gU                  =   'points'
    gL                  =   0
    % Falloff
    fU                  =   'points'
    fS                  =   'cos^2'
    fL                  =   0
    aFL                 =   true
end
properties (Dependent,Hidden,GetAccess=public,SetAccess=?QMol_suite)
    interactionPotential            % Vee
    ghostUnit                       % gU
    ghostLength                     % gL
    falloffUnit                     % fU
    falloffShape                    % fS
    falloffLength                   % fL
    addFalloffLength                % aFL
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Linked objects
    DFT                             % DFT.disc.dom always defines the domain
    % Run-time variables
    V
    N
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
    % Vee ~~~~~~~~~~~~~~
    function set.interactionPotential(obj,val),     obj.Vee     =   val;        end
    function val = get.interactionPotential(obj),   val         =   obj.Vee;    end

    % gU ~~~~~~~~~~~~~~~
    function set.ghostUnit(obj,val),                obj.gU      =   val;        end
    function val = get.ghostUnit(obj),              val         =   obj.gU;     end

    % gL ~~~~~~~~~~~~~~~
    function set.ghostLength(obj,val),              obj.gL      =   val;        end
    function val = get.ghostLength(obj),            val         =   obj.gL;     end

    % fU ~~~~~~~~~~~~~~~
    function set.falloffUnit(obj,val),              obj.fU      =   val;        end
    function val = get.falloffUnit(obj),            val         =   obj.fU;     end

    % fS ~~~~~~~~~~~~~~~
    function set.falloffShape(obj,val),             obj.fS      =   val;        end
    function val = get.falloffShape(obj),           val         =   obj.fS;     end

    % fL ~~~~~~~~~~~~~~~
    function set.falloffLength(obj,val),            obj.fL      =   val;        end
    function val = get.falloffLength(obj),          val         =   obj.fL;     end

    % aFL ~~~~~~~~~~~~~~
    function set.addFalloffLength(obj,val),         obj.aFL     =   val;        end
    function val = get.addFalloffLength(obj),       val         =   obj.aFL;    end

end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.DFT             =   [];
    obj.V               =   [];
    obj.N               =   [];

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
    if isempty(obj.gU),     obj.gU      =   'points';                       end
    if isempty(obj.gL),     obj.gL      =   0;                              end
    if isempty(obj.fU),     obj.fU      =   'points';                       end
    if isempty(obj.fS),     obj.fS      =   'cos^2';                        end
    if isempty(obj.fL),     obj.fL      =   0;                              end
    if isempty(obj.aFL),    obj.aFL     =   true;                           end
end
function initialize(obj,DFT,~)  % support passing a SIC option (ignored)
%initialize initializes the object
    
    % Initialization needed?
    if obj.isInit,   return; end
    
    % Set links
    obj.DFT             =   DFT;
    
    % Extended potential
    obj.N               =   numel(obj.DFT.disc.x);

    switch lower(obj.gU)                                        % Ghost domain
        case {'point','points','grid','vertex','vertices'},     obj.N   =   obj.N + obj.gL;
        case {'a.u.','au','length','distance'},                 obj.N   =   obj.N + ceil(obj.gL/obj.DFT.disc.dx);
        otherwise,  warning('QMol:QMol_DFT_Vh_fft:gU', ['Unknown ghostUnit ' obj.gU '; Ghost domain ignored.'])
    end
    
    if obj.addFalloffLength, switch lower(obj.fU)      %#ok<ALIGN> % Add falloff to ghost domain
        case {'point','points','grid','vertex','vertices'},     obj.N   =   obj.N + 2*obj.fL;
        case {'a.u.','au','length','distance'},                 obj.N   =   obj.N + 2*ceil(obj.fL/obj.DFT.disc.dx);
        otherwise,  warning('QMol:QMol_DFT_Vh_fft:fU', ['Unknown falloffUnit ' obj.fU '; Ghost domain ignored.'])
    end, end
    
    x                   =   (0:obj.N-1) * obj.DFT.disc.dx;
    ind                 =   x > .5*x(end);
    x(ind)              =   x(ind) - x(end) - obj.DFT.disc.dx;  % Need to start exactly at the bottom of Vee potential
    
    obj.V               =   obj.Vee(x);
    
    % Falloff
    switch lower(obj.fU)
        case {'point','points','grid','vertex','vertices'},     L       =   obj.fL * obj.DFT.disc.dx;
        case {'a.u.','au','length','distance'},                 L       =   obj.fL;
        otherwise,                                              L       =   0;
            warning('QMol:QMol_DFT_Vh_fft:fU', ['Unknown falloffUnit ' obj.fU '; Ghost domain ignored.'])
    end
    
    ind_1               =   x < min(x) + L;
    ind_2               =   x > max(x) - L;
    
    switch lower(obj.fS)
        case {'none','off'}
            % Nothing to do
        case {'cos2','cosine2','sin2','sine2','cos^2','cosine^2','sin^2','sine^2'}
            obj.V(ind_1)=   obj.V(ind_1) .* sin((x(ind_1)-min(x))*.5*pi/L).^2;
            obj.V(ind_2)=   obj.V(ind_2) .* sin((x(ind_2)-max(x))*.5*pi/L).^2;
        case {'cos4','cosine4','sin4','sine4','cos^4','cosine^4','sin^4','sine^4'}
            obj.V(ind_1)=   obj.V(ind_1) .* sin((x(ind_1)-min(x))*.5*pi/L).^4;
            obj.V(ind_2)=   obj.V(ind_2) .* sin((x(ind_2)-max(x))*.5*pi/L).^4;
        case 'gaussian'
            obj.V(ind_1)=   obj.V(ind_1) .* exp(-(x(ind_1)-min(x)-L).^2 * .5*9/L^2);
            obj.V(ind_2)=   obj.V(ind_2) .* exp(-(x(ind_2)-max(x)+L).^2 * .5*9/L^2);
        otherwise
            warning('DFT:QMol_DFT_Vh_fft:falloffShape', ...
                ['Unknown falloff shape ' obj.FalloffUnits '. Falloff ignored.'])
    end

    % Fourier-transform convolution kernel
    obj.V           =   fft(obj.V(:));
    
    % Miscellaneous
    obj.isInit          =   true;

end
function initializeDerivative(~)
%initializeDerivative no initialization required for the potential
%   derivative

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_DFT_Vx_XX_fft';
    PropNames           =  {'Vee','gU','gL','fU','fS','fL','aFL', ...
                            'interactionPotential','ghostUnit','ghostLength',...
                            'falloffUnit','falloffShape','falloffLength','addFalloffLength'};     
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
%   derivative derivative associated with the member properties.
%   Optionally, provide the potential object(s), where the discretization
%   should be stored, and if the exchange potential should be added to it.
%   
%   MORE COMMENTS

    % Initialization
    if nargin < 5,  isAdd   =   false;
    if nargin < 4,  DV       =   [];         end, end
    
    if isempty(DV)   ||   ~isAdd  
        DV              =   obj.DFT.disc.DFT_allocatePotential(DV);          % If any potential, reset it
    end

    % Return the handle to the XX potential derivative operator
    if obj.DFT.isSpinPol,   DV.add(@(opt,p,isUp) obj.applyPotentialDerivative(opt,p,isUp));
    else,                   DV.add(@(opt,p)      obj.applyPotentialDerivative(opt,p));       end
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
    n                   =   numel(obj.DFT.disc.x);
    Vp                  =   zeros(n,1);

    % Compute result of XX potential
    if obj.DFT.isSpinPol
        % Spin polarized
        if isUp
            % Local variable
            N           =   obj.DFT.occ{1};                                 %#ok<PROPLC> 

            % Parse exchange
            for k = 1:numel(N), if N(k) >= obj.tol                          %#ok<PROPLC,ALIGN> 
                V       =   ifft(fft(conj(obj.KSOup(:,k)).*p,obj.N).*obj.V);%#ok<PROPLC> 
                Vp      =   Vp - N(k)*obj.KSOup(:,k) .* V(1:n);             %#ok<PROPLC> 
            end, end
            
            % Real potential
            if isreal(p) && isreal(obj.KSOup),  Vp      =   real(Vp);       end
        else
            % Local variable
            N           =   obj.DFT.occ{2};                                 %#ok<PROPLC> 

            % Parse exchange
            for k = 1:numel(N), if N(k) >= obj.tol                          %#ok<PROPLC,ALIGN> 
                V       =   ifft(fft(conj(obj.KSOdw(:,k)).*p,obj.N).*obj.V);%#ok<PROPLC> 
                Vp      =   Vp - N(k)*obj.KSOdw(:,k) .* V(1:n);             %#ok<PROPLC> 
            end, end
            
            % Real potential
            if isreal(p) && isreal(obj.KSOdw),  Vp      =   real(Vp);       end
        end
    else
        % Spin restricted (include spin-channel resolved scailing)
        N               =   0.5*obj.DFT.occ;                                %#ok<PROPLC> 

        % Parse exchange
        for k = 1:numel(N), if N(k) >= obj.tol                              %#ok<PROPLC,ALIGN> 
                V       =   ifft(fft(conj(obj.KSO(:,k))  .*p,obj.N).*obj.V);%#ok<PROPLC> 
                Vp      =   Vp - N(k)*obj.KSO(:,k)   .* V(1:n);             %#ok<PROPLC>
        end, end
            
        % Real potential
        if isreal(p)   && isreal(obj.KSO),      Vp      =   real(Vp);       end

    end

    % Finalize results
    Vp                  =   Vp * obj.DFT.disc.dx;

end
function DVp = applyPotentialDerivative(obj,opt,p,isUp)
%applyPotentialDerivative computes the result of the exact-exchange
%   potential applied to the input one-electron wave function p for the
%   selected type of computation opt
    
    switch lower(opt),                                                      case 'dipacc',                                    if obj.DFT.isSpinPol
        % Dipole acceleration
        DVp             =   2*real(sum(obj.applyPotential(p,isUp) .* ifft(obj.DFT.disc.D.*fft(conj(p))))*obj.DFT.disc.dx);    else
        DVp             =   2*real(sum(obj.applyPotential(p     ) .* ifft(obj.DFT.disc.D.*fft(conj(p))))*obj.DFT.disc.dx);    end
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
    n                   =   numel(obj.DFT.disc.x);

    % Compute XX energy
    if obj.DFT.isSpinPol
        % Up-spin component
        N               =   obj.DFT.occ{1};                                 %#ok<PROPLC> 

        for k = 1:numel(N),     if N(k) >= obj.tol                          %#ok<PROPLC,ALIGN> 
            % "Diagonal" term
            p           =   abs(KSO.KSOup(:,k)).^2;                         %#ok<PROPLC> 
            V           =   ifft(fft(p,obj.N).*obj.V);                      %#ok<PROPLC> 
            E           =   E - 0.5*N(k)^2 * real(sum( p .* V(1:n) ));      %#ok<PROPLC> 

            % "Off-diagonal" terms
            for l=k+1:numel(N), if N(k) >= obj.tol                          %#ok<PROPLC,ALIGN> 
                p       =   conj(KSO.KSOup(:,k)) .* KSO.KSOup(:,l);         %#ok<PROPLC> 
                V       =   ifft(fft(conj(p),obj.N).*obj.V);                %#ok<PROPLC> 
                E       =   E - real( N(k)*N(l) * sum( p .* V(1:n) ) );     %#ok<PROPLC> 
            end, end
        end, end

        % Down-spin component
        N               =   obj.DFT.occ{2}; %#ok<PROPLC> 

        for k = 1:numel(N),     if N(k) >= obj.tol                          %#ok<PROPLC,ALIGN> 
            % "Diagonal" term
            p           =   abs(KSO.KSOdw(:,k)).^2;                         %#ok<PROPLC> 
            V           =   ifft(fft(p,obj.N).*obj.V);                      %#ok<PROPLC> 
            E           =   E - 0.5*N(k)^2 * real(sum( p .* V(1:n) ));      %#ok<PROPLC> 

            % "Off-diagonal" terms
            for l=k+1:numel(N), if N(k) >= obj.tol                          %#ok<PROPLC,ALIGN> 
                p       =   conj(KSO.KSOdw(:,k)) .* KSO.KSOdw(:,l);         %#ok<PROPLC> 
                V       =   ifft(fft(conj(p),obj.N).*obj.V);                %#ok<PROPLC> 
                E       =   E - real( N(k)*N(l) * sum( p .* V(1:n) ) );     %#ok<PROPLC> 
            end, end
        end, end
    else
        % Spin restricted
        N               =   0.5*obj.DFT.occ;                                %#ok<PROPLC> 

        for k = 1:numel(N),     if N(k) >= obj.tol                          %#ok<PROPLC,ALIGN> 
            % "Diagonal" term
            p           =   abs(KSO.KSO(:,k)).^2;                           %#ok<PROPLC> 
            V           =   ifft(fft(p,obj.N).*obj.V);                      %#ok<PROPLC> 
            E           =   E - 0.5*N(k)^2 * real(sum( p .* V(1:n) ));      %#ok<PROPLC> 

            % "Off-diagonal" terms
            for l=k+1:numel(N), if N(k) >= obj.tol                          %#ok<PROPLC,ALIGN> 
                p       =   conj(KSO.KSO(:,k)) .* KSO.KSO(:,l);             %#ok<PROPLC> 
                V       =   ifft(fft(conj(p),obj.N).*obj.V);                %#ok<PROPLC> 
                E       =   E - real( N(k)*N(l) * sum( p .* V(1:n) ) );     %#ok<PROPLC> 
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

