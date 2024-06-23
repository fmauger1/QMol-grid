classdef QMol_DFT_Vh_fft < QMol_suite
%QMol_DFT_Vh_fft implementation of the Hartree potential functional for
%   density-functional theory (DFT) simulations discretized on a Cartesian
%   grid. Hartree-functional components are computed using a fast-Fourier-
%   transform based convolution scheme over an extended domain. It is
%   almost always faster than the convolution based scheme but requires 
%   tuning of the extended domain.

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
    fprintf( '  * QMol_DFT_Vh_fft:\n');
    fprintf(['      > Hartree potential\n',...
             '      > Fast-Fourier-transform convolution scheme\n']); 
    QMol_DFT_Vh_fft.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Initialization
    ref                 =   {};

    % Header
    fprintf('  * Hartree-potential functional        fast-Fourier-transform convolution\n');

    % SIC (if any)
    if obj.SIC == 1
        fprintf('    with average-density self-interaction correct. (ADSIC) [Legrand 2002].\n');
        ref             =   [ref, {'Legrand2002'}];
    end

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

    % Display components (do not trust V, or DV)
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(Nb,'imag');
    if opt
        QMol_DFT_profiler.showMemoryFootprint('Hartree functional',mem,1);
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
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_DFT_Vh_fft,?QMol_DFT})
    % Linked objects
    DFT                             % DFT.disc always defines the domain
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Run-time variables
    V
    N
    SIC                             % 0 = none, 1 = ADSIC
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
    obj.SIC             =   [];
    
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
            warning('QMol:QMol_DFT_Vh_fft:SIC', ...
                    ['Unknown flavor of self-interaction correction (SIC) ' SIC '\n; ' ...
                     'SIC ignored for the Hartree potential.'])
            obj.SIC     =   0;
    end
    
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
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_DFT_Vh_fft';
    PropNames           =  {'Vee','gU','gL','fU','fS','fL','aFL', ...
                            'interactionPotential','ghostUnit','ghostLength',...
                            'falloffUnit','falloffShape','falloffLength','addFalloffLength'};                
end
end
%% Get Hartree potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function V = getPotential(obj,rho,V,isAdd)
%getPotential returns the discretization of the Hartree potential
%   associated with the member properties and input one-body density (rho).
%   Optionally, provide the potential object(s), where the discretization
%   should be stored, and if the Hartree potential should be added to it.
    
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

    % Compute Hartree potential
    if rho.isSpinPol
        Vh              =   ifft(fft(rho.rhoUp+rho.rhoDw,obj.N) .* obj.V);
        Vh              =   real(Vh(1:numel(rho.rhoUp))) * obj.DFT.disc.dx;
    else
        Vh              =   ifft(fft(rho.rho,            obj.N) .* obj.V);
        Vh              =   real(Vh(1:numel(rho.rho  ))) * obj.DFT.disc.dx;
    end
    
    % Self-interaction correction (SIC == 0 for no correction)
    if obj.SIC == 1
        % ADSIC
        N               =   sum(obj.DFT.Ntot);                              %#ok<PROPLC> 
        Vh              =   ( (N-1)/N ) * Vh;                               %#ok<PROPLC> 
    end

    % Return Hartree potential
    V.add(Vh);

end
function DV = getPotentialDerivative(obj,dim,rho,DV,isAdd)
%getPotentialDerivative returns the discretization of the derivative of the
%   Hartree potential associated with the input one-body density (rho) and
%   its derivative (D_rho). Optionally, provide the potential object(s),
%   where the discretization should be stored, and if the exchange
%   potential should be added to it.
%   The potential derivative is simply the Hartree potential of the
%   derivative of the density.

    % Initialization
    if nargin < 5,  isAdd   =   false;                                      %#ok<ALIGN> 
    if nargin < 4,  DV      =   [];
    if nargin < 3,  rho     =   [];     end, end, end

    
    if isempty(DV)   ||   ~isAdd  
        DV              =   obj.DFT.disc.DFT_allocatePotentialGradient(DV); % If any potential, reset it
    end
    
    if isempty(rho)
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);                % Use DFT one-body density storage
        rho             =   obj.DFT.rho;
    end
    % Compute Hartree potential gradient
    switch dim
        case 1
            % Gradient is Hartree potential of gradient of density
            if rho.isSpinPol
                DVh     =   ifft(fft(rho.D_rhoUp+rho.D_rhoDw,obj.N) .* obj.V);
                DVh     =   real(DVh(1:numel(rho.rhoUp))) * obj.DFT.disc.dx;
            else
                DVh     =   ifft(fft(rho.D_rho,              obj.N) .* obj.V);
                DVh     =   real(DVh(1:numel(rho.rho  ))) * obj.DFT.disc.dx;
            end

        otherwise
            error('QMol:QMol_DFT_Vh_conv:getPotentialDerivative', ...
                 ['Unexpected dimension (' num2str(dim) ') for Hartree-potential derivative computation.']);
    end

    % Self-interaction correction (SIC == 0 for no correction)
    if obj.SIC == 1
        % ADSIC
        N               =   sum(obj.DFT.Ntot);                              %#ok<PROPLC> 
        DVh             =   ( (N-1)/N ) * DVh;                              %#ok<PROPLC> 
    end

    % Return Hartree potential
    DV.add(dim,DVh);

end
function E = getEnergy(obj,rho)
%getEnergy returns the Hartree energy associated with the member properties
%   and input one-body density (rho). Empty or missing density uses the 
%   one-body-density of the member DFT object (discouraged).
    
    % Initialization
    if nargin < 2,      rho     =   [];                                     end
    if isempty(rho),    rho     =   obj.DFT.getDensity(obj.DFT.rho);        end

    % Compute Hartree potential
    if rho.isSpinPol
        Vh              =   ifft(fft(rho.rhoUp+rho.rhoDw,obj.N) .* obj.V);
        E               =   0.5*sum((rho.rhoUp+rho.rhoDw) .* real(Vh(1:numel(rho.rhoUp))) ) * obj.DFT.disc.dx^2;
    else
        Vh              =   ifft(fft(rho.rho,            obj.N) .* obj.V);
        E               =   0.5*sum( rho.rho              .* real(Vh(1:numel(rho.rho  ))) ) * obj.DFT.disc.dx^2;
    end
    
    % Self-interaction correction (SIC == 0 for no correction)
    if obj.SIC == 1
        % ADSIC
        N               =   sum(obj.DFT.Ntot);                              %#ok<PROPLC> 
        E               =   ( (N-1)/N ) * E;                                %#ok<PROPLC> 
    end
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

