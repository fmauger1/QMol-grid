classdef QMol_DFT_Vx_LDA_exp < QMol_suite
%QMol_DFT_Vx_LDA_exp implementation of the local-density-approximation
%   (LDA) Slater-exchange potential for an exponential electron-electron
%   interaction potential.
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT_Vx_LDA_exp})
function showInfo
    fprintf( '  * QMol_DFT_Vx_LDA_exp:\n');
    fprintf(['      > LDA Slater exchange\n',...
             '      > Exponential electron interactions\n']); 
    QMol_DFT_Vx_LDA_exp.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Header
    fprintf('  * Slater-exchange functional           local-density approximation (LDA)\n');
    fprintf('    for exponential-potential electron-electron interation [Baker 2015],\n');
    fprintf('    parameterized as:\n');
    fprintf('      Vee(x) = %s * exp( - |x| / %s)',num2str(obj.V0,'%9.3g'),num2str(obj.s,'%9.3g'));
    ref                 =   {'Baker 2015'};

    % SIC (if any)
    if obj.SIC == 0
        fprintf('.\n');
    elseif obj.SIC == 1
        fprintf(',\n    with average-density self-interaction correct. (ADSIC) [Legrand 2002].\n');
        ref             =   [ref, {'Legrand 2002'}];
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
        QMol_DFT_profiler.showMemoryFootprint('Exchange functional (LDA exponential)',mem,1);
    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Expotential potential parameters
    V0                  =   1/sqrt(2)
    s                   =   5
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    potentialHeight                 % V0
    potentialWidth                  % s
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_DFT_Vx_LDA_exp,?QMol_DFT})
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
    tol                 =   1e-10   % Threshold density for zero potential
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % V0 ~~~~~~~~~~~~~~~
    function set.potentialHeight(obj,val),          obj.V0      =   val;    end
    function val = get.potentialHeight(obj),        val         =   obj.V0; end
    % s ~~~~~~~~~~~~~~~~
    function set.potentialWidth(obj,val),           obj.s       =   val;    end
    function val = get.potentialWidth(obj),         val         =   obj.s;  end
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
    if isempty(obj.V0),     obj.V0      =   1/sqrt(2);                      end
    if isempty(obj.s),      obj.s       =   5;                              end
end
function initialize(obj,DFT,SIC)
%initialize initializes the object

    % Initialization
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
            warning('QMol:QMol_Vx_LDA_exp:SIC', ...
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
    
    ClassName           =   'QMol_Vx_LDA_exp';
    PropNames           =   {'V0','s','potentialHeight','potentialWidth'};
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

    % Compute the LDA-exchange potential
    if rho.isSpinPol
        % Spin-polarized computation
        V.add(obj.atanFact(2*rho.rhoUp),obj.atanFact(2*rho.rhoDw));
    else
        % Spin-restricted computation
        V.add(obj.atanFact(rho.rho));
    end

    % Self-interaction correction (SIC == 0 for no correction)
    if obj.SIC == 1
        % ADSIC
        if rho.isSpinPol
            % Spin polarized
            V.add(-obj.atanFact(2/obj.DFT.Ntot(1) * rho.rhoUp),...
                  -obj.atanFact(2/obj.DFT.Ntot(2) * rho.rhoDw));
        else
            % Spin restricted
            V.add(-obj.atanFact(2/obj.DFT.Ntot * rho.rho));
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
    
    if isempty(rho)
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);                % Use DFT one-body density storage
        rho             =   obj.DFT.rho;
    end

    % Compute the LDA-exchange potential derivative
    switch dim
        case 1
            % Potential gradient
            if rho.isSpinPol
                % Spin-polarized computation
                DV.add(1,real(ifft(obj.DFT.disc.D .* fft(obj.atanFact(2*rho.rhoUp)))), ...
                         real(ifft(obj.DFT.disc.D .* fft(obj.atanFact(2*rho.rhoDw)))) );
            else
                % Spin-restricted computation
                DV.add(1,real(ifft(obj.DFT.disc.D .* fft(obj.atanFact(rho.rho    )))) );
            end

            % SIC
            if obj.SIC == 1
                % ADSIC
                if rho.isSpinPol
                    % Spin polarized
                    DV.add(1,-real(ifft(obj.DFT.disc.D .* fft(obj.atanFact(2/obj.DFT.Ntot(1) * rho.rhoUp)))),...
                             -real(ifft(obj.DFT.disc.D .* fft(obj.atanFact(2/obj.DFT.Ntot(2) * rho.rhoDw)))) );
                else
                    % Spin restricted
                    DV.add(1,-real(ifft(obj.DFT.disc.D .* fft(obj.atanFact(2/obj.DFT.Ntot    * rho.rho  )))) );
                end
            end

        otherwise
            error('QMol:QMol_DFT_Vx_LDA_exp:getPotentialDerivative', ...
                 ['Unexpected dimension (' num2str(dim) ') for Hartree-potential derivative computation.']);
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
            E           =   sum(rho.rhoUp.*(obj.fullFact(2*rho.rhoUp)-obj.fullFact(2/obj.DFT.Ntot(1)*rho.rhoUp))) + ...
                            sum(rho.rhoDw.*(obj.fullFact(2*rho.rhoDw)-obj.fullFact(2/obj.DFT.Ntot(2)*rho.rhoDw)));
        else
            % Spin restricted
            E           =   sum(rho.rho.*(obj.fullFact(rho.rho)-obj.fullFact(2/sum(obj.DFT.Ntot)*rho.rho)));
        end
    else
        % No SIC
        if rho.isSpinPol
            % Spin polarized
            E           =   sum(rho.rhoUp.*obj.fullFact(2*rho.rhoUp)) + ...
                            sum(rho.rhoDw.*obj.fullFact(2*rho.rhoDw));
        else
            % Spin restricted
            E           =   sum(rho.rho.*obj.fullFact(rho.rho));
        end
    end

    E                   =   E * obj.DFT.disc.dx;
end
end
%% LDA exchange functional components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=private)
function V = atanFact(obj,rho)
%atanFact atan factor for the exchange functional definition
    
    V                   =  -obj.V0/pi * atan(pi*obj.s * rho);
end
function Ex = fullFact(obj,rho)
%fullFact full (exchange-energy per particle) factor for the exchange
%   functional definition
    
    Ex                  =   obj.V0/obj.s*.5/pi^2*log(1+(pi*obj.s*rho).^2)./rho ...
                           -obj.V0/pi * atan(pi*obj.s * rho);
    Ex(rho<obj.tol)     =   0;

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

