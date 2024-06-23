classdef QMol_DFT_density < QMol_suite
%QMol_DFT_density density-functional-theory (DFT) one-body density on
%   Cartesian-grid domain discretization.
%
%   NOTES:
%   > One can change the data-equality criterion through the member
%     property eqTol. Side effects to doing so are untested, though (the
%     code is developed with eqTol = 1e-10).
   
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT_density})
function showInfo
    fprintf('  * QMol_DFT_density:\n      > One-body density\n'); 
    QMol_DFT_density.version;
end
end
methods (Access=public)
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components (do not trust rho, rhoUp, or rhoDw, etc. they might
    % have been initialized empty; fetch size from domain and parent DFT)
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.disc.x),'real');

    if opt
        QMol_DFT_profiler.showMemoryFootprint('One-body density', 0,1);     if obj.isSpinPol
        QMol_DFT_profiler.showMemoryFootprint('density'     ,2*mem,2);
        QMol_DFT_profiler.showMemoryFootprint('gradient'    ,2*mem,2);      else
        QMol_DFT_profiler.showMemoryFootprint('density'     ,  mem,2);
        QMol_DFT_profiler.showMemoryFootprint('gradient'    ,  mem,2);      end
    end

    % Return memory footprint
    if obj.isSpinPol,   mem     =   4*mem;
    else,               mem     =   2*mem;      end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Transient,Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Linked objects
    disc                            % discretization object (for gradient)
end
properties (GetAccess=public,SetAccess=?QMol_suite)
    isSpinPol
end
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    rho
    rhoUp
    rhoDw
end
properties (Hidden,Transient,GetAccess=public,SetAccess=?QMol_suite)
    % Gradient
    isGrad              =   false
    D_rho
    D_rhoUp
    D_rhoDw
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    density                         % rho
    densityUp                       % rhoUp
    densityDown                     % rhoDw
end
properties (Constant,Access=private)
    eqTol               =   1e-10   % tolerance for equality of data
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % rho ~~~~~~~~~~~~~~
    function val = get.density(obj),            val         =   obj.rho;    end
    function set.density(obj,val),              obj.rho     =   val;        end
    % rhoUp ~~~~~~~~~~~~
    function val = get.densityUp(obj),          val         =   obj.rhoUp;  end
    function set.densityUp(obj,val),            obj.rhoUp   =   val;        end
    % rhoDw ~~~~~~~~~~~~
    function val = get.densityDown(obj),        val         =   obj.rhoDw;  end
    function set.densityDown(obj,val),          obj.rhoDw   =   val;        end
end
%% Accessing derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
function D = get.D_rho(obj)
% getting D_rho (no spin-polarization check is performed)
    
    % Do we have a correct gradient?
    if ~obj.isGrad,     obj.setGradient(1);     end

    % Return gradient
    D                   =   obj.D_rho;
end
function D = get.D_rhoUp(obj)
% getting D_rhoUp (no spin-polarization check is performed)
    
    % Do we have a correct gradient?
    if ~obj.isGrad,     obj.setGradient(1);     end

    % Return gradient
    D                   =   obj.D_rhoUp;
end
function D = get.D_rhoDw(obj)
% getting D_rhoDw (no spin-polarization check is performed)
    
    % Do we have a correct gradient?
    if ~obj.isGrad,     obj.setGradient(1);     end

    % Return gradient
    D                   =   obj.D_rhoDw;
end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % streamlined clear all
    if nargin == 1
        % Clear gradient
        obj.resetGradient(false);

        % Clear all other properties
        obj.reset;

        obj.disc        =   [];
        obj.isSpinPol   =   [];
        obj.rho         =   [];
        obj.rhoUp       =   [];
        obj.rhoDw       =   [];

    % parent-class selected clear
    else
        clear@QMol_suite(obj,varargin{:});
    end
end
function reset(obj)
%reset resets all temporary (transient) properties of the object

    % Soft reset the derivative(s)
    obj.resetGradient(false);

    % Initialization status
    obj.isInit          =   false;
end
function resetGradient(obj,opt)
%resetDerivative resets the derivative(s) of the density without resetting
%   the object. Can be used to force recomputing the derivative(s) or
%   clearing memory
    
    % Initialization
    if nargin < 2,      opt     =   false;      end

    % Reset the derivatives
    obj.isGrad          =   false;

    if opt,     obj.D_rho   =   [];     obj.D_rhoUp =   [];     obj.D_rhoDw =   []; end
end
function setGradient(obj,dim)
%setGradient computes the density gradient and stores it in the proper
%   D_rho(Up/Dw) property
    
    % Initialization (for support with higher dimension)
    if nargin < 2,      dim     =   1;          end

    % Set gradient
    if numel(dim) ~= 1   ||   dim ~= 1
        warning('QMol:QMol_DFT_density:setGradient', ...
                ['Unexpected dimension (' num2str(dim) ') when setting the density gradient. No gradient set.'])
    else,                                                                   if obj.isSpinPol
        obj.D_rhoUp     =   real(ifft(obj.disc.D .* fft(obj.rhoUp)));
        obj.D_rhoDw     =   real(ifft(obj.disc.D .* fft(obj.rhoDw)));       else
        obj.D_rho       =   real(ifft(obj.disc.D .* fft(obj.rho  )));       end
    end
end
function D = getGradient(obj,dim,isUp)
%getGradient uniform interface for accessing the gradient derivative
    
    if dim == 1,                                                            if ~obj.isSpinPol
            D               =   obj.D_rho;                                  elseif isUp
            D               =   obj.D_rhoUp;                                else
            D               =   obj.D_rhoDw;                                end
    else
        error('QMol:QMol_DFT_density:getGradient', ...
                ['Unexpected dimension (' num2str(dim) ') when accessing the density gradient. No gradient set.'])
    end
end
function initialize(obj,disc)
%initialize initializes the object
    
    % Reinitialize any component
    obj.reset;

    % Update discretization
    if nargin == 1,     obj.disc    =   [];
    else,               obj.disc    =   disc;                               end
    
    % Update initialization status
    obj.isInit          =   true;
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_DFT_density';
    PropNames           =  {'isSpinPol','rho','rhoUp','rhoDw', ...
                            'density','densityUp','densityDown'};
end
end
%% Operator overload %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function t = eq(data1,data2)
%eq overloads the == operator 

    % Test domain
    if      isempty(data1.disc)   ||   isempty(data2.disc),     t   =   true;            % not enough information on discretizations
    else,   t   =   data1.disc == data2.disc;                               end

    % Test actual data
    if data1.isSpinPol
        t               =   t   &&   data2.isSpinPol                                && ... both spin polarized
                            all(size(data1.rhoUp) == size(data2.rhoUp))             && ... same number of data
                            all(size(data1.rhoDw) == size(data2.rhoDw))             && ...
                            all(abs(data1.rhoUp-data2.rhoUp) <= data1.eqTol,'all')  && ... all elements are the same (within tolerance)
                            all(abs(data1.rhoDw-data2.rhoDw) <= data1.eqTol,'all');
    else
        t               =   t   &&  ~data2.isSpinPol                                && ... both spin restricted
                            all(size(data1.rho) == size(data2.rho))                 && ... same number of data
                            all(abs(data1.rho-data2.rho) <= data1.eqTol,'all');          % all elements are the same (within tolerance)
    end
end
function t = ne(data1,data2)
%ne overloads the ~= operator
    t                   =   ~(data1.eq(data2));
end
end
%% QMol-grid package methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_suite)
function getDensity(obj,occ,KSO)
%getDensity builds the one-body density for KOhn-Sham orbitals and
%   occupation coefficients

    % For basis-set models, reconstruct the density
    if KSO.isBasis
        KSO.disc.KSO    =   KSO.disc.DFT_reconstructOrbital(KSO,KSO.disc.KSO);
        KSO             =   KSO.disc.KSO;
    end

    if KSO.isSpinPol
        % Compute one-body densities
        obj.rhoUp       =   occ{1}(1) * abs(KSO.KSOup(:,1)).^2;
        obj.rhoDw       =   occ{2}(1) * abs(KSO.KSOdw(:,1)).^2;
        for k = 2:numel(occ{1})
            obj.rhoUp   =   obj.rhoUp + occ{1}(k) * abs(KSO.KSOup(:,k)).^2;
        end
        for k = 2:numel(occ{2})
            obj.rhoDw   =   obj.rhoDw + occ{2}(k) * abs(KSO.KSOdw(:,k)).^2;
        end

        % House keeping
        obj.isSpinPol   =   true;
        obj.rho         =   [];
    else
        % Compute one-body density
        obj.rho         =   occ(1) * abs(KSO.KSO(:,1)).^2;
        for k = 2:numel(occ)
            obj.rho     =   obj.rho + occ(k) * abs(KSO.KSO(:,k)).^2;
        end

        % House keeping
        obj.isSpinPol   =   false;
        obj.rhoUp       =   [];
        obj.rhoDw       =   [];
    end

    % Initialize object (making sure the density and orbital disc match)
    obj.initialize(KSO.disc);
end
function N = getCharge(obj)
%getCharge computes the electronic charge associated with the one-body 
%   density object
    
    if obj.isSpinPol,   N   =  [sum(obj.rhoUp),sum(obj.rhoDw)]*obj.disc.dx;
    else,               N   =   sum(obj.rho  )*obj.disc.dx;                 end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

