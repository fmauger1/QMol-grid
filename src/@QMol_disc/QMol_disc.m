classdef QMol_disc < QMol_suite
%QMol_disc Cartesian-grid domain discretization. It defines
%   > The domain discretization grid
%   > Common differential operators on that grid
%
% PROPERTIES
%     x (xspan) | QM | dx | D | T | Tv | nV (basisSize)
%     isBasis
   
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release
%   01.22.000   04/20/2025  F. Mauger
%       Integrale CI module (+ new version number)

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.22.000','04/20/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_disc})
function showInfo
    fprintf('  * QMol_disc: Domain discretization\n      > Cartesian grid discretization\n'); 
    QMol_disc.version;
end
end
methods (Access=?QMol_suite)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Initialization
    fprintf('  * Domain discretization                                   Cartesian grid\n');

    % Domain setting
    if ~obj.isInit
        fprintf('    Domain discretization has not be initialized\n')
    else
        N               =   factor(numel(obj.x));

        fprintf('    Grid = %s:%s:%s\n',num2str(obj.x(1),3),num2str(obj.dx,3),num2str(obj.x(end),3));
        fprintf('    Size = %u ',numel(obj.x));
        if isscalar(N)
            fprintf('(prime) points\n');
        else
            fprintf('(%u',N(1)); fprintf(' x %u',N(2:end)); fprintf(') points\n')
        end
    end

    % Version
    obj.version;

    % References
    ref                 =   {};
end
end 
methods (Access=public)
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components
    mem_real            =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.x),'real');
    mem_imag            =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.x),'imag');
    
    if opt
        QMol_DFT_profiler.showMemoryFootprint('Domain discretization'       ,0       ,1);
        QMol_DFT_profiler.showMemoryFootprint('domain'                      ,mem_real,2);
        QMol_DFT_profiler.showMemoryFootprint('gradient'                    ,mem_imag,2);
        QMol_DFT_profiler.showMemoryFootprint('Laplacian'                   ,mem_real,2);
        QMol_DFT_profiler.showMemoryFootprint('Laplacian (velocity gauge)'  ,mem_imag,2);
    end

    % Return size
    mem                 =   2*mem_real + 2*mem_imag;

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    x
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    xspan                  % x
end
properties (Transient,GetAccess=public,SetAccess={?QMol_disc,?QMol_DFT,?QMol_SEq,?QMol_CI})
    % Linked objects
    QM                     % Parent QMol (theory) object
end
properties (Transient,GetAccess=public,SetAccess=?QMol_disc)
    % Domain
    dx
    % Differential operators
    D                       % Gradient
    T                       % Laplacian
    Tv                      % Laplacian in velocity gauge
end
properties (Hidden,Transient,GetAccess=public,SetAccess=?QMol_disc)
    % Basis size
    nV
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_disc)
    basisSize               % nV
end
methods (Access=public,Static)
    function b = isBasis,   b   =   false;  end
end
properties (Constant,Access=private)
    eqTol               =   1e-10   % tolerance for equality of data
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % x ~~~~~~~~~~~~~~~~
    function set.xspan(obj,val),            obj.x       =   val;            end
    function val = get.xspan(obj),          val         =   obj.x;          end
    % nV ~~~~~~~~~~~~~~~
    function set.basisSize(obj,val),        obj.nV      =   val;            end
    function val = get.basisSize(obj),      val         =   obj.nV;         end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.QM              =   [];
    obj.dx              =   [];
    obj.D               =   [];
    obj.T               =   [];
    obj.Tv              =   [];
    obj.nV              =   [];
    
    % Initialization status
    obj.isInit          =   false;
end
function initialize(obj,QM)
%initialize initializes the object
    
    % Update links
    if nargin < 2,  QM  =   [];                 end
    obj.QM              =   QM;

    % Initialization 
    if obj.isInit,   return; end

    % Domain
    obj.dx              =   obj.x(2)-obj.x(1);

    % Differential operators
    obj.D               =   fourierTool.fftGrid(obj.x).';
    obj.T               =   .5*obj.D.^2;        % Tv is time dependent
    obj.D               =   1i*obj.D;
    obj.nV              =   numel(obj.x);
    
    % Miscellaneous
    obj.isInit          =   true;
end
function setTv(obj,A)
%setTv sets the Laplacian (kinetic) operator Tv, in the velocity gauge, for
%   the input potential vector A. Empty A clears Tv (reverting back to the
%   default T operator).
    
    % Initialization
    if nargin < 2,  A       =   [];                         end

    % Set Tv
    if isempty(A),  obj.Tv  =   [];
    else,           obj.Tv  =   0.5*(-1i*obj.D+A).^2;       end

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_disc';
    PropNames           =   {'x','xspan'};
end
end
%% Operator overload %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function t = eq(disc1,disc2)
%eq overloads the == operator
    t                   =   numel(disc1.x) == numel(disc2.x)                        && ... same number of elements
                            all(abs(disc1.x(:)-disc2.x(:)) <= disc1.eqTol,'all');        % all elements are the same (within tolerance)
end
function t = ne(disc1,disc2)
%ne overloads the ~= operator
    t                   =   ~(disc1.eq(disc2));
end
end
%% DFT specific methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public,Static)
    function s = DFT_classOrbital,              s   =   'QMol_DFT_orbital';     end
    function s = DFT_classDensity,              s   =   'QMol_DFT_density';     end
    function s = DFT_classPotential,            s   =   'QMol_DFT_Vks';         end
    function s = DFT_classPotentialGradient,    s   =   'QMol_DFT_Vks_grad';    end
end
methods (Access=public)
    % Data allocation
    rho = DFT_allocateDensity(obj,rho,isFull)
    KSO = DFT_allocateOrbital(obj,N,KSO,randStr,isR)
    Vks = DFT_allocatePotential(obj,Vks)
    DVks = DFT_allocatePotentialGradient(obj,DVks)
    % Miscellaneous
    d = DFT_dist(obj,X1,X2)
    X = DFT_mix(obj,X,varargin)
    p = DFT_normalizeOrbital(obj,p)
    p = DFT_randomOrbital(obj,randStr,S)
    b = DFT_SCF_AndersonMixCoeff(obj,X1i,X2i,X1o,X2o)
    s = DFT_sizeOrbital(obj)
    % Operators and functionals
    E = DFT_energyKinetic(obj,occ,KSO)
    [E,DE] = DFT_energyOrbital(obj,V,KSO)
    Hp = DFT_operatorHamiltonian(obj,V,varargin)
    Hp = DFT_operatorKinetic(obj,p,S)
    Hp = DFT_operatorPotential(obj,V,varargin)
end
methods (Access=private)
function p = DFT_applySymmetry(~,S,p)
%DFT_applySymmetry applies the proper symmetry on object
    
    if S == -1      % Antisymmetrization
        p               =   .5*(p - flip(p,1));
    elseif S == 1   % Symmetrization
        p               =   .5*(p + flip(p,1));
    end
end
end
%% Schrodinger-equation specific methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public,Static)
    function s = SE_classWaveFunction,          s   =   'QMol_SE_wfcn';         end
    function s = SE_classPotential,             s   =   'QMol_SE_V';            end
end
methods (Access=public)
    % Data allocation
    wfcn    =   SE_allocateWaveFunction(obj,N,wfcn,randStr,isR)
    % Miscellaneous
    p       =   SE_normalizeWaveFunction(obj,p)
    p       =   SE_randomWaveFunction(obj,randStr,S)
    s       =   SE_sizeWaveFunction(obj)
    % Operators and functionals
    E       =   SE_energyKinetic(obj,wfcn)
    [E,DE]  =   SE_energyWaveFunction(obj,V,wfcn)
    Hp      =   SE_operatorHamiltonian(obj,V,p,S)
    Hp      =   SE_operatorKinetic(obj,p,S)
    Hp      =   SE_operatorPotential(obj,V,p,S)
end
methods (Access=private)
function p = SE_applySymmetry(~,S,p)
%SE_applySymmetry applies the proper symmetry on object
    
    if S == -1      % Antisymmetrization
        p               =   .5*(p - flip(p,1));
    elseif S == 1   % Symmetrization
        p               =   .5*(p + flip(p,1));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

