classdef QMol_DFT_Vks_basis < QMol_DFT_Vks
%QMol_DFT_Vks Kohn-Sham (KS) potential operator for basis-set 
%   discretization. It defines
%   > Common interface for orbitals-dependent and independent KS potentials
%   > The potential discretization itself
%    
%   NOTES (for developers)
%   > The discretization object is accessible through the DFT member
%     property (DFT.disc). This is to avoid discretization mismatch if the
%     DFT object is updated.
%   > Only the explicit (with respect to the one-body density) part of the
%     KS potential is discretized. Implicit components (e.g., exact
%     exchange) only exist as an (linear) operator and cannot be
%     discretized on the grid.

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT_Vks})
function showInfo
    fprintf('  * QMol_DFT_Vks_basis:\n      > Kohn-Sham potential\n      > Basis-set discretization\n'); 
    QMol_DFT_Vks_basis.version;
end
end
methods (Access=public)
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components (do not trust V, Vup, or Vdw, they might
    % have been initialized empty; fetch size from domain and parent DFT)
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(size(obj.disc.v,2)^2,'real');

    if opt,                                                                     if obj.isSpinPol
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham potential',2*mem,1);   else
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham potential',  mem,1);   end
    end

    % Return memory footprint
    if obj.isSpinPol,   mem     =   2*mem;      end

end
end
% %% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Transient,GetAccess=public,SetAccess=private)
    % Matrix representation for the (explicit part of the) KS potential
    mV
    mVup
    mVdw
end
methods (Access=public,Static)
    function b = isBasis,   b   =   true;  end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object
    
    % Run-time variables
    obj.mV              =   [];
    obj.mVup            =   [];
    obj.mVdw            =   [];

    % Initialization status
    obj.isInit          =   false;
end
function initialize(obj,disc)
%initialize initializes the object
    
    % Call parent-class initialization (sets isInit to true)
    initialize@QMol_DFT_Vks(obj,disc);
    
    % Matrix representation for the potential
    if obj.isSpinPol
        % Spin polarized
        obj.mV          =   [];
        obj.mVup        =   NaN(disc.nV,disc.nV);
        obj.mVdw        =   NaN(disc.nV,disc.nV);

        if isreal(disc.v), for k = 1:disc.nV, for l = 1:k                   %#ok<ALIGN> 
            obj.mVup(k,l)=  sum(     disc.v(:,k) .* obj.Vup.*disc.v(:,l)) * disc.dx;   
            obj.mVdw(k,l)=  sum(     disc.v(:,k) .* obj.Vdw.*disc.v(:,l)) * disc.dx;   
            obj.mVup(l,k)=  obj.mVup(k,l);      obj.mVdw(l,k)   =   obj.mVdw(k,l);
        end, end, else ,    for k = 1:disc.nV, for l = 1:k                  %#ok<ALIGN> 
            obj.mVup(k,l)=  sum(conj(disc.v(:,k)).* obj.Vup.*disc.v(:,l)) * disc.dx;   
            obj.mVdw(k,l)=  sum(conj(disc.v(:,k)).* obj.Vdw.*disc.v(:,l)) * disc.dx;   
            obj.mVup(l,k)=  conj(obj.mVup(k,l));    obj.mVdw(l,k)   =   conj(obj.mVdw(k,l));
        end, end, end
    else
        % Spin restricted
        obj.mV          =   NaN(disc.nV,disc.nV);
        obj.mVup        =   [];     obj.mVdw    =   [];

        if isreal(disc.v), for k = 1:disc.nV, for l = 1:k                   %#ok<ALIGN> 
            obj.mV(k,l) =   sum(     disc.v(:,k) .* obj.V.*disc.v(:,l)) * disc.dx;  obj.mV(l,k) =   obj.mV(k,l);
        end, end, else ,    for k = 1:disc.nV, for l = 1:k                  %#ok<ALIGN> 
            obj.mV(k,l) =   sum(conj(disc.v(:,k)).* obj.V.*disc.v(:,l)) * disc.dx;  obj.mV(l,k) =   conj(obj.mV(k,l));
        end, end, end
    end
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    [~,PropNames]       =   QMol_DFT_Vks.propertyNames;
    ClassName           =   'QMol_DFT_Vks_basis';
end
end
% %% Arithmetic on potentials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
% function add(obj,varargin)
function Hp = applyPotential(obj,p,isUp)
%applyPotential potential operator V | psi >
%   The potential is assumed to have been properly initialized, but no
%   check for that is performed (try applying it no matter what).
 
    % Explicit part
    if obj.isSpinPol                                                        %#ok<ALIGN> 
        if isUp,    if isempty(obj.Vup),    Hp  =   0*p;    else,   Hp  =   obj.mVup * p;   end
        else,       if isempty(obj.Vdw),    Hp  =   0*p;    else,   Hp  =   obj.mVdw * p;   end, end
    else,           if isempty(obj.V),      Hp  =   0*p;    else,   Hp  =   obj.mV   * p;   end, end

    % Implicit part
    if ~isempty(obj.Vimp)
        % Reconstruct wave function
        P               =   zeros(numel(obj.disc.x),1);                         % zeros(size(obj.V));
        for k = 1:obj.disc.nV,  P   =   P + p(k) * obj.disc.v(:,k);     end

        % Apply implicit potential(s)
        H               =   zeros(numel(obj.disc.x),1);                         % zeros(size(obj.V));
        for k = 1:numel(obj.Vimp)
            if obj.isSpinPol,               H   =   H  + obj.Vimp{k}(P,isUp);
            else,                           H   =   H  + obj.Vimp{k}(P);                    end
        end

        % Project implicit bit
        if isreal(obj.disc.v),  for k = 1:obj.disc.nV,   Hp(k)   =   Hp(k) + sum(     obj.disc.v(:,k) .*H)*obj.disc.dx;  end
        else,                   for k = 1:obj.disc.nV,   Hp(k)   =   Hp(k) + sum(conj(obj.disc.v(:,k)).*H)*obj.disc.dx;  end, end
    else
        
    end

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

