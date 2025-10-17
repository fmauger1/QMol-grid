classdef QMol_TDDFT_abs_mask < QMol_suite
%QMol_TDDFT_abs_mask mask absorbing boundaries for TDDFT simulations

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDDFT_abs_mask})
function showInfo
    fprintf( '  * QMol_TDDFT_abs_mask:\n');
    fprintf(['      > Absorbing boundaries\n',...
             '      > Mask function\n']); 
    QMol_TDDFT_abs_mask.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Mask properties
    fprintf('  * Absorbing boundaries                                     mask function\n');
    if isscalar(obj.L),         fprintf('    Length = %s a.u. on both ends\n',num2str(obj.L));
    elseif numel(obj.L) == 2,   fprintf('    Length = %s (left) and %s (right) a.u.\n',num2str(obj.L(1)),num2str(obj.L(2))); end

    if ischar(obj.shape),                       fprintf('    Shape  = %s\n',obj.shape);
    elseif isa(obj.shape,'function_handle'),    fprintf('    Shape  = %s\n',func2str(obj.shape)); end
    ref                 =   {};

    % Version
    obj.version;
 
end
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the absorber with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end
 
    % Evaluate (and display) estimate of memory footprint
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.DFT.disc.x),'real');

    if opt, QMol_DFT_profiler.showMemoryFootprint('Absorbing boundaries (mask)', mem,1);    end
 
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    L                   =   10
end
properties (GetAccess=public,SetAccess=?QMol_suite)
    shape               =   'cos^1/8'
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    length                  % L
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_TDDFT_abs_mask,?QMol_TDDFT})
    % Linked objects
    DFT                     % DFT.disc always defines the domain
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Run-time variables
    W                       % Discretization of the mask
    nKSO                    % number of KSO
end
properties (Constant)
    isCAP               =   false
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % L ~~~~~~~~~~~~~~~~
    function set.length(obj,val),           obj.L       =   val;            end
    function val = get.length(obj),         val         =   obj.L;          end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.DFT             =   [];
    obj.W               =   [];
    obj.nKSO            =   [];
    
    % Initialization status
    obj.isInit          =   false;
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.L),      obj.L       =   10;                             end
    if isempty(obj.shape),  obj.shape   =   'cos^1/8';                      end
end
function initialize(obj,DFT,~)
%initialize initializes the object
    
    % Reinitialize any component
    obj.reset;
    
    % Set links
    obj.DFT             =   DFT;

    % Build the mask
    if ischar(obj.shape), switch lower(obj.shape)                           %#ok<ALIGN> 
        case {'cos^1/8','cos^0.125','cos^.125','sin^1/8','sin^0.125','sin^.125'}
                                            fcn =   @(x) cos(.5*pi*x).^.125;
        case {'cos^2','sin^2'},             fcn =   @(x) cos(.5*pi*x).^2;
        otherwise,                          fcn =   @(x) cos(.5*pi*x).^.125;
            warning('QMol:TDDFT:absorberMaskShape', ...
                ['Unknown absorbing-boundary mask shape ' obj.shape '; Default cos^1/8 used instead.']), end
    elseif isa(obj.shape,'function_handle'),fcn =   obj.shape;
    else,                                   fcn =   @(x) cos(.5*pi*x).^.125;
        warning('QMol:TDDFT:absorberMaskShape', ...
            'Unknown/unsupported absorbing-boundary mask shape; Default cos^1/8 used instead.')
    end

    obj.W               =   ones(numel(obj.DFT.disc.x),1);
    if isscalar(obj.L),     lg  =   obj.L * [1 1];
    else,                   lg  =   obj.L;              end
    ind                 =   obj.DFT.disc.x < obj.DFT.disc.x(1) + lg(1);
    obj.W(ind)          =   fcn(1 - (obj.DFT.disc.x(ind)-obj.DFT.disc.x(1))/lg(1));
    ind                 =   obj.DFT.disc.x > obj.DFT.disc.x(end) - lg(2);
    obj.W(ind)          =   fcn((obj.DFT.disc.x(ind)-obj.DFT.disc.x(end)+lg(2))/lg(2));

    % Miscellaneous components
    if obj.DFT.isSpinPol,   obj.nKSO    =   [numel(obj.DFT.occ{1}),numel(obj.DFT.occ{2})];
    else,                   obj.nKSO    =   numel(obj.DFT.occ);             end
    
    % Update initialization status
    obj.isInit          =   true;

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_TDDFT_abs_mask';
    PropNames           =  {'L','shape','length'};                
end
end
%% Absorber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function applyMask(obj,KSO)
%applyMask
    
    if obj.DFT.isSpinPol
        for k = 1:obj.nKSO(1),  KSO.KSOup(:,k)  =   KSO.KSOup(:,k) .* obj.W;    end
        for k = 1:obj.nKSO(2),  KSO.KSOdw(:,k)  =   KSO.KSOdw(:,k) .* obj.W;    end
    else
        for k = 1:obj.nKSO,     KSO.KSO(  :,k)  =   KSO.KSO(  :,k) .* obj.W;    end
    end
end
function V = getPotential(~,V)
%getPotential no complex potential, return the object untouched
 
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

