classdef QMol_TDDFT_abs_CAP < QMol_suite
%QMol_TDDFT_abs_CAP complex-absorbing potential for TDDFT simulations

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDDFT_abs_CAP})
function showInfo
    fprintf( '  * QMol_TDDFT_abs_CAP:\n');
    fprintf(['      > Absorbing boundaries\n',...
             '      > Complex absorbing potential\n']); 
    QMol_TDDFT_abs_CAP.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Mask properties
    fprintf('  * Absorbing boundaries                       complex absorbing potential\n');
    if isscalar(obj.L),         fprintf('    Length = %s a.u. on both ends\n',num2str(obj.L));
    elseif numel(obj.L) == 2,   fprintf('    Length = %s (left) and %s (right) a.u.\n',num2str(obj.L(1)),num2str(obj.L(2))); end
                                fprintf('    Amplitude = %s a.u.',num2str(obj.amplitude));

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
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.DFT.disc.x),'complex');

    if opt, QMol_DFT_profiler.showMemoryFootprint('Absorbing boundaries (mask)', mem,1);    end
 
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    L                   =   10
    V0                  =   .5
end
properties (GetAccess=public,SetAccess=?QMol_suite)
    shape               =   'sin^1/8'
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    length                  % L
    amplitude               % V0
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_TDDFT_abs_CAP,?QMol_TDDFT})
    % Linked objects
    DFT                     % DFT.disc always defines the domain
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Run-time variables
    V                       % Discretization of the mask
    nKSO                    % number of KSO
end
properties (Constant)
    isCAP               =   true
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % L ~~~~~~~~~~~~~~~~
    function set.length(obj,val),           obj.L       =   val;            end
    function val = get.length(obj),         val         =   obj.L;          end
    % V0 ~~~~~~~~~~~~~~~
    function set.amplitude(obj,val),        obj.V0      =   val;            end
    function val = get.amplitude(obj),      val         =   obj.V0;         end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.DFT             =   [];
    obj.V               =   [];
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
    if isempty(obj.L),      obj.V0      =   .5;                             end
    if isempty(obj.shape),  obj.shape   =   'sin^1/8';                      end
end
function initialize(obj,DFT,isFwd)
%initialize initializes the object
    
    % Reinitialize any component
    obj.reset;
    
    % Set links
    obj.DFT             =   DFT;

    % Build the CAP
    if ischar(obj.shape), switch lower(obj.shape)                           %#ok<ALIGN> 
        case {'sin^1/8','sin^0.125','sin^.125','cos^1/8','cos^0.125','cos^.125'}
                                            fcn =   @(x) 1-cos(.5*pi*x).^.125;
        case {'sin^2','cos^2'},             fcn =   @(x) 1-cos(.5*pi*x).^2;
        case {'sin','cos'},                 fcn =   @(x) 1-cos(.5*pi*x);
        otherwise,                          fcn =   @(x) 1-cos(.5*pi*x).^.125;
            warning('QMol:TDDFT:absorberMaskShape', ...
                ['Unknown absorbing-boundary mask shape ' obj.shape '; Default cos^1/8 used instead.']), end
    elseif isa(obj.shape,'function_handle'),fcn =   obj.shape;
    else,                                   fcn =   @(x) 1-cos(.5*pi*x).^.125;
        warning('QMol:TDDFT:absorberMaskShape', ...
            'Unknown/unsupported absorbing-boundary mask sahep; Default cos^1/8 used instead.')
    end

    obj.V               =   zeros(numel(obj.DFT.disc.x),1);
    if isscalar(obj.L),     lg  =   obj.L * [1 1];
    else,                   lg  =   obj.L;                                  end
    ind                 =   obj.DFT.disc.x < obj.DFT.disc.x(1) + lg(1);
    obj.V(ind)          =   fcn(1 - (obj.DFT.disc.x(ind)-obj.DFT.disc.x(1))/lg(1));
    ind                 =   obj.DFT.disc.x > obj.DFT.disc.x(end) - lg(2);
    obj.V(ind)          =   fcn((obj.DFT.disc.x(ind)-obj.DFT.disc.x(end)+lg(2))/lg(2));

    % Scale the CAP
    if isFwd,   obj.V   =  -1i*obj.V0*obj.V;
    else,       obj.V   =   1i*obj.V0*obj.V;                                end


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
    
    ClassName           =   'QMol_TDDFT_abs_CAP';
    PropNames           =  {'L','V0','shape','length','amplitude'};
end
end
%% Absorber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function applyMask(~,~)
%applyMask
    
    %nothing to do
end
function V = getPotential(obj,V)
%getPotential no complex potential, return the object untouched
    
    % Add CAP to V
    V.add(obj.V);
 
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

