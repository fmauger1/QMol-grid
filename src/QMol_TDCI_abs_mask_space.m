classdef QMol_TDCI_abs_mask_space < QMol_TDCI_abs_mask
%QMol_TDCI_abs_mask_space TDCI mask-type absorbing-boundry damping matrix
%   Use QMol_TDCI_abs_mask_space to build mask damping matrix for TDCI 
%   calculations based on the spatial overlap between the CI molecular-
%   orbital basis and a spatial mask. The overlap mask is characterized by
%   (i) the distance from the origin at which it starts, (ii) its width, 
%   and (iii) its shape.
%
%   After each time step in the TDCI propagation, the CI wave function is
%   multiplied by the damping matrix. 
%
%   Editable properties:
%     * distance, length, shape
%
%   Methods:
%     * Changing class properties: set, reset, clear
%     * Run-time documentation: showDocumentation, getMemoryProfile
%     * Damping matrix: getDampingMatrix
%
%   ABS = QMol_TDCI_abs_mask_space('name1','value1',___) creates a mask-
%     type damping matrix with the name properties set to the specified
%     values.
%
%   See also QMol_TDCI

%   Version     Date        Author
%   01.23.000   06/07/2025  F. Mauger
%       Creation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.000','06/07/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDCI_abs_mask})
function showInfo
    fprintf(['  * QMol_TDCI_abs_mask_space:\n', ...
             '      > Damping matrix for mask-type absorbing boundary conditions\n', ...
             '      > Spatial criterion for the damping strength\n']); 
    QMol_TDCI_abs_mask_space.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Mask properties
    fprintf('  * Damping matrix (absorbing boundaries)                    mask function\n');
    if isscalar(obj.D),         fprintf('    Distance  = %s a.u. on both ends\n',num2str(obj.D));
    elseif numel(obj.D) == 2,   fprintf('    Distance  = %s (left) and %s (right) a.u.\n',num2str(obj.D(1)),num2str(obj.D(2))); end
    if isscalar(obj.L),         fprintf('    Length    = %s a.u. on both ends\n',num2str(obj.L));
    elseif numel(obj.L) == 2,   fprintf('    Length    = %s (left) and %s (right) a.u.\n',num2str(obj.L(1)),num2str(obj.L(2))); end

    if ischar(obj.shape),                       fprintf('    Shape  = %s\n',obj.shape);
    elseif isa(obj.shape,'function_handle'),    fprintf('    Shape  = %s\n',func2str(obj.shape)); end
    ref                 =   {};

    % Version
    obj.version;
 
end
function mem = getMemoryProfile(~,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the absorber with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end
 
    % No memory footpring (no large bits of memory stored in the object)
    mem                 =   0;
    if opt, QMol_profiler.showMemoryFootprint('Mask-type absorbing boundary conditions', mem,1);    end
 
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    D
    L                   =   10
end
properties (GetAccess=public,SetAccess=?QMol_suite)
    shape               =   'sin^1/8'   % shape of the mask [ 'sin^1/8' (default) | 'sin^2' | 'sin' | function handle ]
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    distance                % (D) distance from the origin for the beginning of the mask [ positive scalar | positive vector [distLeft distRight] | (default []) ]
    length                  % (L) length of the mask [ positive scalar (default 10) | positive vector [lengthLeft lengthRight] ]
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_TDCI_abs_CAP,?QMol_TDCI})
    % Linked objects
    CI                      % CI.disc always defines the domain
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % D ~~~~~~~~~~~~~~~~
    function set.distance(obj,val),         obj.D       =   val;            end
    function val = get.distance(obj),       val         =   obj.D;          end
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
    obj.CI              =   [];
    
    % Initialization status
    obj.isInit          =   false;
end
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.L),      obj.L       =   10;                             end
    if isempty(obj.shape),  obj.shape   =   'sin^1/8';                      end
end
function initialize(obj,CI)
%initialize initializes the object

    % Initialization
    obj.CI              =   CI;

    % Update initialization status
    obj.isInit          =   true;

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_TDCI_abs_mask_space';
    PropNames           =  {'D','L','distance','shape','length'};
end
end
%% Damping matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function M = getDampingMatrix(obj)
%getDampingMatrix get the damping matrix
%
%   M = obj.getDampingMatrix() returns the damping matrix.

    % Build the CAP mask
    if ischar(obj.shape), switch lower(obj.shape)                           %#ok<ALIGN> 
        case {'sin^1/8','sin^0.125','sin^.125','cos^1/8','cos^0.125','cos^.125'}
                                            fcn =   @(x) cos(.5*pi*x).^.125;
        case {'sin^2','cos^2'},             fcn =   @(x) cos(.5*pi*x).^2;
        case {'sin','cos'},                 fcn =   @(x) cos(.5*pi*x);
        otherwise,                          fcn =   @(x) cos(.5*pi*x).^.125;
            warning('QMol:TDCI:CAPMaskShape', ...
                ['Unknown CAP mask shape ' obj.shape '; Default cos^1/8 used instead.']), end
    elseif isa(obj.shape,'function_handle'),fcn =   obj.shape;
    else,                                   fcn =   @(x) 1-cos(.5*pi*x).^.125;
        warning('QMol:TDCI:CAPMaskShape', ...
            'Unknown/unsupported CAP mask shape; Default cos^1/8 used instead.')
    end

    V                   =   ones(numel(obj.CI.disc.x),1);
    if isscalar(obj.D),     d   =   obj.D * [1 1];
    else,                   d   =   obj.D;                                  end
    if isscalar(obj.L),     lg  =   obj.L * [1 1];
    else,                   lg  =   obj.L;                                  end
    ind                 =   obj.CI.disc.x > -(d(1)+lg(1))   &   obj.CI.disc.x < -d(1);
    V(ind)              =   fcn((-d(1)-obj.CI.disc.x(ind))/lg(1));
    ind                 =   obj.CI.disc.x >  d(2)   &   obj.CI.disc.x <  (d(2)+lg(2));
    V(ind)              =   fcn(( obj.CI.disc.x(ind)-d(2))/lg(2));
    ind                 =   obj.CI.disc.x <= -(d(1)+lg(1))   |   obj.CI.disc.x >=  (d(2)+lg(2));
    V(ind)              =   0;

    % Calculate the damping matrix
    M                   =   sum(V.* (obj.CI.SOB.^2),1)*obj.CI.disc.dx;
    M                   =   prod(M(abs(obj.CI.CSB)),2);
    M                   =   spdiags(M,0,size(obj.CI.CSB,1),size(obj.CI.CSB,1));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

