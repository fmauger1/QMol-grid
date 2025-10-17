classdef QMol_TDCI_abs_CAP_space < QMol_TDCI_abs_CAP
%QMol_TDCI_abs_CAP_space CAP-type absorbing-boundry damping matrix for TDCI
%   Use QMol_TDCI_abs_CAP_space to build complex-absorbing potential (CAP)
%   damping matrix for TDCI calculations based on the spatial overlap
%   between the CI molecular-orbital basis and a spatial mask. The overlap
%   mask is characterized by (i) the distance from the origin at which it
%   starts, (ii) its width, (iii) its maximum amplitude, and (iv) its
%   shape.
%
%   The CAP damping matrix is added to the CI matrix and thus automatically
%   scales with the propagation time step. 
%
%   Editable properties:
%     * distance, length, amplitude, shape
%
%   Methods:
%     * Changing class properties: set, reset, clear
%     * Run-time documentation: showDocumentation, getMemoryProfile
%     * Damping matrix: getDampingMatrix
%
%   ABS = QMol_TDCI_abs_CAP_space('name1','value1',___) creates a CAP-type
%     damping matrix with the name properties set to the specified values.
%
%   See also QMol_TDCI

%   Version     Date        Author
%   01.23.000   06/06/2025  F. Mauger
%       Creation
%   01.23.001   06/07/2025  F. Mauger
%       Fix overlap with density (not orbital)

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.001','06/07/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDCI_abs_CAP})
function showInfo
    fprintf(['  * QMol_TDCI_abs_CAP_space:\n', ...
             '      > Damping matrix for CAP-type absorbing boundary conditions\n', ...
             '      > Spatial criterion for the damping strength\n']); 
    QMol_TDCI_abs_CAP_space.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Mask properties
    fprintf('  * Damping matrix (absorbing boundaries)      complex-absorbing potential\n');
    if isscalar(obj.D),         fprintf('    Distance  = %s a.u. on both ends\n',num2str(obj.D));
    elseif numel(obj.D) == 2,   fprintf('    Distance  = %s (left) and %s (right) a.u.\n',num2str(obj.D(1)),num2str(obj.D(2))); end
    if isscalar(obj.L),         fprintf('    Length    = %s a.u. on both ends\n',num2str(obj.L));
    elseif numel(obj.L) == 2,   fprintf('    Length    = %s (left) and %s (right) a.u.\n',num2str(obj.L(1)),num2str(obj.L(2))); end
                                fprintf('    Amplitude = %s a.u.',num2str(obj.amplitude));

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
    if opt, QMol_profiler.showMemoryFootprint('CAP-type absorbing boundary conditions', mem,1);    end
 
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    D
    L                   =   10
    V0                  =   .5
end
properties (GetAccess=public,SetAccess=?QMol_suite)
    shape               =   'sin^1/8'   % shape of the CAP mask [ 'sin^1/8' (default) | 'sin^2' | 'sin' | function handle ]
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    distance                % (D) distance from the origin for the beginning of the CAP mask [ positive scalar | positive vector [distLeft distRight] | (default []) ]
    length                  % (L) length of the CAP mask [ positive scalar (default 10) | positive vector [lengthLeft lengthRight] ]
    amplitude               % (V0) maximum amplitude of the CAP mask [ positive scalar (default 0.5) ]
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
    if isempty(obj.V0),     obj.V0      =   .5;                             end
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
    
    ClassName           =   'QMol_TDCI_abs_CAP_space';
    PropNames           =  {'D','L','V0','distance','shape','length','amplitude'};
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
                                            fcn =   @(x) 1-cos(.5*pi*x).^.125;
        case {'sin^2','cos^2'},             fcn =   @(x) 1-cos(.5*pi*x).^2;
        case {'sin','cos'},                 fcn =   @(x) 1-cos(.5*pi*x);
        otherwise,                          fcn =   @(x) 1-cos(.5*pi*x).^.125;
            warning('QMol:TDCI:CAPMaskShape', ...
                ['Unknown CAP mask shape ' obj.shape '; Default cos^1/8 used instead.']), end
    elseif isa(obj.shape,'function_handle'),fcn =   obj.shape;
    else,                                   fcn =   @(x) 1-cos(.5*pi*x).^.125;
        warning('QMol:TDCI:CAPMaskShape', ...
            'Unknown/unsupported CAP mask shape; Default cos^1/8 used instead.')
    end

    V                   =   zeros(numel(obj.CI.disc.x),1);
    if isscalar(obj.D),     d   =   obj.D * [1 1];
    else,                   d   =   obj.D;                                  end
    if isscalar(obj.L),     lg  =   obj.L * [1 1];
    else,                   lg  =   obj.L;                                  end
    ind                 =   obj.CI.disc.x > -(d(1)+lg(1))   &   obj.CI.disc.x < -d(1);
    V(ind)              =   fcn((-d(1)-obj.CI.disc.x(ind))/lg(1));
    ind                 =   obj.CI.disc.x >  d(2)   &   obj.CI.disc.x <  (d(2)+lg(2));
    V(ind)              =   fcn(( obj.CI.disc.x(ind)-d(2))/lg(2));
    ind                 =   obj.CI.disc.x <= -(d(1)+lg(1))   |   obj.CI.disc.x >=  (d(2)+lg(2));
    V(ind)              =   1;

    V                   =   obj.V0*V;

    % Calculate the damping matrix
    M                   =   sum(V.*(obj.CI.SOB.^2),1)*obj.CI.disc.dx;
    M                   =   sum(M(abs(obj.CI.CSB)),2);
    M                   =   spdiags(M,0,size(obj.CI.CSB,1),size(obj.CI.CSB,1));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

