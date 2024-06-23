classdef QMol_Va_Gaussian < QMol_suite
%QMol_Va_Gaussian atomic pseudopotential with an Gaussian shape (short
%   range) parameterized as -V0*exp(-(x-X0).^2/2/s^2). It defines
%   > parameters
%     name      name of the atomic center (string)
%               AKA: atom
%     V0        potential depth
%               AKA: potentialDepth
%     s         potential width
%               AKA: potentialWidth
%     X0        (central) position of the atomic center
%               AKA: position
%     m         mass (for model including nuclear degrees of freedom)
%               AKA: mass
%     p         momentum (for model including nuclear degrees of freedom)
%               AKA: momentum
   
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_Va_Gaussian})
function showInfo
    fprintf( '  * QMol_Va_Gaussian:\n');
    fprintf(['      > Gaussian-shape atomic pseudopotential\n'...
             '      > Short range potential\n']); 
    QMol_Va_Gaussian.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj,opt)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Initialization
    if nargin == 1,     opt     =   [];         end
    if isempty(opt),    opt     =   'full';     end

    ref                 =   {};
    
    % Display (appropriate) documentation
    switch lower(opt)
        case {'full','all'}
            % Display both potential shape and parameters
            fprintf('  * Atomic center                                     (Gaussian potential)\n');
            fprintf('    Parameterized as V(x) = -V0 * exp( -(x-X0).^2 / (2*s^2) ), with\n');
            if ~isempty(obj.name),  fprintf('    name     = %s\n',obj.name); end
            fprintf('    V0       = %-5.2f \n    s        = %-5.2f \n    X0       = %-6.2f\n',obj.V0,obj.s,obj.X0);
            if ~isempty(obj.m),     fprintf('    mass     = %-9.2f\n',obj.m); end
            if ~isempty(obj.p),     fprintf('    momentum = %-9.2f\n',obj.p); end
            obj.version;
        case {'param','parameters'}
            % Compile name name and nuclear properties
            msg         =   [];
            if ~isempty(obj.m),     msg     =   [msg ' | m = ' num2str(obj.m,'%8.2f')]; end
            if ~isempty(obj.p),     msg     =   [msg ' | p = ' num2str(obj.p,'%8.2f')]; end
            
            if isempty(obj.name), a = '???'; else, a = obj.name; end

            if isempty(msg)
                msg     =   a;
            else
                msg     =   [a ' (' msg(4:end) ')'];
            end

            fprintf('    > %-57s (Gaussian)\n',[msg ', parameterized as'])
            fprintf('      V0 = %5.2f | s = %5.2f | X0 = %6.2f\n',obj.V0,obj.s,obj.X0);

            ref         =   'Gaussian';
        case {'pot','potential'}
            fprintf('  * Gaussian-shape potential                                    (Gaussian)\n')
            fprintf('    Parameterized as V(x) = -V0 * exp( -(x-X0).^2 / (2*s^2) )\n');
            obj.version;
        otherwise
            warning('QMolGrid:Va_Gaussian:documentationCase', ...
                ['Unknown documentation case ' opt ' for Gaussian atomic potential.\nSkipped the documentation output.'])
    end
    
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    name
    V0                  =   1
    s                   =   2
    X0                  =   0
    m
    p
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    atom                            % name
    potentialDepth                  % V0
    potentialWidth                  % s
    position                        % X0
    mass                            % m
    momentum                        % p
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % name ~~~~~~~~~~~~~
    function set.atom(obj,val),                 obj.name    =   val;        end
    function val = get.atom(obj),               val         =   obj.name;   end
    % V0 ~~~~~~~~~~~~~~~
    function set.potentialDepth(obj,val),       obj.V0      =   val;        end
    function val = get.potentialDepth(obj),     val         =   obj.V0;     end
    % s ~~~~~~~~~~~~~~~~
    function set.potentialWidth(obj,val),       obj.s       =   val;        end
    function val = get.potentialWidth(obj),     val         =   obj.s;      end
    % X0 ~~~~~~~~~~~~~~
    function set.position(obj,val),             obj.X0      =   val;        end
    function val = get.position(obj),           val         =   obj.X0;     end
    % m ~~~~~~~~~~~~~~~
    function set.mass(obj,val),                 obj.m       =   val;        end
    function val = get.mass(obj),               val         =   obj.m;      end
    % p ~~~~~~~~~~~~~~~~
    function set.momentum(obj,val),             obj.p       =   val;        end
    function val = get.momentum(obj),           val         =   obj.p;      end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.V0),         obj.V0      =   1;      end
    if isempty(obj.s),          obj.s       =   2;      end
    if isempty(obj.x0),         obj.X0      =   0;      end
end
function initialize(obj,~)
%initialize nothing to initialize. For good measure update isInit property
    obj.isInit       =   true;
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_Va_Gaussian';
    PropNames           =  {'name','V0','s','X0','m','p', ...
                            'atom','potentialDepth','potentialWidth', ...
                            'position','mass','momentum'};
end
end
%% Potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function V = getPotential(obj,x)
%getPotential returns the (Gaussian) potential value at querry points
    V                   =  -obj.V0*exp(-.5*(x-obj.X0).^2/obj.s^2);
end
function V = getPotentialDerivative(obj,~,x)
%getPotentialDerivative returns the derivative of the (Gaussian)
%   potential at querry points
    V                   =   obj.V0/obj.s^2*(x-obj.X0).*exp(-.5*(x-obj.X0).^2/obj.s^2);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

