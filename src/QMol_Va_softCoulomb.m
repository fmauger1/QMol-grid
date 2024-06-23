classdef QMol_Va_softCoulomb < QMol_suite
%QMol_Va_softCoulomb soft-Coulomb atomic pseudopotential (long range)
%   parameterized as -Z./sqrt((x-X0).^2+a^2). It defines
%   > parameters
%     name      name of the atomic center (string)
%               AKA: atom
%     Z         (effective) charge
%               AKA: charge
%     a         softening parameters
%               AKA: softeningParameter
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
methods (Static,Access={?QMol_doc,?QMol_Va_softCoulomb})
function showInfo
    fprintf( '  * QMol_Va_softCoulomb:\n');
    fprintf(['      > Soft-Coulomb atomic pseudopotential\n' ...
             '      > Long range potential\n']); 
    QMol_Va_softCoulomb.version;
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
            fprintf('  * Atomic center                                 (soft-Coulomb potential)\n');
            fprintf('    Parameterized as V(x) = -Z ./ sqrt( (x-X0).^2 + a^2 ), with\n');
            if ~isempty(obj.name),  fprintf('    name     = %s\n',obj.name); end
            fprintf('    Z        = %-5.2f \n    a        = %-5.2f \n    X0       = %-6.2f\n',obj.Z,obj.a,obj.X0);
            if ~isempty(obj.m),     fprintf('    mass     = %-9.2f\n',obj.m); end
            if ~isempty(obj.p),     fprintf('    momentum = %-9.2f\n',obj.p); end
            obj.version;
        case {'param','parameters'}
            % Compile name name and nuclear properties
            msg         =   [];
            if ~isempty(obj.m),     msg     =   [msg ' | m = ' num2str(obj.m,'%8.2f')]; end
            if ~isempty(obj.p),     msg     =   [msg ' | p = ' num2str(obj.p,'%8.2f')]; end
            
            if isempty(obj.name), b = '???'; else, b = obj.name; end

            if isempty(msg)
                msg     =   b;
            else
                msg     =   [b ' (' msg(4:end) ')'];
            end

            fprintf('    > %-53s (soft Coulomb)\n',[msg ', parameterized as']);
            fprintf('      Z = %5.2f | a = %5.2f | X0 = %6.2f\n',obj.Z,obj.a,obj.X0);

            ref         =   'softCoulomb';
        case {'pot','potential'}
            fprintf('  * Soft-Coulomb potential [Javanainen 1988]                (soft Coulomb)\n');
            fprintf('    Parameterized as V(x) = -Z ./ sqrt( (x-X0).^2 + a^2 ). \n'); 
            obj.version;

            ref         =   {'Javanainen 1988'};
        otherwise
            warning('QMolGrid:Va_softCoulomb:documentationCase', ...
                ['Unknown documentation case ' opt ' for soft-Coulomb atomic potential.\nSkipped the documentation output.'])
    end
    
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    name
    Z                   =   1
    a                   =   1
    X0                  =   0
    m
    p
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    atom                            % name
    charge                          % Z
    softeningParameter              % a
    position                        % X0
    mass                            % m
    momentum                        % p
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % name ~~~~~~~~~~~~~
    function set.atom(obj,val),                 obj.name    =   val;        end
    function val = get.atom(obj),               val         =   obj.name;   end
    % Z ~~~~~~~~~~~~~~~~
    function set.charge(obj,val),               obj.Z       =   val;        end
    function val = get.charge(obj),             val         =   obj.Z;      end
    % a ~~~~~~~~~~~~~~~~
    function set.softeningParameter(obj,val),   obj.a       =   val;        end
    function val = get.softeningParameter(obj), val         =   obj.a;      end
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
    if isempty(obj.Z),          obj.Z       =   1;      end
    if isempty(obj.a),          obj.a       =   1;      end
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
    
    ClassName           =   'QMol_Va_softCoulomb';
    PropNames           =  {'name','Z','a','X0','m','p', ...
                            'atom','charge','softeningParameter', ...
                            'position','mass','momentum'};
end
end
%% Potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function V = getPotential(obj,x)
%getPotential returns the (soft-Coulomb) potential value at querry points
    V                   =  -obj.Z./sqrt((x-obj.X0).^2+obj.a^2);
end
function V = getPotentialDerivative(obj,~,x)
%getPotentialDerivative returns the derivative of the (soft-Coulomb)
%   potential at querry points
    V                   =   obj.Z*(x-obj.X0).*((x-obj.X0).^2+obj.a^2).^-1.5;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

