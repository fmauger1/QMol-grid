classdef QMol_extField_dipole < QMol_suite
%QMol_extField_dipole external field, in the dipole approximation

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_extField_dipole})
function showInfo
    fprintf('  * QMol_extField_dipole:\n');
    fprintf('      > External driving field\n');
    fprintf('      > Dipole approximation\n');
    QMol_extField_dipole.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Header
    fprintf('  * External driving field                            dipole approximation\n');

    % Potential vector
    fprintf('    Potential vector:          ');
    if isa(obj.A,'function_handle'),        fprintf('function\n');          %#ok<ALIGN> 
    elseif isa(obj.A,'griddedInterpolant'), fprintf('interpolation\n')
    else,                                   fprintf('N/A\n');               end

    % Electric field
    fprintf('    Electric field:            ');
    if isa(obj.E,'function_handle'),        fprintf('function\n');          %#ok<ALIGN> 
    elseif isa(obj.E,'griddedInterpolant'), fprintf('interpolation\n')
    else,                                   fprintf('N/A\n');               end

    % Electric field derivative
    fprintf('    Electric field derivative: ');
    if isa(obj.DE,'function_handle'),       fprintf('function\n');          %#ok<ALIGN> 
    elseif isa(obj.DE,'griddedInterpolant'),fprintf('interpolation\n')
    else,                                   fprintf('N/A\n');               end

    % Finalize
    ref                 =   {};
    obj.version;

end
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    mem                 =   0;
    if opt
        QMol_DFT_profiler.showMemoryFootprint('External driving field (dipole approx.)',0,1);
    end

    % Potential vector
    if isa(obj.A,'griddedInterpolant')
        m               =   numel(obj.A.GridVectors{1}) + numel(obj.A.Values);
        m               =   QMol_DFT_profiler.getMemoryFootprint(m,'real');
        mem             =   mem + m;
        QMol_DFT_profiler.showMemoryFootprint('Potential vector',m,2);
    end

    % Electric field
    if isa(obj.E,'griddedInterpolant')
        m               =   numel(obj.E.GridVectors{1}) + numel(obj.E.Values);
        m               =   QMol_DFT_profiler.getMemoryFootprint(m,'real');
        mem             =   mem + m;
        QMol_DFT_profiler.showMemoryFootprint('Electric field',m,2);
    end

    % Electric field derivative
    if isa(obj.DE,'griddedInterpolant')
        m               =   numel(obj.DE.GridVectors{1}) + numel(obj.DE.Values);
        m               =   QMol_DFT_profiler.getMemoryFootprint(m,'real');
        mem             =   mem + m;
        QMol_DFT_profiler.showMemoryFootprint('Electric field derivative',m,2);
    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    A
    E
    DE
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    potentialVector                 % A
    electricField                   % E
    electricFieldDerivative         % DE
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_extField_dipole,?QMol_TDDFT,?QMol_TDSE})
    % Linked objects
    disc
end
properties(Constant,Access=public)
    type                =   'dipole'
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % A ~~~~~~~~~~~~~~~~
    function set.potentialVector(obj,val),          obj.A       =   val;    end
    function val = get.potentialVector(obj),        val         =   obj.A;  end
    % E ~~~~~~~~~~~~~~~~
    function set.electricField(obj,val),            obj.E       =   val;    end
    function val = get.electricField(obj),          val         =   obj.E;  end
    % DE ~~~~~~~~~~~~~~~
    function set.electricFieldDerivative(obj,val),  obj.DE      =   val;    end
    function val = get.electricFieldDerivative(obj),val         =   obj.DE; end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.disc            =   [];
    
    % Initialization status
    obj.isInit          =   false;
end
function initialize(obj,disc)
%initialize initializes the object

    % Initialization
    obj.reset;
    
    % Set links
    if nargin == 2,         obj.disc    =   disc;                           end
    
    % Miscellaneous
    obj.isInit          =   true;

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_extField_dipole';
    PropNames           =   {'A','E','DE','potentialVector','electricField','electricFieldDerivative'};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

