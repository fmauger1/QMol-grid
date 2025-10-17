classdef (Abstract,Hidden) QMol_suite < handle
%QMol_suite reference class for the QMol_suite package suite: all classes
%   in the QMol-grid package derive from this class. It defines package-
%   wide features:
%     * Constructors with name-value pair assignment
%     * The set, reset, and clear methods
%
%   Most end-users will never directly interact with this class but rather
%   the various quantum-mechanical computation frameworks the QMol-grid
%   package supports. Those include:
%     * Ground- and excited-state electronic structure, both using a 
%       Cartesian-grid or basis-set discretization:
%       * Density-functional theory (DFT):
%         * spin restricted DFT, with QMol_DFT_spinRes.
%         * spin polarized DFT, with QMol_DFT_spinPol.
%       * Hartree-Fock (HF), obtained by running a DFT calculation with an 
%         exact exchange and no correlation functional.
%       * Schrodinger equation (SE), with QMol_SE.
%     * Time-dependent electron dynamics, currently restricted to 
%       Cartesian-grid discretization:
%       * Real-time time-dependent DFT (TDDFT), with the QMol_TDDFT suite
%         of propagators.
%         Note: only local (LDA and GGA type) fuctionals are currently
%         supported by QMol-grid's TDDFT propagators
%       * Time-dependent Schrodinger equation (TDSE), with the QMol_TDSE
%       suite of propagators.
%
%   Editable properties: N/A
%
%   Methods:
%     * Changing the class properties: set, reset, clear
%
%   See also QMol_DFT_spinRes, QMol_DFT_spinPol, QMol_SE, QMol_TDDFT, 
%       QMol_TDSE

% 2-Clause BSD License
%
% Copyright (c) 2024, Francois Mauger, all right reserved
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release
%   01.21.001   12/08/2024  F. Mauger
%       Clean docstring and hide unused methods inherited from handle

%% Display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Hidden,Static,Access=?QMol_doc)
function showKernelHeader
%showKernelHeader for the documentation, displays the kernel info in the
%   header
    
    QMol_doc.showVersion('01.21.001','12/08/2024','F. Mauger')
end
end
%% Object initialization status %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,Transient,GetAccess=public,SetAccess=?QMol_suite)
    % Housekeeping
    isInit              =   false
end
properties (Dependent,GetAccess=public)
    isInitialized                   % (isInit) whether the object is initialized (true) or not (false)
end
methods
    % isInit ~~~~~~~~~~~
    function val = get.isInitialized(obj),      val         =   obj.isInit;     end
    function set.isInitialized(obj,val),        obj.isInit  =   val;            end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access = public)
function obj = QMol_suite(varargin)         % =============================
%QMol_suite default constructor for the QMol-grid package. Optionally,
%   specify name properties with their corresponding values. Several name-
%   value pairs can be specified consecutively, and names are case 
%   insensitive
%       QMol_suite()
%       QMol_suite('name',value)                   
%       QMol_suite(name1,value1,name2,value2,___)
%   See the list of editable properties in the class description for valid
%   names (again, all case insensitive).

%   NOTES:
%   * Subclass inherit the name-value pair assignment, provided they define
%     a propertyName method with the list of assignable properties.
    
    % If any property input, pass them to the set member method
    if nargin > 1,  obj.set(varargin{:});   end
end
function reset(obj)                         % =============================
%reset changes the object isInitialized property to false.
%   Internally, the reset method also clears all temporary (transient) 
%   properties of the object.

%   NOTES:
%   * Should be overloaded by each subclasses to perform proper reset
%     actions.
%   * Don't forget to update the isInit property to false
    
    obj.isInit          =   false;
end
function set(obj,varargin)                  % =============================
%set updates the name properties of a QMol-grid object to the specified 
%   value. Several name-value pairs can be specified consecutively and
%   names are case insensitive
%       obj.set('name',value)                   
%       obj.set(name1,value1,name2,value2,___)
%   See the list of editable properties in the class description for valid
%   names (again, all case insensitive)

    % Initialization
    obj.reset;
    [CC,CN]             =   obj.propertyNames;

    % Test number of inputs
    if mod(nargin,2) == 0
        warning(['QMolGrid:' CC ':missingPropertyValue'], ...
                'Inputs should be of name-value pair type. Lase entry ignored.')
    end
    
    % Set properties
    for k = 1:floor(.5*nargin)
        % Identify property
        l               =   find(strcmpi(varargin{2*k-1},CN));
        
        % Set property to proper value
        if ~isempty(l)
            obj.(CN{l}) =   varargin{2*k};
        else
            warning(['QMolGrid:' CC ':propertyName'], ...
                ['Unknown property ' varargin{2*k-1} ' -- entry ignored'])
        end
    end
    
end
function clear(obj,varargin)                % =============================
%clear clears all or selected properties.
%   As a side effect, it also runs the object's reset method.
%   
%   obj.clear() clears all the editable properties listed in the class 
%   description
%
%   obj.clear('name') and obj.clear('name1','name2',___) selectively clears
%   the specified name properties. Names can be any of the editable
%   properties listed in the class description and are case insensitive.
    
    % Initialization
    obj.reset;
    [CC,CN]             =   obj.propertyNames;

    if nargin == 1
        % Clear all member properties
        cl              =   CN;
    else
        % Clear listed member properties
        cl              =   varargin;
    end

    % Clear properties
    for k = 1:numel(cl)
        % Identify property
        l               =   find(strcmpi(cl{k},CN));
        
        % Set property to proper value
        if ~isempty(l)
            obj.(CN{l}) =   [];
        else
            warning(['QMolGrid:' CC ':propertyName'], ...
                ['Unknown property ' cl{k} ' -- entry ignored'])
        end
    end

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that are editable
%   through name-value assignment
    
    ClassName           =   [];
    PropNames           =   {};
end
end
%% Remove methods inherited from handle from the help/documentation %%%%%%%
methods (Hidden)
    addlistener                                                             %#ok<NOIN>
    eq                                                                      %#ok<NOIN>
    findobj                                                                 %#ok<NOIN>
    findprop                                                                %#ok<NOIN>
    ge                                                                      %#ok<NOIN>
    gt                                                                      %#ok<NOIN>
    le                                                                      %#ok<NOIN>
    listener                                                                %#ok<NOIN>
    lt                                                                      %#ok<NOIN>
    ne                                                                      %#ok<NOIN>
    notify                                                                  %#ok<NOIN>
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

