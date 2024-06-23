classdef (Abstract,Hidden) QMol_suite < handle
%QMol_suite reference class for the QMol_suite package suite. It defines
%   > Name-value pair assignment 
%   > Set access control
%   > Header, footer, and funding

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

%% Display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Hidden,Static,Access=?QMol_doc)
function showKernelHeader
%showKernelHeader for the documentation, display the kernel info in the
%   header
    
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
%% Object initialization status %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,Transient,GetAccess=public,SetAccess=?QMol_suite)
    % Housekeeping
    isInit              =   false
end
properties (Dependent,GetAccess=public)
    isInitialized                   % isInit
end
methods
    % isInit ~~~~~~~~~~~
    function val = get.isInitialized(obj),      val         =   obj.isInit;     end
    function set.isInitialized(obj,val),        obj.isInit  =   val;            end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access = public)
function obj = QMol_suite(varargin)         % =============================
%QMol_suite default constructor for the QMol_suite toolbox. All subclass 
%   inherit the name-value pair assignation, provided they define a 
%   propertyName method with the list of assignable properties
    
    % If any property input, pass them to the set member method
    if nargin > 1,  obj.set(varargin{:});   end
end
function reset(obj)                         % =============================
%reset clears all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
%
%   NOTE: Don't forget to update the isInit property to false
    
    obj.isInit          =   false;
end
function set(obj,varargin)                  % =============================
%set sets named member properties to defined values

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
%clear clears all or selected member properties
    
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
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   [];
    PropNames           =   {};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

