classdef QMol_progressBar < QMol_suite
%QMol_progressBar management and display of progress bar
%
%   Use QMol_progressBar to manage the display of a progress bar in
%   calculations. This is used internally to provide a unified handling of
%   progress bars and end users are not expected to directly interact with
%   the class.

%   Version     Date        Author
%   01.23.000   05/22/2025  F. Mauger
%       Creation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.000','05/22/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_progressBar})
function showInfo
    fprintf('  * QMol_progressBar:\n      > Progress bar\n'); 
    QMol_progressBar.version;
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,Transient,GetAccess=public,SetAccess=?QMol_suite)
    N
    vMin            =   0
    vMax            =   1
end
properties (GetAccess=public,SetAccess=?QMol_suite)
    motif           =   '|'     % motif for the progress bar [ single character (default '|') ]
    showZero        =   true    % whether to show a bar at 0% progress [ true (default) | false ]
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    barLength               % (N) number of increments in the progress bar [ integer (default []) ]
    startValue              % (vMin) starting value for the progress bar [ scalar (default 0) ]
    endValue                % (vMax) ending value for the progress bar [ scalar (default 1) ]
end
properties (Transient,Access=private)
    iBar                    % index of the last displayed bar
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % N ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function val = get.barLength(obj),              val     =   obj.N;      end
    function set.barLength(obj,val),                obj.N   =   val;        end
    % vMin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function val = get.startValue(obj),             val     =   obj.vMin;   end
    function set.startValue(obj,val),               obj.vMin=   val;        end
    % vMax ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function val = get.endValue(obj),               val     =   obj.vMax;   end
    function set.endValue(obj,val),                 obj.vMax=   val;        end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % Run parent clear
    clear@QMol_suite(obj,varargin{:});

    % Reset default parameters (if needed)
    if isempty(obj.vMin),       obj.vMin        =   0;                      end
    if isempty(obj.vMax),       obj.vMax        =   1;                      end
    if isempty(obj.motif),      obj.motif       =   '|';                    end
    if isempty(obj.showZero),   obj.showZero    =   true;                   end
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    % Parent-class components
    ClassName           =   'QMol_progressBar';
    PropNames           =   {'N','barLength','vMin','startValue','vMax','endValue','motif','showZero'};
end
end
%% Progress bar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function start(obj,v)
%start starts a new progress bar
%
%   obj.start() or obj.start(0) initiates an empty progress bar.
%
%   obj.start(v) initiates and displays the progress associated with the
%   value v along the bar.

    % Initialization
    if nargin == 1,     v           =   0;          end

    % Create the progress bar
    if obj.showZero,    obj.iBar    =   1;          fprintf('%s',obj.motif);
    else,               obj.iBar    =   0;          end

    % Display initial progress
    if v > 0,           obj.update(v);      end
end
function update(obj,v)
%update update the progress bar
%
%   obj.update(v) updates the bar to the progress associated with the value
%   v.

    ind                 =   min(round(obj.N*(v-obj.vMin)/(obj.vMax-obj.vMin)),obj.N);
    while ind > obj.iBar
        fprintf('%s',obj.motif);
        obj.iBar        =   obj.iBar + 1;  
    end
end
function finish(obj,newLine)
%finish finishes the progress bar
%
%   obj.finish() or obj.finish(true) finishes the progress bar by filling
%   it up completely and starts a new line.
%
%   obj.finish(false) finishes the progress bar without starting a new
%   line.

    % Fill up the bar
    while obj.iBar < obj.N
        fprintf('%s',obj.motif);
        obj.iBar        =   obj.iBar + 1;
    end

    % New line
    if nargin == 1   ||   newLine
        fprintf('\n');
    end

    % Delete progress bar counter
    obj.iBar            =   [];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

