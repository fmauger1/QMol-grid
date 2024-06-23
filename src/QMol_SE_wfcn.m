classdef QMol_SE_wfcn < QMol_suite
%QMol_SE_wfcn Schrodinger-equation (SE) wave function(s) on Cartesian-grid 
%   domain discretization.
%
%   NOTES:
%   > One can change the data-equality criterion through the member
%     property eqTol. Side effects to doing so are untested, though (the
%     code is developed with eqTol = 1e-10).
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_SE_wfcn})
function showInfo
    fprintf('  * QMol_SE_wfcn:\n      > SE wave function(s)\n'); 
    QMol_SE_wfcn.version;
end
end
methods (Access=public)
function mem = getMemoryProfile(obj,opt,N)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 3,  N   =   1;
    if nargin < 2,  opt =   false;  end, end

    % Display components (do not trust wfcn, it might have been initialized 
    % empty; fetch size from domain and input number of wave functions)
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.disc.x),'imag'); if opt
    
    QMol_SE_profiler.showMemoryFootprint('Wave function(s)',N*mem,1);                       end
    mem                 =   N*mem;
    
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Transient,Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Linked objects
    disc                            % discretization object (for possible future development)
end
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    wfcn
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    waveFunction                    % wfcn
end
methods (Access=public,Static)
    function b = isBasis,   b   =   false;  end
end
properties (Constant,Access=private)
    eqTol               =   1e-10   % tolerance for equality of data
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % wfcn ~~~~~~~~~~~~~
    function val = get.waveFunction(obj),       val         =   obj.wfcn;    end
    function set.waveFunction(obj,val),         obj.wfcn    =   val;        end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function clear(obj,varargin)
%clear clears all or selected member properties
    
    % streamlined clear all
    if nargin == 1
        % Clear all properties
        obj.reset;

        obj.disc        =   [];
        obj.wfcn        =   [];

    % parent-class selected clear
    else
        clear@QMol_suite(obj,varargin{:});
    end
end
function initialize(obj,disc)
%initialize initializes the object
    
    % Reinitialize any component
    obj.reset;

    % Update discretization
    if nargin == 1,     obj.disc    =   [];
    else,               obj.disc    =   disc;                               end
    
    % Update initialization status
    obj.isInit          =   true;
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_SE_wfcn';
    PropNames           =  {'wfcn','waveFunction'};
end
end
%% Operator overload %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function t = eq(data1,data2)
%eq overloads the == operator  

    % Test domain
    if      isempty(data1.disc)   ||   isempty(data2.disc),     t   =   true;   % not enough information on discretizations
    else,   t   =   data1.disc == data2.disc;                               end

    % Test actual data
    t           =   t   &&  all(size(data1.wfcn) == size(data2.wfcn))       && ... same number of data
                            all(abs(data1.wfcn-data2.wfcn) <= data1.eqTol,'all');% all elements are the same (within tolerance)

end
function t = ne(data1,data2)
%ne overloads the ~= operator
    t                   =   ~(data1.eq(data2));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

