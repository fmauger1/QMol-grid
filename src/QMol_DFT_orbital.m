classdef QMol_DFT_orbital < QMol_suite
%QMol_DFT_orbital density-functional-theory (DFT) orbital(s) on
%   Cartesian-grid domain discretization.
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
methods (Static,Access={?QMol_doc,?QMol_DFT_orbital})
function showInfo
    fprintf('  * QMol_DFT_orbital:\n      > Kohn-Sham orbital\n'); 
    QMol_DFT_orbital.version;
end
end
methods (Access=public)
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components (do not trust KSO, KSOup, or KSOdw, they might
    % have been initialized empty; fetch size from domain and parent DFT)
    mem                 =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.disc.x),'imag');

    if obj.isSpinPol,                                                                               if opt
        % Spin polarized
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham orbitals', 0,1);
        QMol_DFT_profiler.showMemoryFootprint('spin up'     ,numel(obj.disc.QM.occ{1})*mem,2);
        QMol_DFT_profiler.showMemoryFootprint('spin down'   ,numel(obj.disc.QM.occ{2})*mem,2);      end
        
        mem             =   (numel(obj.disc.QM.occ{1}) + numel(obj.disc.QM.occ{2})) * mem;
    else,                                                                                           if opt
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham orbitals',numel(obj.disc.QM.occ)*mem,1);   end
        
        mem             =   numel(obj.disc.QM.occ)*mem;
    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Transient,Hidden,GetAccess=public,SetAccess=?QMol_suite)
    % Linked objects
    disc                            % discretization object (for possible future development)
end
properties (GetAccess=public,SetAccess=?QMol_suite)
    isSpinPol
end
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    KSO
    KSOup
    KSOdw
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    orbital                         % KSO
    orbitalUp                       % KSOup
    orbitalDown                     % KSOdw
end
methods (Access=public,Static)
    function b = isBasis,   b   =   false;  end
end
properties (Constant,Access=private)
    eqTol               =   1e-10   % tolerance for equality of data
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % KSO ~~~~~~~~~~~~~~
    function val = get.orbital(obj),            val         =   obj.KSO;    end
    function set.orbital(obj,val),              obj.KSO     =   val;        end
    % KSOup ~~~~~~~~~~~~
    function val = get.orbitalUp(obj),          val         =   obj.KSOup;  end
    function set.orbitalUp(obj,val),            obj.KSOup   =   val;        end
    % KSOdw ~~~~~~~~~~~~
    function val = get.orbitalDown(obj),        val         =   obj.KSOdw;  end
    function set.orbitalDown(obj,val),          obj.KSOdw   =   val;        end
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
        obj.isSpinPol   =   [];
        obj.KSO         =   [];
        obj.KSOup       =   [];
        obj.KSOdw       =   [];

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
    
    ClassName           =   'QMol_DFT_orbital';
    PropNames           =  {'isSpinPol','KSO','KSOup','KSOdw',...
                            'orbital','orbitalUp','orbitalDown'};
end
end
%% Operator overload %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function t = eq(data1,data2)
%eq overloads the == operator  

    % Test domain
    if      isempty(data1.disc)   ||   isempty(data2.disc),     t   =   true;            % not enough information on discretizations
    else,   t   =   data1.disc == data2.disc;                               end

    % Test actual data
    if data1.isSpinPol
        t               =   t   &&   data2.isSpinPol                                && ... both spin polarized
                            all(size(data1.KSOup) == size(data2.KSOup))             && ... same number of data
                            all(size(data1.KSOdw) == size(data2.KSOdw))             && ...
                            all(abs(data1.KSOup-data2.KSOup) <= data1.eqTol,'all')  && ... all elements are the same (within tolerance)
                            all(abs(data1.KSOdw-data2.KSOdw) <= data1.eqTol,'all');
    else
        t               =   t   &&   ~data2.isSpinPol                               && ... both spin restricted
                            all(size(data1.KSO) == size(data2.KSO))                 && ... same number of data
                            all(abs(data1.KSO-data2.KSO) <= data1.eqTol,'all');          % all elements are the same (within tolerance)
    end
end
function t = ne(data1,data2)
%ne overloads the ~= operator
    t                   =   ~(data1.eq(data2));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

