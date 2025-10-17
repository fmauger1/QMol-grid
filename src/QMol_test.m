classdef (Abstract) QMol_test < QMol_suite
%QMol_test test suite for the QMol-grid package
%   By convention, all classes in the QMol-grid package implement a
%   corresponding unit test, even if no actual test is performed.
%
%   QMol_test.test() runs the entire test suite (depending on the system,
%   this may take a while). The test suite prints out the results of each
%   test, together with a summary of the number of tests that passed.
%
%   QMol_test.test('test1','test2',___) selectively runs the set of named
%   tests. The input test names correspond to the QMol-grid classes to be 
%   tested, omitting the initial 'QMol_' part of the name. For instance, 
%   'test1' runs the unit tests for the QMol_test1 component. Names of 
%   tests are case insensitive.
%
%   Activate the summary mode, where only a summary of the test results is
%   displayed, by including '-summary' as the first argument
%       QMol_test.test('-summary')                          all tests
%       QMol_test.test('-summary','test1','test2',___)      selected tests
%
%   Editable properties: N/A
%
%   Methods: test
%
%   See also QMol_suite, QMol_doc 
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release
%   01.21.001   12/08/2024  F. Mauger
%       Clean docstring

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.001','12/08/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test:\n      > Test suite interface\n'); QMol_test.version;
end
end
methods (Static,Access=?QMol_doc)
function showComponents
%showComponents displays the list of available unit tests in the QMol-grid
%   package. 

%   NOTE: 
%   * Test components associated with the selected discretization framework
%     are listed in the QMol_info class.
%     Warning: this will probably change when support for 2 (and 3)
%     dimensional calculations is included.
    
    % Get list of unit tests
    l                   =   QMol_testInfo.testList;

    % Parse list of unit tests
    QMol_doc.showSection('Test suite');

    for k = 1:numel(l)
        eval(['QMol_test_' l{k} '.showInfo']);
    end
    
    fprintf('\n')

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Access=private)
    isSummary           =   false
    nbError             =   0
end
%% Display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_suite)
function showSection(S)
%showSection formats and diplays a section break
    fprintf('  * %s\n',S);
end
end
methods (Access=?QMol_test)
function showResult(obj,S,R)
%showResult formats and displays the result of a test
    
    % Test result
    if R
        msg             =   ' PASS';
    else
        obj.nbError     =   obj.nbError + 1;
        msg             =   '*FAIL';
    end

    % Display result
    if ~obj.isSummary   ||   any(~R,'all')
        fprintf('    > %-62s %+5s\n',S,msg);
    end
end
end
methods (Access=public)
function obj = QMol_test(S)
%QMol_test constructor, to handle summary test runs
    
    % Initialization
    if nargin < 1,  S   =   false;  end

    % Update summary mode
    obj.isSummary       =   S;
end
end
%% Run unit tests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=public)
function test(varargin)
%test runs the unit tests for the elected input components
%   QMol_test.test()                        run all tests
%   QMol_test.test('test1','test2',___)     run selected tests
%   QMol_test.test('-summary'___)           run in summary mode

    % List of tests
    LU                  =   QMol_testInfo.testList;
    if nargin == 0                                              % Full test suite, no summary
        TU              =   LU;
        isSummary       =   false;                                          %#ok<NASGU> 
    elseif any(strcmpi(varargin{1},{'-summary','-short','-'}))  % Summarized test
        isSummary       =   true;                                           %#ok<NASGU> 
        if nargin == 1
            TU          =   LU;
        else
            TU          =   varargin(2:end);
        end
    else                                                        % Selected tests, no summary
        isSummary       =   false;                                          %#ok<NASGU> 
        TU              =   varargin;
    end

    % Test suite header
    fprintf(['########################### QMol-grid toolbox ############################\n'...
             '  * Test suite (release)\n']);
    QMol_test.version;
    QMol_testInfo.showHeader;
    QMol_doc.showFooter; fprintf('\n');

    QMol_doc.showSection('External components');
    fprintf('  * convertUnit\n'); QMol_doc.showVersion(convertUnit.version,convertUnit.lastMod,convertUnit.author);
    fprintf('  * fourierTool\n'); QMol_doc.showVersion(convertUnit.version,convertUnit.lastMod,convertUnit.author);
    fprintf('\n')

    QMol_doc.showSection('License');
    QMol_doc.showLicense;
    fprintf('\n')

    QMol_doc.showFunding;
    
    % Run individual unit tests
    FT                  =   {};

    for k = 1:numel(TU)
        % Identify unit test
        l               =   find(strcmpi(TU{k},LU));
        if isempty(l)
            warning('QMolGrid:test:testName',['Unknown unit test ' TU{k} ' -- entry ignored']);
            continue
        end

        QMol_doc.showSection(LU{l});

        % Run unit test
        T               =   eval(['QMol_test_' LU{l} '(isSummary)']);
        T.testUnit;

        fprintf('    ------------------------\n')
        eval(['QMol_test_' LU{l} '.version;'])

        % Finalize unit test
        if T.nbError > 0
            fprintf('    %+70s\n\n',[num2str(T.nbError) ' test(s) FAILED'])
            FT = [FT, LU(l)];                                               %#ok<AGROW> 
        else
            fprintf('    %+70s\n\n','All test(s) PASSED')
        end

        clear T
    end

    % Summary
    QMol_doc.showSection('Summary')
    
    TU                  =   numel(TU);

    if isempty(FT)
        fprintf('  * PASSED all (%u) test(s)\n\n',TU)
    else
        fprintf('  * PASSED %u and FAILED %u tests\n',TU-numel(FT),numel(FT));
        fprintf('    Failed: %s',FT{1});
        if numel(FT) > 1, fprintf(', %s',FT{2:end}); end
        fprintf('\n\n')
    end

    % Test suite final footer
    QMol_doc.showFooter;
end
end
%% DFT methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods(Static,Hidden,Access=protected)
function Ex = DFT_1D_LDA_eps_x(funV,rho)
%DFT_1D_LDA_eps_x(funV,rho) computes the 1D local-density approximation
%   (LDA) exchange-energy per particle. funV is a handle function
%   describing the electron-electron interaction potential and rho is an
%   array containing the values of the one-body density at which the
%   exchange-energy per particle should be computed.
    
    % Initialization
    Ex                  =   NaN(size(rho));

    % Compute energy
    for k = 1:numel(rho)
        a               =   2/(pi*rho(k));
        Ex(k)           =   integral(@(u) -sin(u).^2./u.^2.*funV(u*a)/pi,0,Inf);
    end
end
function Ex = DFT_1D_LDA_d_eps_x(funDV,rho)
%DFT_1D_LDA_d_eps_x(funDV,rho) computes the derivative of the 1D local-
%   density approximation (LDA) exchange-energy per particle. funDV is a
%   handle function describing the derivative of the electron-electron
%   interaction potential and rho is an array containing the values of the
%   one-body density at which the exchange-energy per particle should be
%   computed.
    
    % Initialization
    Ex                  =   NaN(size(rho));

    % Compute energy
    for k = 1:numel(rho)
        a               =   2/(pi*rho(k))  ;
        b               =   2/(pi*rho(k))^2;
        Ex(k)           =   integral(@(u)  sin(u).^2./u.*funDV(u*a)*b,0,Inf);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

