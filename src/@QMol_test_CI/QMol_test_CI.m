classdef QMol_test_CI < QMol_test
%QMol_test_CI_conv suite of unit tests for QMol_CI_conv

%   Version     Date        Author
%   01.22.000   05/20/2025  F. Mauger
%       Creation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_test)
function version
    QMol_doc.showVersion('01.22.000','05/20/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test_CI\n'); 
    QMol_test_CI.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Object initialization
    obj.test_initialization;

    % Set configuration state basis
    obj.test_commonConfigurationBasis;
    obj.test_setConfigurationBasis_CIS;
    obj.test_setConfigurationBasis_CISD;
    obj.test_setConfigurationBasis_RAS;

    % Ground state and energy
    obj.test_getEnergy;
    obj.test_computeGroundState;
    
end
end
%% Test components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=private)
    function test_commonConfigurationBasis(obj)

    % Initialization
    obj.showSection('Common configuration state basis features');

    % cleanConfigurationBasis
    CI                  =   QMol_CI_conv('configurationBasis',[-1 -2 1 2; -1 -3 1 2; -2 -1 1 2;-1 -2 1 3]);
    CI.cleanConfigurationBasis;

    obj.showResult('cleanConfigurationBasis',all(CI.configurationBasis == [-1 -2 1 2; -1 -3 1 2; -1 -2 1 3]));
end
% Methods in separate files
    test_initialization(obj)
    test_setConfigurationBasis_CIS(obj)
    test_setConfigurationBasis_CISD(obj)
    test_setConfigurationBasis_RAS(obj)
    test_getEnergy(obj)
    test_computeGroundState(obj)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

