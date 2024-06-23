classdef QMol_test_disc < QMol_test
%QMol_test_disc suite of unit tests for QMol_disc

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_test)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test_disc\n'); 
    QMol_test_disc.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Object initialization
    obj.test_initialization;
    obj.test_operator;
    % DFT methods
    obj.test_DFT_classNames;
    obj.test_DFT_allocation;
    obj.test_DFT_miscellaneous;
    obj.test_DFT_operator_spinRes;
    obj.test_DFT_operator_spinPol;
    % Schrodinger-equation methods
    obj.test_SE_classNames;
    obj.test_SE_allocation;
    obj.test_SE_miscellaneous;
    obj.test_SE_operator;
end
end
methods (Access=private)
    % Object initializetion
    test_initialization(obj)
    test_operator(obj)
    % DFT methods
    test_DFT_classNames(obj)
    test_DFT_allocation(obj)
    test_DFT_miscellaneous(obj)
    test_DFT_operator_spinRes(obj)
    test_DFT_operator_spinPol(obj)
    % Schrodinger-equation methods
    test_SE_classNames(obj)
    test_SE_allocation(obj)
    test_SE_miscellaneous(obj)
    test_SE_operator(obj)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

