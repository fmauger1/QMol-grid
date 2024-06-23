classdef QMol_test_disc_basis < QMol_test
%QMol_test_disc_basis suite of unit tests for QMol_disc_basis

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
    fprintf('  * QMol_test_disc_basis\n'); 
    QMol_test_disc_basis.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Common stuff
    obj.test_initialization;
    obj.test_operator;
    % (TD)DFT stuff
    obj.test_DFT_classNames;
    obj.test_DFT_allocation;
    obj.test_DFT_orbital;
    obj.test_DFT_miscellaneous;
    obj.test_DFT_operator_spinRes;
    obj.test_DFT_operator_spinPol;
    % (TD)SE stuff
    obj.test_SE_classNames;
    obj.test_SE_allocation;
    obj.test_SE_waveFunction;
    obj.test_SE_miscellaneous;
    obj.test_SE_operator;
end
end
methods (Access=private)
    % Common stuff
    test_initialization(obj)
    test_operator(obj)
    % (TD)DFT stuff
    test_DFT_classNames(obj)
    test_DFT_allocation(obj)
    test_DFT_orbital(obj)
    test_DFT_miscellaneous(obj)
    test_DFT_operator_spinRes(obj)
    test_DFT_operator_spinPol(obj)
    % (TD)SE stuff
    test_SE_classNames(obj)
    test_SE_allocation(obj)
    test_SE_waveFunction(obj)
    test_SE_miscellaneous(obj)
    test_SE_operator(obj)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

