classdef QMol_test_SE_wfcn_basis < QMol_test
%QMol_test_SE_wfcn_basis suite of unit tests for QMol_SE_wfcn_basis

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
    fprintf('  * QMol_test_SE_wfcn_basis\n'); 
    QMol_test_SE_wfcn_basis.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test components
    obj.test_operator;
    
end
end
%% Test components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=private)
function test_operator(obj) %==============================================
%test_operator unit tests for operator overload in the class
    
    % Initialization
    obj.showSection('Object comparison');
    
    % Spin restricted
    R                   =   rand(15,7);
    V1                  =   QMol_SE_wfcn('wavefunction',R);
    V2                  =   QMol_SE_wfcn('wavefunction',R + .5e-10);
    obj.showResult('== (eq operator)' ,V1 == V2);

    V2.clear;       V2.set('wfcn',R + 1.01e-10);
    obj.showResult('~= (value mismatch)' ,V1 ~= V2);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

