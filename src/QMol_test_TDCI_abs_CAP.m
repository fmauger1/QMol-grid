classdef QMol_test_TDCI_abs_CAP < QMol_test
%QMol_test_TDCI_abs_CAP of unit tests for QMol_TDCI_abs_CAP

%   Version     Date        Author
%   01.23.000   06/04/2025  F. Mauger
%       Creation
%   01.23.001   06/05/2025  F. Mauger
%       Add test for getDampingMatrix

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_test)
function version
    QMol_doc.showVersion('01.23.001','06/05/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test_TDCI_abs_CAP\n'); 
    QMol_test_TDCI_abs_CAP.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class

    % Initialization
    obj.showSection('CAP-type damping matrix');

    % Test damping matrix
    randStr             =   RandStream('dsfmt19937','Seed',0);
    M                   =   rand(randStr,10,10);

    ABS                 =   QMol_TDCI_abs_CAP('dampingMatrix',M);
    ABS.initialize();

    obj.showResult('getDampingMatrix',all(M == ABS.getDampingMatrix,'all'));

    % Test isCAP flag
    obj.showResult('isCAP',ABS.isCAP);
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

