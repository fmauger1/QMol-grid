classdef QMol_test_TDCI_abs_mask < QMol_test
%QMol_test_TDCI_abs_mask of unit tests for QMol_TDCI_abs_mask

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
    fprintf('  * QMol_test_TDCI_abs_mask\n'); 
    QMol_test_TDCI_abs_mask.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class

    % Initialization
    obj.showSection('mask-type damping matrix');

    % Test damping matrix
    randStr             =   RandStream('dsfmt19937','Seed',0);
    M                   =   rand(randStr,10,10);

    ABS                 =   QMol_TDCI_abs_mask('dampingMatrix',M);
    ABS.initialize();

    obj.showResult('getDampingMatrix',all(M == ABS.getDampingMatrix,'all'));

    % Test isCAP flag
    obj.showResult('isCAP',~ABS.isCAP);
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

