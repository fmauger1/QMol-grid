classdef QMol_test_TDDFT_ESSO_4BM < QMol_test
%QMol_test_TDDFT_ESSO_4BM suite of unit tests for QMol_TDDFT_ESSO_4BM

%   Version     Date        Author
%   01.23.000   07/23/2025  F. Mauger
%       Creation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_test)
function version
    QMol_doc.showVersion('01.23.000','07/23/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test_TDDFT_ESSO_4BM\n'); 
    QMol_test_TDDFT_ESSO_4BM.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(~)
%testUnit run all unit tests on the class

    fprintf(['*** QMol_test_TDDFT_ESSO_4BM does not define any unit tests. Run the\n' ...
             '    TDDFT documentation examples to test the Forest-Ruth propagator. ****\n']);
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

