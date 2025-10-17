classdef QMol_test_TDDFT_ESSO_2S < QMol_test
%QMol_test_TDDFT_ESSO_2S suite of unit tests for QMol_TDDFT_ESSO_2S

%   Version     Date        Author
%   01.23.000   07/22/2025  F. Mauger
%       Creation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_test)
function version
    QMol_doc.showVersion('01.23.000','07/22/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test_TDDFT_ESSO_2S\n'); 
    QMol_test_TDDFT_ESSO_2S.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(~)
%testUnit run all unit tests on the class

    fprintf(['*** QMol_test_TDDFT_ESSO_2S does not define any unit tests. Run the TDDFT\n' ...
             '    documentation examples to test the Strang-splitting propagator.  ****\n']);
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

