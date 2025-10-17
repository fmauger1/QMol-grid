classdef QMol_test_TDDFT_extSymp < QMol_test
%QMol_test_TDDFT_extSymp suite of unit tests for QMol_TDDFT_extSymp

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
    fprintf('  * QMol_test_TDDFT_extSymp\n'); 
    QMol_test_TDDFT_extSymp.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(~)
%testUnit run all unit tests on the class
    
    fprintf(['*** QMol_TDDFT_extSymp is an abstract class of the kernel of the \n' ...
             '    QMol-grid package. It has no unit test associated to it.          ****\n']);
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

