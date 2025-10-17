classdef QMol_test_TDCI_sympSplitOp < QMol_test
%QMol_test_TDCI_sympSplitOp suite of unit tests for QMol_TDCI_sympSplitOp

%   Version     Date        Author
%   01.23.000   06/04/2025  F. Mauger
%       Creation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_test)
function version
    QMol_doc.showVersion('01.23.000','06/04/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test_TDCI_sympSplitOp\n'); 
    QMol_test_TDCI_sympSplitOp.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(~)
%testUnit run all unit tests on the class
    
    fprintf(['*** QMol_TDCI_sympSplitOp is an abstract class of the kernel of the \n' ...
             '    QMol-grid package. It has no unit test associated to it.          ****\n']);
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

