classdef QMol_test_TDCI_SSO < QMol_test
%QMol_test_TDCI_SSO suite of unit tests for QMol_TDCI_SSO

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
    fprintf('  * QMol_test_TDCI_SSO\n'); 
    QMol_test_TDCI_SSO.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(~)
%testUnit run all unit tests on the class
    
    fprintf(['*** QMol_test_TDCI_SSO is an abstract class of the QMol-grid package. It\n' ...
             '    has no unit test associated to it.                                ****\n']);
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

