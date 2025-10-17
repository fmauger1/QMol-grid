classdef QMol_test_TDCI_SSO_4BM < QMol_test
%QMol_test_TDCI_SSO_4BM suite of unit tests for QMol_TDCI_SSO_4BM

%   Version     Date        Author
%   01.23.000   06/07/2025  F. Mauger
%       Creation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_test)
function version
    QMol_doc.showVersion('01.23.000','06/07/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test_TDCI_SSO_4BM\n'); 
    QMol_test_TDCI_SSO_4BM.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(~)
%testUnit run all unit tests on the class

    fprintf(['*** QMol_test_TDCI_SSO_4BM does not define any unit tests. Run the TDCI-\n' ...
             '    documentation examples to test the Blanes and Moan 4th order\n' ...
             '    propagator.                                                       ****\n']);
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

