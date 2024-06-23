classdef QMol_test_DFT < QMol_test
%QMol_test_DFT suite of unit tests for QMol_DFT

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
    fprintf('  * QMol_test_DFT\n'); 
    QMol_test_DFT.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(~)
%testUnit run all unit tests on the class
    
    fprintf(['*** QMol_DFT is an abstract class defining the common interface for DFT\n' ...
             '    simulations in the QMol-grid package. It has no associated unit tests.\n' ...
             '    > QMol_DFT_spinPol implements spin-polarized models\n' ...
             '    > QMol_DFT_spinRes implements spin-restricted models\n' ...
             '    see their respective suite of unit tests.                         ****\n']);
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

