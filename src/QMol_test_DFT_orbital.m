classdef QMol_test_DFT_orbital < QMol_test
%QMol_test_DFT_orbital suite of unit tests for QMol_DFT_orbital

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
    fprintf('  * QMol_test_DFT_orbital\n'); 
    QMol_test_DFT_orbital.version;
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
    V1                  =   QMol_DFT_orbital('isSpinPol',false,'orbital',R);
    V2                  =   QMol_DFT_orbital('isSpinPol',false,'orbital',R + .5e-10);
    obj.showResult('== (eq operator, spin restricted)' ,V1 == V2);

    V2.clear;       V2.set('isSpinPol',false,'orbital',R + 1.01e-10);
    obj.showResult('~= (value mismatch, spin restricted)' ,V1 ~= V2);

    V2.clear;       V2.set('isSpinPol',false,'orbital',[]);
    obj.showResult('~= (size mismatch, spin restricted)' ,V1 ~= V2);

    % Spin polarized
    V1.clear;       V1.set('isSpinPol',true,'orbitalUp',2*R       ,'orbitalDown',R-1);
    V2.clear;       V2.set('isSpinPol',true,'orbitalUp',2*R-.5e-10,'orbitalDown',R-1+.5e-10);
    obj.showResult('== (eq operator, spin polarized)' ,V1 == V2);

    V2.clear;       V2.set('isSpinPol',true,'orbitalUp',2*R-1.01e-10,'orbitalDown',R-1+.5e-10);
    obj.showResult('~= (up-spin value mismatch, spin polarized)' ,V1 ~= V2);

    V2.clear;       V2.set('isSpinPol',true,'orbitalUp',2*R-.5e-10,'orbitalDown',R-1+1.01e-10);
    obj.showResult('~= (down-spin value mismatch, spin polarized)' ,V1 ~= V2);

    V2.clear;       V2.set('isSpinPol',true,'orbitalUp',[],'orbitalDown',R-1+.5e-10);
    obj.showResult('~= (up-spin size mismatch, spin polarized)' ,V1 ~= V2);

    V2.clear;       V2.set('isSpinPol',true,'orbitalUp',2*R-.5e-10,'orbitalDown',[]);
    obj.showResult('~= (down-spin size mismatch, spin polarized)' ,V1 ~= V2);


    V1.clear;       V1.set('isSpinPol',true,'orbital',R,'orbitalUp',2*R       ,'orbitalDown',R-1);
    V2.clear;       V2.set('isSpinPol',false,'orbital',R,'orbitalUp',2*R-.5e-10,'orbitalDown',R-1+.5e-10);
    obj.showResult('~= (spin polarized/restricted mismatch)' ,V1 ~= V2);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

