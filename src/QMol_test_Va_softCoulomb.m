classdef QMol_test_Va_softCoulomb < QMol_test
%QMol_test_Va_softCoulomb suite of unit tests for QMol_Va_softCoulomb

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
    fprintf('  * QMol_test_Va_softCoulomb\n'); 
    QMol_test_Va_softCoulomb.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Initialization

    Z                   =   15;
    a                   =   pi;
    X0                  =   log(15);
    V                   =   @(x) -Z ./ sqrt((x-X0).^2 + a^2 );

    Va                  =   QMol_Va_softCoulomb('charge',Z,'softeningParameter',a,'position',X0);
    Va.initialize;
    
    x                   =  -10:.08:16;
    obj.showResult('getPotential' ,max(abs( Va.getPotential(x)-V(x) )) < 1e-10);

    dx                  =   1e-5;
    obj.showResult('getPotentialDerivative' ,max(abs( Va.getPotentialDerivative(1,x)-(V(x+.5*dx)-V(x-.5*dx))/dx )) < 5e-10);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

