classdef QMol_test_Va_Gaussian < QMol_test
%QMol_test_Va_Gaussian suite of unit tests for QMol_Va_Gaussian

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
    fprintf('  * QMol_test_Va_Gaussian\n'); 
    QMol_test_Va_Gaussian.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Initialization
    obj.showSection('Potential');

    V0                  =   pi;
    s                   =   5/3;
    X0                  =   exp(1);
    V                   =   @(x) -V0 * exp(-(x-X0).^2 / (2 * s^2) );

    Va                  =   QMol_Va_Gaussian('potentialDepth',V0,'potentialWidth',s,'position',X0);
    Va.initialize;
    
    x                   =  -15:.1:20;
    obj.showResult('getPotential' ,max(abs( Va.getPotential(x)-V(x) )) < 1e-10);

    dx                  =   1e-5;
    obj.showResult('getPotentialDerivative' ,max(abs( Va.getPotentialDerivative(1,x)-(V(x+.5*dx)-V(x-.5*dx))/dx )) < 5e-10);
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

