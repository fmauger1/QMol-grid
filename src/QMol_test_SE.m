classdef QMol_test_SE < QMol_test
%QMol_test_SE suite of unit tests for QMol_SE

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
    fprintf('  * QMol_test_SE\n'); 
    QMol_test_SE.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test units
    obj.test_initialization;
    obj.test_SE_energy;
    
end
end
methods (Access=private)
function test_initialization(obj) %========================================
    
    % DFT model
    obj.showSection('Initialization');
    
    X                   =   -20:.1:15;

    Vse                 =   QMol_SE_V('inputPotential',@(x)-5./sqrt(x.^2+.7));

    SE                  =   QMol_SE( ...
                                'xspan',                        X,          ...
                                'numberWaveFunction',           3,          ...
                                'potential',                    Vse);
    
    SE.initialize;

    % discretization 
    obj.showResult('discretization',SE.discretization == QMol_disc('xspan',X));

    % wave functions
    obj.showResult('waveFunction', ...
        all(size(SE.waveFunction.waveFunction) == [numel(X) 3]));           % Wave function is allocated

    % isInitialized
    obj.showResult('isInitialized',SE.isInitialized);

    % dim
    obj.showResult('dim',SE.dim == 1);
    
end
function test_SE_energy(obj) %=============================================

    % DFT model
    obj.showSection('Energy components');
    
    X                   =   -20:.1:15;

    Vse                 =   QMol_SE_V('inputPotential',@(x)-5./sqrt(x.^2+.7));

    SE                  =   QMol_SE( ...
                                'xspan',                        X,  ...
                                'numberWaveFunction',           3,  ...
                                'potential',                    Vse);
    
    SE.initialize;

    dx                  =   X(2)-X(1);      X = X.';
    p                   =  [      exp(-X.^2*.5    ) / sqrt(sum( (      exp(-X.^2*.5    )).^2 )*dx), ...
                            X   .*exp(-X.^2*.5/2^2) / sqrt(sum( (X   .*exp(-X.^2*.5/2^2)).^2 )*dx), ...
                            X.^2.*exp(-X.^2*.5/3^2) / sqrt(sum( (X.^2.*exp(-X.^2*.5/3^2)).^2 )*dx)];
    SE.wfcn.wfcn        =   p;

    % SE energy
    [Etot,Ekin,Epot]    =   SE.getEnergy('SE');

    obj.showResult('getEnergy (SE -- total)',Etot == Ekin+Epot);

    E                   =   SE.disc.SE_energyKinetic(SE.wfcn);
    obj.showResult('getEnergy (SE -- kinetic)', E == Ekin);

    E                   =   Vse.getEnergy(SE.wfcn);
    obj.showResult('getEnergy (SE -- external)', E == Epot);

    % Wave function energy
    [E_1,DE_1]          =   SE.getEnergy('wfcn');
    [E_2,DE_2]          =   SE.disc.SE_energyWaveFunction(Vse,SE.wfcn);
    obj.showResult('getEnergy (wave function)', ...
        all(E_1 == E_2)   &&  all(DE_1 == DE_2));

    V                   =   QMol_SE_V('inputPotential',1.2*Vse.V);
    V.initialize(SE);

    [E_1,DE_1]          =   SE.getEnergy('wfcn',V);
    [E_2,DE_2]          =   SE.disc.SE_energyWaveFunction(V,SE.wfcn);
    obj.showResult('getEnergy (wave function -- potential supplied)', ...
        all(E_1 == E_2)   &&  all(DE_1 == DE_2));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

