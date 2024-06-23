classdef QMol_test_DFT_Vks_grad < QMol_test
%QMol_test_DFT_Vks_grad suite of unit tests for QMol_DFT_Vks_grad

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
    fprintf('  * QMol_test_DFT_Vks_grad\n'); 
    QMol_test_DFT_Vks_grad.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test units
    obj.runArithmetic_spinRes;
    obj.runArithmetic_spinPol;
    
end
end
%% Test components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=private)
function runArithmetic_spinRes(obj) %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%runTest unit tests for the class
    
    % Initialization
    obj.showSection('Arithmetic with potential gradients (spin restricted)');

    x                   =   QMol_disc('x',-10:.1:20);
    DFT                 =   QMol_DFT_spinRes('discretization',x);
    x.initialize(DFT);

    % Test add
    DVks                 =   x.DFT_allocatePotentialGradient;

    DV                  =   rand(numel(x.x),1);
    DVks.add(1,DV);
    R                   =   all(size(DVks.potentialGradient) == size(DV))   &&  ... same size
                            all(abs(DVks.DV - DV) < 1e-14);                       % same values
    obj.showResult('add (explicit, to empty potential)',R);

    DVks.add(1,DV);
    R                   =   all(size(DVks.potentialGradient) == size(DV))   &&  ... same size
                            all(abs(DVks.DV - 2*DV) < 1e-14);                     % same values
    obj.showResult('add (explicit, to existing potential)',R);

    % Test applyPotential
    p                   =   rand(numel(x.x),1);
    DVp                 =   DVks.applyPotentialGradient(1,p);
    R                   =   all(size(DVks.potentialGradient) == size(DVp))  &&  ... same size
                            all(abs(DVks.DV.*p - DVp) < 1e-14);                   % same values
    obj.showResult('applyPotentialGradient (explicit)',R);

    DVks.add(1,@(p) p);   DVks.add(1,@(p) p.^2);
    DVp                 =   DVks.applyPotentialGradient(1,p);
    R                   =   all(size(DVks.potentialGradient) == size(DVp))  &&  ... same size
                            all(abs(DVks.DV.*p+p+p.^2 - DVp) < 1e-14);            % same values
    obj.showResult('applyPotentialGradient (implicit)',R);

end
function runArithmetic_spinPol(obj) %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%runTest unit tests for the class
    
    % Initialization
    obj.showSection('Arithmetic with potential gradients (spin polarized)');

    x                   =   QMol_disc('x',-10:.1:20);
    DFT                 =   QMol_DFT_spinPol('discretization',x);
    x.initialize(DFT);

    % Test add
    DVks                =   x.DFT_allocatePotentialGradient;

    DV                  =   rand(numel(x.x),1);
    DVks.add(1,DV);
    R                   =   all(size(DVks.potentialGradientUp) == size(DV))     &&  ... same size
                            all(size(DVks.potentialGradientDown) == size(DV))   &&  ... 
                            all(abs(DVks.DVup - DV) < 1e-14)                    &&  ... same values
                            all(abs(DVks.DVdw - DV) < 1e-14);
    obj.showResult('add(1,DV) (explicit, to empty potential)',R);

    DVks.add(1,DV);
    R                   =   all(size(DVks.potentialGradientUp) == size(DV))     &&  ... same size
                            all(size(DVks.potentialGradientDown) == size(DV))   &&  ... 
                            all(abs(DVks.DVup - 2*DV) < 1e-14)                  &&  ... same values
                            all(abs(DVks.DVdw - 2*DV) < 1e-14);
    obj.showResult('add(1,DV) (explicit, to existing potential)',R);
    
    x.DFT_allocatePotentialGradient(DVks);
    DVks.add(1,DV,0*DV);
    R                   =   all(size(DVks.potentialGradientUp) == size(DV))     &&  ... same size
                            all(size(DVks.potentialGradientDown) == size(DV))   &&  ... 
                            all(abs(DVks.DVup - DV) < 1e-14)                    &&  ... same values
                            all(DVks.DVdw == 0);           % same values
    obj.showResult('add(1,V_up,V_down) (explicit, to empty potential)',R);

    DVks.add(1,DV,-DV);
    R                   =   all(size(DVks.potentialGradientUp) == size(DV))     &&  ... same size
                            all(size(DVks.potentialGradientDown) == size(DV))   &&  ... 
                            all(abs(DVks.DVup - 2*DV) < 1e-14)                  &&  ... same values
                            all(abs(DVks.DVdw + DV) < 1e-14);
    obj.showResult('add(V_up,V_down) (explicit, to existing potential)',R);

    % Test applyPotential
    p                   =   rand(numel(x.x),1);

    DVp                 =   DVks.applyPotentialGradient(1,p,true);
    R                   =   all(size(DVks.potentialGradientUp) == size(DVp))    &&  ... same size
                            all(abs(DVks.DVup.*p - DVp) < 1e-14);                     % same values
    obj.showResult('applyPotentialGradient (explicit, spin up)',R);
    
    DVp                 =   DVks.applyPotentialGradient(1,p,false);
    R                   =   all(size(DVks.potentialGradientDown) == size(DVp))  &&  ... same size
                            all(abs(DVks.DVdw.*p - DVp) < 1e-14);                     % same values
    obj.showResult('applyPotentialGradient (explicit, spin down)',R);

    DVks.add(1,@(p,isUp) p+isUp*p.^2);      DVks.add(1,@(p,isUp) p.^3+2*(~isUp)*p);
    DVp                 =   DVks.applyPotentialGradient(1,p,true);
    R                   =   all(size(DVks.potentialGradientUp) == size(DVp))    &&  ... same size
                            all(abs(DVks.DVup.*p+p+p.^2+p.^3 - DVp) < 1e-14);         % same values
    obj.showResult('applyPotentialGradient (implicit, spin up)',R);
    
    DVp                 =   DVks.applyPotentialGradient(1,p,false);
    R                   =   all(size(DVks.potentialGradientDown) == size(DVp))  &&  ... same size
                            all(abs(DVks.DVdw.*p+3*p+p.^3 - DVp) < 1e-14);           % same values
    obj.showResult('applyPotentialGradient (implicit, spin down)',R);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

