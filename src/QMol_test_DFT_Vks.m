classdef QMol_test_DFT_Vks < QMol_test
%QMol_test_DFT_Vks suite of unit tests for QMol_DFT_Vks

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
    fprintf('  * QMol_test_DFT_Vks\n'); 
    QMol_test_DFT_Vks.version;
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
    obj.showSection('Arithmetic with potentials (spin restricted)');

    x                   =   QMol_disc('x',-10:.1:20);
    DFT                 =   QMol_DFT_spinRes('discretization',x);
    x.initialize(DFT);

    % Test add
    Vks                 =   x.DFT_allocatePotential;

    V                   =   rand(numel(x.x),1);
    Vks.add(V);
    R                   =   all(size(Vks.potential) == size(V))             &&  ... same size
                            all(abs(Vks.V - V) < 1e-14);                          % same values
    obj.showResult('add (explicit, to empty potential)',R);

    Vks.add(V);
    R                   =   all(size(Vks.potential) == size(V))             &&  ... same size
                            all(abs(Vks.V - 2*V) < 1e-14);                        % same values
    obj.showResult('add (explicit, to existing potential)',R);

    % Test applyPotential
    p                   =   rand(numel(x.x),1);
    Vp                  =   Vks.applyPotential(p);
    R                   =   all(size(Vks.potential) == size(Vp))            &&  ... same size
                            all(abs(Vks.V.*p - Vp) < 1e-14);                      % same values
    obj.showResult('applyPotential (explicit)',R);

    Vks.add(@(p) p);    Vks.add(@(p) p.^2);
    Vp                  =   Vks.applyPotential(p);
    R                   =   all(size(Vks.potential) == size(Vp))            &&  ... same size
                            all(abs(Vks.V.*p+p+p.^2 - Vp) < 1e-14);               % same values
    obj.showResult('applyPotential (implicit)',R);

end
function runArithmetic_spinPol(obj) %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%runTest unit tests for the class
    
    % Initialization
    obj.showSection('Arithmetic with potentials (spin polarized)');

    x                   =   QMol_disc('x',-10:.1:20);
    DFT                 =   QMol_DFT_spinPol('discretization',x);
    x.initialize(DFT);

    % Test add
    Vks                 =   x.DFT_allocatePotential;

    V                   =   rand(numel(x.x),1);
    Vks.add(V);
    R                   =   all(size(Vks.potentialUp) == size(V))           &&  ... same size
                            all(size(Vks.potentialDown) == size(V))         &&  ... 
                            all(abs(Vks.Vup - V) < 1e-14)                   &&  ... same values
                            all(abs(Vks.Vdw - V) < 1e-14);
    obj.showResult('add(V) (explicit, to empty potential)',R);

    Vks.add(V);
    R                   =   all(size(Vks.potentialUp) == size(V))           &&  ... same size
                            all(size(Vks.potentialDown) == size(V))         &&  ... 
                            all(abs(Vks.Vup - 2*V) < 1e-14)                 &&  ... same values
                            all(abs(Vks.Vdw - 2*V) < 1e-14);
    obj.showResult('add(V) (explicit, to existing potential)',R);
    
    x.DFT_allocatePotential(Vks);
    Vks.add(V,0*V);
    R                   =   all(size(Vks.potentialUp) == size(V))           &&  ... same size
                            all(size(Vks.potentialDown) == size(V))         &&  ... 
                            all(abs(Vks.Vup - V) < 1e-14)                   &&  ... same values
                            all(Vks.Vdw == 0);           % same values
    obj.showResult('add(V_up,V_down) (explicit, to empty potential)',R);

    Vks.add(V,-V);
    R                   =   all(size(Vks.potentialUp) == size(V))           &&  ... same size
                            all(size(Vks.potentialDown) == size(V))         &&  ... 
                            all(abs(Vks.Vup - 2*V) < 1e-14)                 &&  ... same values
                            all(abs(Vks.Vdw + V) < 1e-14);
    obj.showResult('add(V_up,V_down) (explicit, to existing potential)',R);

    % Test applyPotential
    p                   =   rand(numel(x.x),1);

    Vp                  =   Vks.applyPotential(p,true);
    R                   =   all(size(Vks.potentialUp) == size(Vp))          &&  ... same size
                            all(abs(Vks.Vup.*p - Vp) < 1e-14);                    % same values
    obj.showResult('applyPotential (explicit, spin up)',R);
    
    Vp                  =   Vks.applyPotential(p,false);
    R                   =   all(size(Vks.potentialDown) == size(Vp))        &&  ... same size
                            all(abs(Vks.Vdw.*p - Vp) < 1e-14);                    % same values
    obj.showResult('applyPotential (explicit, spin down)',R);

    Vks.add(@(p,isUp) p+isUp*p.^2);     Vks.add(@(p,isUp) p.^3+2*(~isUp)*p);
    Vp                  =   Vks.applyPotential(p,true);
    R                   =   all(size(Vks.potentialUp) == size(Vp))          &&  ... same size
                            all(abs(Vks.Vup.*p+p+p.^2+p.^3 - Vp) < 1e-14);        % same values
    obj.showResult('applyPotential (implicit, spin up)',R);
    
    Vp                  =   Vks.applyPotential(p,false);
    R                   =   all(size(Vks.potentialDown) == size(Vp))        &&  ... same size
                            all(abs(Vks.Vdw.*p+3*p+p.^3 - Vp) < 1e-14);           % same values
    obj.showResult('applyPotential (implicit, spin down)',R);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

