classdef QMol_test_DFT_Vks_basis < QMol_test
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
    fprintf('  * QMol_test_DFT_Vks_basis\n'); 
    QMol_test_DFT_Vks_basis.version;
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

    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];

    x                   =   QMol_disc_basis('xspan',x,'basis',v);           x.orthonormalizeBasis;
    DFT                 =   QMol_DFT_spinRes('discretization',x);
    x.initialize(DFT);

    % Test add
    Vks                 =   x.DFT_allocatePotential;

    V                   =   rand(randStr,[numel(x.x),1]);
    Vks.add(V);
    R                   =   all(size(Vks.potential) == size(V))             &&  ... same size
                            all(abs(Vks.V - V) < 1e-14);                          % same values
    obj.showResult('add (explicit, to empty potential)',R);

    Vks.add(V);
    R                   =   all(size(Vks.potential) == size(V))             &&  ... same size
                            all(abs(Vks.V - 2*V) < 1e-14);                        % same values
    obj.showResult('add (explicit, to existing potential)',R);

    % Test applyPotential
    Vks.initialize(x);
    p                   =   rand(randStr,[x.basisSize,1]);
    Vp                  =   Vks.applyPotential(p);
    P                   =   zeros(size(Vks.potential));
    for k = 1:x.nV, P   =   P + p(k) * x.basis(:,k);                        end
    for k = 1:x.nV, Vp(k)=  Vp(k) - sum(x.basis(:,k).*Vks.potential.*P)*x.dx;   end
    obj.showResult('applyPotential (real basis)',all(abs(Vp) < 1e-10,'all'));

    v                   =   v.*[exp( 2i*x.x),exp( 3i*x.x),exp(-4i*x.x),exp(-2i*x.x)];
    x.set('basis',v);   x.orthonormalizeBasis;  x.initialize(DFT);
    Vks.add(zeros(size(Vks.potential)));        Vks.initialize(x);
    p                   =   rand(randStr,[x.basisSize,1]) + 1i*rand(randStr,[x.basisSize,1]);
    Vp                  =   Vks.applyPotential(p);
    P                   =   zeros(size(Vks.potential));
    for k = 1:x.nV, P   =   P + p(k) * x.basis(:,k);                        end
    for k = 1:x.nV, Vp(k)=  Vp(k) - sum(conj(x.basis(:,k)).*Vks.potential.*P)*x.dx;   end
    obj.showResult('applyPotential (complex basis)',all(abs(Vp) < 1e-10,'all'));

    Vks.add(@(p) p);    Vks.add(@(p) p.^2);     Vks.initialize(x);
    p                   =   rand(randStr,[x.basisSize,1]) + 1i*rand(randStr,[x.basisSize,1]);
    Vp                  =   Vks.applyPotential(p);
    P                   =   zeros(size(Vks.potential));
    for k = 1:x.nV, P   =   P + p(k) * x.basis(:,k);                        end
    for k = 1:x.nV, Vp(k)=  Vp(k) - sum(conj(x.basis(:,k)).*(Vks.potential.*P+P+P.^2))*x.dx;   end
    obj.showResult('applyPotential (implicit)',all(abs(Vp) < 1e-10,'all'));

end
function runArithmetic_spinPol(obj) %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%runTest unit tests for the class
    
    % Initialization
    obj.showSection('Arithmetic with potentials (spin polarized)');

    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];

    x                   =   QMol_disc_basis('x',x,'basis',v);               x.orthonormalizeBasis;
    DFT                 =   QMol_DFT_spinPol('discretization',x);
    x.initialize(DFT);

    % Test add
    Vks                 =   x.DFT_allocatePotential;

    V                   =   rand(randStr,[numel(x.x),1]);
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
    Vks.initialize(x);
    p                   =   rand(randStr,[x.basisSize,1]);
    Vp                  =   Vks.applyPotential(p,true);
    P                   =   zeros(size(Vks.potentialUp));
    for k = 1:x.nV, P   =   P + p(k) * x.basis(:,k);                        end
    for k = 1:x.nV, Vp(k)=  Vp(k) - sum(x.basis(:,k).*Vks.potentialUp.*P)*x.dx;   end
    obj.showResult('applyPotential (real basis, spin up)',all(abs(Vp) < 1e-10,'all'));

    Vp                  =   Vks.applyPotential(p,false);
    P                   =   zeros(size(Vks.potentialDown));
    for k = 1:x.nV, P   =   P + p(k) * x.basis(:,k);                        end
    for k = 1:x.nV, Vp(k)=  Vp(k) - sum(x.basis(:,k).*Vks.potentialDown.*P)*x.dx;   end
    obj.showResult('applyPotential (real basis, spin down)',all(abs(Vp) < 1e-10,'all'));

    v                   =   v.*[exp( 2i*x.x),exp( 3i*x.x),exp(-4i*x.x),exp(-2i*x.x)];
    x.set('basis',v);   x.orthonormalizeBasis;  x.initialize(DFT);
    Vks.add(zeros(size(Vks.potentialUp)));      Vks.initialize(x);
    p                   =   rand(randStr,[x.basisSize,1]) + 1i*rand(randStr,[x.basisSize,1]);
    Vp                  =   Vks.applyPotential(p,true);
    P                   =   zeros(size(Vks.potentialUp));
    for k = 1:x.nV, P   =   P + p(k) * x.basis(:,k);                        end
    for k = 1:x.nV, Vp(k)=  Vp(k) - sum(conj(x.basis(:,k)).*Vks.potentialUp.*P)*x.dx;   end
    obj.showResult('applyPotential (complex basis, spin up)',all(abs(Vp) < 1e-10,'all'));

    Vp                  =   Vks.applyPotential(p,false);
    P                   =   zeros(size(Vks.potentialDown));
    for k = 1:x.nV, P   =   P + p(k) * x.basis(:,k);                        end
    for k = 1:x.nV, Vp(k)=  Vp(k) - sum(conj(x.basis(:,k)).*Vks.potentialDown.*P)*x.dx;   end
    obj.showResult('applyPotential (complex basis, spin down)',all(abs(Vp) < 1e-10,'all'));

    Vks.add(@(p,isUp) p+isUp*p.^2);     Vks.add(@(p,isUp) p.^3+2*(~isUp)*p);    Vks.initialize(x);
    p                   =   rand(randStr,[x.basisSize,1]) + 1i*rand(randStr,[x.basisSize,1]);
    Vp                  =   Vks.applyPotential(p,true);
    P                   =   zeros(size(Vks.potentialUp));
    for k = 1:x.nV, P   =   P + p(k) * x.basis(:,k);                        end
    for k = 1:x.nV, Vp(k)=  Vp(k) - sum(conj(x.basis(:,k)).*(Vks.potentialUp.*P+P+P.^2+P.^3))*x.dx;   end
    obj.showResult('applyPotential (implicit, spin up)',all(abs(Vp) < 1e-10,'all'));

    Vp                  =   Vks.applyPotential(p,false);
    P                   =   zeros(size(Vks.potentialDown));
    for k = 1:x.nV, P   =   P + p(k) * x.basis(:,k);                        end
    for k = 1:x.nV, Vp(k)=  Vp(k) - sum(conj(x.basis(:,k)).*(Vks.potentialDown.*P+3*P+P.^3))*x.dx;   end
    obj.showResult('applyPotential (implicit, spin down)',all(abs(Vp) < 1e-10,'all'));

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

