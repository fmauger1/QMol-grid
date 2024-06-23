classdef QMol_test_SE_V < QMol_test
%QMol_test_SE_Vt suite of unit tests for QMol_SE_V

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
    fprintf('  * QMol_test_SE_V\n'); 
    QMol_test_SE_V.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test units
    obj.runDefinePotential;
    obj.runUsePotential;
    
end
end
methods (Access=private)
function runDefinePotential(obj) %=========================================
%runDefinePotential unit tests defining the external potential
    
    % Initialization
    obj.showSection('Potential discretization');

    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v_r                 =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];
    v_c                 =   [v_r.*[exp( 2i*x),exp( 3i*x),exp(-4i*x),exp(-2i*x)], exp(-(x-4).^2/2.5)];


    disc                =   QMol_disc(      'xspan',x);
    d_re                =   QMol_disc_basis('xspan',x,'basis',v_r);         d_re.orthonormalizeBasis;
    d_im                =   QMol_disc_basis('xspan',x,'basis',v_c);         d_im.orthonormalizeBasis;

    SE                  =   QMol_SE;
    Vse                 =   QMol_SE_V;

    % User-defined function
    V                   =   @(x) -3*x./sqrt((x-sqrt(2)).^2 + exp(2));
    Vse.set('inputPotential',V);
    Vin                 =   V(x);

    SE.set('discretization',disc);      disc.initialize(SE);        Vse.initialize(SE);
    obj.showResult('function handle (grid)' ,max(abs( Vse.V-Vin )) < 1e-10);

    SE.set('discretization',d_re);      d_re.initialize(SE);        Vse.initialize(SE);
    p                   =   rand(randStr,[d_re.basisSize,1],'like',1i);
    Vp                  =   Vse.mV*p;
    P                   =   zeros(size(Vse.potential));
    for k = 1:d_re.nV, P    =   P + p(k) * d_re.basis(:,k);                         end
    for k = 1:d_re.nV, Vp(k)=  Vp(k) - sum(d_re.basis(:,k).*Vin.*P)*d_re.dx;        end
    obj.showResult('function handle (real basis)',all(abs(Vp) < 1e-10,'all'));

    SE.set('discretization',d_im);      d_im.initialize(SE);        Vse.initialize(SE);
    p                   =   rand(randStr,[d_im.basisSize,1],'like',1i);
    Vp                  =   Vse.mV*p;
    P                   =   zeros(size(Vse.potential));
    for k = 1:d_im.nV, P    =   P + p(k) * d_im.basis(:,k);                         end
    for k = 1:d_im.nV, Vp(k)=  Vp(k) - sum(conj(d_im.basis(:,k)).*Vin.*P)*d_im.dx;  end
    obj.showResult('function handle (complex basis)',all(abs(Vp) < 1e-10,'all'));

    % User-defined discretization
    V                   =   V(disc.x(:));
    Vse.set('inputPotential',V);
    Vin                 =   V;

    SE.set('discretization',disc);      disc.initialize(SE);        Vse.initialize(SE);
    obj.showResult('user-defined discretization (grid)' ,max(abs( Vse.V-Vin )) < 1e-10);

    SE.set('discretization',d_re);      d_re.initialize(SE);        Vse.initialize(SE);
    p                   =   rand(randStr,[d_re.basisSize,1],'like',1i);
    Vp                  =   Vse.mV*p;
    P                   =   zeros(size(Vse.potential));
    for k = 1:d_re.nV, P    =   P + p(k) * d_re.basis(:,k);                         end
    for k = 1:d_re.nV, Vp(k)=  Vp(k) - sum(d_re.basis(:,k).*Vin.*P)*d_re.dx;        end
    obj.showResult('user-defined discretization (real basis)',all(abs(Vp) < 1e-10,'all'));

    SE.set('discretization',d_im);      d_im.initialize(SE);        Vse.initialize(SE);
    p                   =   rand(randStr,[d_im.basisSize,1],'like',1i);
    Vp                  =   Vse.mV*p;
    P                   =   zeros(size(Vse.potential));
    for k = 1:d_im.nV, P    =   P + p(k) * d_im.basis(:,k);                         end
    for k = 1:d_im.nV, Vp(k)=  Vp(k) - sum(conj(d_im.basis(:,k)).*Vin.*P)*d_im.dx;  end
    obj.showResult('user-defined discretization (complex basis)',all(abs(Vp) < 1e-10,'all'));

    % List of atomic centers
    Va                  =  {QMol_Va_Gaussian('V0',1,'s',.5,'X0',-3), ...
                            QMol_Va_Gaussian('V0',sqrt(2),'s',exp(1),'X0',log(5)),...
                            QMol_Va_softCoulomb('Z',1/3,'a',pi,'X0',5)};
    Vse.set('inputPotential',[],'atom',Va);
    Vin                 =   Va{1}.getPotential(disc.x(:)) + ...
                            Va{2}.getPotential(disc.x(:)) + ...
                            Va{3}.getPotential(disc.x(:));

    SE.set('discretization',disc);      disc.initialize(SE);        Vse.initialize(SE);
    obj.showResult('list of atomic centers (grid)' ,max(abs( Vse.V-Vin )) < 1e-10);

    SE.set('discretization',d_re);      d_re.initialize(SE);        Vse.initialize(SE);
    p                   =   rand(randStr,[d_re.basisSize,1],'like',1i);
    Vp                  =   Vse.mV*p;
    P                   =   zeros(size(Vse.potential));
    for k = 1:d_re.nV, P    =   P + p(k) * d_re.basis(:,k);                         end
    for k = 1:d_re.nV, Vp(k)=  Vp(k) - sum(d_re.basis(:,k).*Vin.*P)*d_re.dx;        end
    obj.showResult('list of atomic centers (real basis)',all(abs(Vp) < 1e-10,'all'));

    SE.set('discretization',d_im);      d_im.initialize(SE);        Vse.initialize(SE);
    p                   =   rand(randStr,[d_im.basisSize,1],'like',1i);
    Vp                  =   Vse.mV*p;
    P                   =   zeros(size(Vse.potential));
    for k = 1:d_im.nV, P    =   P + p(k) * d_im.basis(:,k);                         end
    for k = 1:d_im.nV, Vp(k)=  Vp(k) - sum(conj(d_im.basis(:,k)).*Vin.*P)*d_im.dx;  end
    obj.showResult('list of atomic centers (complex basis)',all(abs(Vp) < 1e-10,'all'));

    % Both user-defined and list of atoms
    V                   =   @(x) -3*x./sqrt((x-sqrt(2)).^2 + exp(2));
    Vse.set('inputPotential',V,'atom',Va{1});
    Vin                 =   V(disc.x(:)) + Va{1}.getPotential(disc.x(:));

    SE.set('discretization',disc);      disc.initialize(SE);        Vse.initialize(SE);
    obj.showResult('user-defined + llist of atomic centers (grid)' ,max(abs( Vse.V-Vin )) < 1e-10);

    SE.set('discretization',d_re);      d_re.initialize(SE);        Vse.initialize(SE);
    p                   =   rand(randStr,[d_re.basisSize,1],'like',1i);
    Vp                  =   Vse.mV*p;
    P                   =   zeros(size(Vse.potential));
    for k = 1:d_re.nV, P    =   P + p(k) * d_re.basis(:,k);                         end
    for k = 1:d_re.nV, Vp(k)=  Vp(k) - sum(d_re.basis(:,k).*Vin.*P)*d_re.dx;        end
    obj.showResult('user-defined + llist of atomic centers (real basis)',all(abs(Vp) < 1e-10,'all'));

    SE.set('discretization',d_im);      d_im.initialize(SE);        Vse.initialize(SE);
    p                   =   rand(randStr,[d_im.basisSize,1],'like',1i);
    Vp                  =   Vse.mV*p;
    P                   =   zeros(size(Vse.potential));
    for k = 1:d_im.nV, P    =   P + p(k) * d_im.basis(:,k);                         end
    for k = 1:d_im.nV, Vp(k)=  Vp(k) - sum(conj(d_im.basis(:,k)).*Vin.*P)*d_im.dx;  end
    obj.showResult('user-defined + llist of atomic centers (complex basis)',all(abs(Vp) < 1e-10,'all'));

    % Potential derivative ~~~~~~~~~~~~~~~~~
    obj.showSection('External potential derivative discretization');

    V                   =   @(x) -3*x./sqrt((x-sqrt(2)).^2 + exp(2));
    DV                  =   @(x) -3./sqrt((x-sqrt(2)).^2 + exp(2)) + 3*x.*(x-sqrt(2)).*((x-sqrt(2)).^2 + exp(2)).^-1.5;

    Vse.set('inputPotential',V,'atom',[]);
    Vse.initialize(SE);     Vse.setDerivative;
    obj.showResult('finite difference approximation',max(abs(Vse.DV-DV(disc.x(:)))) < 5e-10);
    
    Vse.set('inputPotentialDerivative',DV);
    Vse.initialize(SE);     Vse.setDerivative;
    obj.showResult('function handle',max(abs(Vse.DV-DV(disc.x(:)))) < 1e-10);

    DV                  =   DV(disc.x(:));
    Vse.set('inputPotentialDerivative',DV);
    Vse.initialize(SE);     Vse.setDerivative;
    obj.showResult('user-defined discretization',max(abs(Vse.DV-DV)) < 1e-10);
    
    Vse.set('inputPotential',[],'atom',Va);
    Vse.initialize(SE);     Vse.setDerivative;

    DV                  =   DV + Va{1}.getPotentialDerivative(1,disc.x(:)) + ...
                                 Va{2}.getPotentialDerivative(1,disc.x(:)) + ...
                                 Va{3}.getPotentialDerivative(1,disc.x(:));
    R                   =   max(abs( Vse.DV-DV )) < 1e-10;
    obj.showResult('list of atomic centers' ,R);
    if ~R, fprintf('      (check QMol_Va_Gaussian and QMol_Va_softCoulomb components)\n');    end

    
    SE.set('discretization',d_re);      d_re.initialize(SE);        Vse.initialize(SE);
    Vse.setDerivative;      DVin    =   Vse.DV;
    p                   =   rand(randStr,[d_re.basisSize,1],'like',1i);
    DVp                 =   Vse.mDV*p;
    P                   =   zeros(size(Vse.potentialDerivative));
    for k = 1:d_re.nV, P    =   P + p(k) * d_re.basis(:,k);                         end
    for k = 1:d_re.nV, DVp(k)=  DVp(k) - sum(d_re.basis(:,k).*DVin.*P)*d_re.dx;        end
    obj.showResult('potentialDerivativeMatrix (real basis)',all(abs(DVp) < 1e-10,'all'));
    
    SE.set('discretization',d_im);      d_im.initialize(SE);        Vse.initialize(SE);
    Vse.setDerivative;      DVin    =   Vse.DV;
    p                   =   rand(randStr,[d_im.basisSize,1],'like',1i);
    DVp                 =   Vse.mDV*p;
    P                   =   zeros(size(Vse.potentialDerivative));
    for k = 1:d_im.nV, P    =   P + p(k) * d_im.basis(:,k);                         end
    for k = 1:d_im.nV, DVp(k)=  DVp(k) - sum(conj(d_im.basis(:,k)).*DVin.*P)*d_im.dx;        end
    obj.showResult('potentialDerivativeMatrix (complex basis)',all(abs(DVp) < 1e-10,'all'));

end
function runUsePotential(obj) %=========================================
%runUsePotential unit tests using the external potential
    
    % Initialization

    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v_r                 =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];
    v_c                 =   [v_r.*[exp( 2i*x),exp( 3i*x),exp(-4i*x),exp(-2i*x)], exp(-(x-4).^2/2.5)];


    disc                =   QMol_disc(      'xspan',x);
    d_re                =   QMol_disc_basis('xspan',x,'basis',v_r);         d_re.orthonormalizeBasis;
    d_im                =   QMol_disc_basis('xspan',x,'basis',v_c);         d_im.orthonormalizeBasis;

    SE                  =   QMol_SE;
    Vse                 =   QMol_SE_V;

    % Potential energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('Potential energy');

    V                   =   @(x) -3*x./sqrt((x-sqrt(2)).^2 + exp(2));
    Vse.set('inputPotential',V);
    Vin                 =   V(x);

    SE.set('discretization',disc);      disc.initialize(SE);        Vse.initialize(SE);
    x0                  =   5;
    s                   =   2;
    psi                 =   QMol_SE_wfcn('wavefunction',...
                               [exp(-(x-x0  ).^2 * .5/s^2), ...
                                exp(-(x-2*x0).^2 * .8/s^2), ...
                                exp(-(x+x0  ).^2 * .3/s^2)]);
    E                   =   sum(Vin.*sum(psi.wfcn.^2,2))*(x(2)-x(1));
    obj.showResult('getEnergy (grid)',abs(E-Vse.getEnergy(psi)) < 1e-10);

    SE.set('discretization',d_re);      d_re.initialize(SE);        Vse.initialize(SE);
    p                   =   rand(randStr,[d_re.basisSize,8],'like',1i);
    psi                 =   QMol_SE_wfcn_basis('wavefunction',p);
    P                   =   d_re.basis*p;
    E                   =   sum(Vin.*sum(abs(P).^2,2))*(x(2)-x(1));
    obj.showResult('getEnergy (real basis)',abs(E-Vse.getEnergy(psi)) < 1e-10);

    SE.set('discretization',d_im);      d_im.initialize(SE);        Vse.initialize(SE);
    p                   =   rand(randStr,[d_im.basisSize,7],'like',1i);
    psi                 =   QMol_SE_wfcn_basis('wavefunction',p);
    P                   =   d_im.basis*p;
    E                   =   sum(Vin.*sum(abs(P).^2,2))*(x(2)-x(1));
    obj.showResult('getEnergy (complex basis)',abs(E-Vse.getEnergy(psi)) < 1e-10);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

