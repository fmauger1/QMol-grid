classdef QMol_test_DFT_eig_basis < QMol_test
%QMol_test_DFT_eig_basis unit test for the QMol_DFT_eig_basis class in 1D

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
    fprintf('  * QMol_test_DFT_eig_basis\n'); 
    QMol_test_DFT_eig_basis.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test units
    obj.test_spin_restricted;
    obj.test_spin_polarized;
    
end
end
%% Test components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=private)
    function test_spin_restricted(obj) %===================================
%test_spin_restricted
    
    % Initialization
    obj.showSection('Spin-restricted DFT model');
    
    % Real basis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];
    disc                =   QMol_disc_basis('xspan',x,'basis',v);   disc.orthonormalizeBasis;

    Vext                =   QMol_DFT_Vext('atom', ...
                               {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5), ...
                                QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3), ...
                                QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)});
    Vh                  =   QMol_DFT_Vh_conv;
    Vxc                 =   QMol_DFT_Vx_LDA_exp;

    DFT                 =   QMol_DFT_spinRes(...
                                'disc',                             disc,  ...
                                'selfInteractionCorrection',        'ADSIC',    ...
                                'occupation',                       [2 1 1],    ...
                                'externalPotential',                Vext,       ...
                                'HartreePotential',                 Vh,         ...
                                'exchangeCorrelationPotential',     Vxc);
    DFT.initialize;

    ES                  =   QMol_DFT_eig_basis;
    ES.initialize(DFT);
    
    x                   =   DFT.x(:);
    KSO                 =   QMol_DFT_orbital('KSO',[exp(-x.^2), x.*exp(-(x-3).^2/2), cos(2*x).*exp(-(x-1).^2/1.5)]);
    DFT.KSO             =   disc.DFT_projectOrbital(KSO);
    DFT.KSO.KSO         =   disc.DFT_normalizeOrbital(DFT.KSO.KSO);
    DFT.initialize;
    
    Vks                 =   DFT.getPotential;
    E                   =   ES.computeEigenstates;

    R_size              =   all(size(DFT.KSO.KSO) == [disc.basisSize 3]);

    R_norm              =   abs( sum(abs(DFT.KSO.KSO(:,1)).^2) - 1) < 1e-10  &&  ... % KSO are normalized
                            abs( sum(abs(DFT.KSO.KSO(:,2)).^2) - 1) < 1e-10  &&  ...
                            abs( sum(abs(DFT.KSO.KSO(:,3)).^2) - 1) < 1e-10;

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSO(:,1));
    E1                  =   sum( DFT.KSO.KSO(:,1) .* Hp );
    DE1                 =   sqrt(sum( abs(Hp - E(1)*DFT.KSO.KSO(:,1)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSO(:,2));
    E2                  =   sum( DFT.KSO.KSO(:,2) .* Hp );
    DE2                 =   sqrt(sum( abs(Hp - E(2)*DFT.KSO.KSO(:,2)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSO(:,3));
    E3                  =   sum( DFT.KSO.KSO(:,3) .* Hp );
    DE3                 =   sqrt(sum( abs(Hp - E(3)*DFT.KSO.KSO(:,3)).^2 ));

    R_energy            =   max(abs(E - [E1;E2;E3])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3]) < 1e-10;
    
    obj.showResult('computeEigenstates (real basis)',R_size && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of Kohn-Sham orbitals\n'); end
    if ~R_norm,     fprintf('      - Kohn-Sham orbitals are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Kohn-Sham orbitals are not converged\n'); end

    % Complex basis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)] .* ...
                            [exp( 2i*x),    exp( 3i*x),       exp(-4i*x),exp(-2i*x)];
    disc.set('basis',v);    disc.orthonormalizeBasis;
    DFT.reset();            DFT.initialize;

    ES.reset;               ES.initialize(DFT);
    KSO.KSO             =   KSO.KSO .* [exp(-1i*x), exp(2i*x), exp(1i*x)];
    DFT.KSO             =   disc.DFT_projectOrbital(KSO);
    DFT.KSO.KSO         =   disc.DFT_normalizeOrbital(DFT.KSO.KSO);
    DFT.initialize;
    
    Vks                 =   DFT.getPotential;
    E                   =   ES.computeEigenstates;

    R_size              =   all(size(DFT.KSO.KSO) == [disc.basisSize 3]);

    R_norm              =   abs( sum(abs(DFT.KSO.KSO(:,1)).^2) - 1) < 1e-10  &&  ... % KSO are normalized
                            abs( sum(abs(DFT.KSO.KSO(:,2)).^2) - 1) < 1e-10  &&  ...
                            abs( sum(abs(DFT.KSO.KSO(:,3)).^2) - 1) < 1e-10;

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSO(:,1));
    E1                  =   sum( conj(DFT.KSO.KSO(:,1)) .* Hp );
    DE1                 =   sqrt(sum( abs(Hp - E(1)*DFT.KSO.KSO(:,1)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSO(:,2));
    E2                  =   sum( conj(DFT.KSO.KSO(:,2)) .* Hp );
    DE2                 =   sqrt(sum( abs(Hp - E(2)*DFT.KSO.KSO(:,2)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSO(:,3));
    E3                  =   sum( conj(DFT.KSO.KSO(:,3)) .* Hp );
    DE3                 =   sqrt(sum( abs(Hp - E(3)*DFT.KSO.KSO(:,3)).^2 ));

    R_energy            =   max(abs(E - [E1;E2;E3])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3]) < 1e-10;

    obj.showResult('computeEigenstates (complex basis)',R_size && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of Kohn-Sham orbitals\n'); end
    if ~R_norm,     fprintf('      - Kohn-Sham orbitals are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Kohn-Sham orbitals are not converged\n'); end

end
function test_spin_polarized(obj) %===================================
%test_spin_polarized
    
    % Initialization
    obj.showSection('Spin-polarized DFT model');
    
    % Real basis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];
    disc                =   QMol_disc_basis('xspan',x,'basis',v);   disc.orthonormalizeBasis;

    Vext                =   QMol_DFT_Vext('atom', ...
                               {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5), ...
                                QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3), ...
                                QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)});
    Vh                  =   QMol_DFT_Vh_conv;
    Vxc                 =   QMol_DFT_Vx_LDA_exp;

    DFT                 =   QMol_DFT_spinPol(...
                                'disc',                             disc,  ...
                                'selfInteractionCorrection',        'ADSIC',    ...
                                'occupation',                       {[1 .5 .5],[1 1]},    ...
                                'externalPotential',                Vext,       ...
                                'HartreePotential',                 Vh,         ...
                                'exchangeCorrelationPotential',     Vxc);
    DFT.initialize;

    ES                  =   QMol_DFT_eig_basis;
    ES.initialize(DFT);
    
    x                   =   DFT.x(:);
    KSO                 =   QMol_DFT_orbital( ...
                                'KSOup',[exp(-x.^2), x.*exp(-(x-3).^2/2), cos(2*x).*exp(-(x-1).^2/1.5)], ...
                                'KSOdw',[exp(-(x+1).^2/2), x.*exp(-x.^2/1.5)]);
    DFT.KSO             =   disc.DFT_projectOrbital(KSO);
    DFT.KSO.KSO         =   disc.DFT_normalizeOrbital(DFT.KSO.KSO);
    DFT.initialize;
    
    Vks                 =   DFT.getPotential;
    E                   =   ES.computeEigenstates;
    Eup                 =   E{1};
    Edw                 =   E{2};

    R_size              =   all(size(DFT.KSO.KSOup) == [disc.basisSize 3]) && ...
                            all(size(DFT.KSO.KSOdw) == [disc.basisSize 2]);

    R_norm              =   abs( sum(abs(DFT.KSO.KSOup(:,1)).^2) - 1) < 1e-10  &&  ... % KSO are normalized
                            abs( sum(abs(DFT.KSO.KSOup(:,2)).^2) - 1) < 1e-10  &&  ...
                            abs( sum(abs(DFT.KSO.KSOup(:,3)).^2) - 1) < 1e-10  &&  ...
                            abs( sum(abs(DFT.KSO.KSOdw(:,1)).^2) - 1) < 1e-10  &&  ...
                            abs( sum(abs(DFT.KSO.KSOdw(:,2)).^2) - 1) < 1e-10;

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSOup(:,1),true);
    E1                  =   sum( DFT.KSO.KSOup(:,1) .* Hp );
    DE1                 =   sqrt(sum( abs(Hp - Eup(1)*DFT.KSO.KSOup(:,1)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSOup(:,2),true);
    E2                  =   sum( DFT.KSO.KSOup(:,2) .* Hp );
    DE2                 =   sqrt(sum( abs(Hp - Eup(2)*DFT.KSO.KSOup(:,2)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSOup(:,3),true);
    E3                  =   sum( DFT.KSO.KSOup(:,3) .* Hp );
    DE3                 =   sqrt(sum( abs(Hp - Eup(3)*DFT.KSO.KSOup(:,3)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSOdw(:,1),false);
    E4                  =   sum( DFT.KSO.KSOdw(:,1) .* Hp );
    DE4                 =   sqrt(sum( abs(Hp - Edw(1)*DFT.KSO.KSOdw(:,1)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSOdw(:,2),false);
    E5                  =   sum( DFT.KSO.KSOdw(:,2) .* Hp );
    DE5                 =   sqrt(sum( abs(Hp - Edw(2)*DFT.KSO.KSOdw(:,2)).^2 ));

    R_energy            =   max(abs(Eup - [E1;E2;E3])) < 1e-10  && ...
                            max(abs(Edw - [E4;E5   ])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3 DE4 DE5]) < 1e-10;
    
    obj.showResult('computeEigenstates (real basis)',R_size && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of Kohn-Sham orbitals\n'); end
    if ~R_norm,     fprintf('      - Kohn-Sham orbitals are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Kohn-Sham orbitals are not converged\n'); end

    % Complex basis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)] .* ...
                            [exp( 2i*x),    exp( 3i*x),       exp(-4i*x),exp(-2i*x)];
    disc.set('basis',v);    disc.orthonormalizeBasis;
    DFT.reset();            DFT.initialize;

    ES.reset;               ES.initialize(DFT);
    KSO.KSOup           =   KSO.KSOup .* [exp(-1i*x), exp(2i*x), exp(1i*x)];
    KSO.KSOdw           =   KSO.KSOdw .* [exp(2i*x), exp(1i*x)];
    DFT.KSO             =   disc.DFT_projectOrbital(KSO);
    DFT.KSO.KSO         =   disc.DFT_normalizeOrbital(DFT.KSO.KSO);
    DFT.initialize;
    
    Vks                 =   DFT.getPotential;
    E                   =   ES.computeEigenstates;
    Eup                 =   E{1};
    Edw                 =   E{2};

    R_size              =   all(size(DFT.KSO.KSOup) == [disc.basisSize 3]) && ...
                            all(size(DFT.KSO.KSOdw) == [disc.basisSize 2]);

    R_norm              =   abs( sum(abs(DFT.KSO.KSOup(:,1)).^2) - 1) < 1e-10  &&  ... % KSO are normalized
                            abs( sum(abs(DFT.KSO.KSOup(:,2)).^2) - 1) < 1e-10  &&  ...
                            abs( sum(abs(DFT.KSO.KSOup(:,3)).^2) - 1) < 1e-10  &&  ...
                            abs( sum(abs(DFT.KSO.KSOdw(:,1)).^2) - 1) < 1e-10  &&  ...
                            abs( sum(abs(DFT.KSO.KSOdw(:,2)).^2) - 1) < 1e-10;

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSOup(:,1),true);
    E1                  =   sum( conj(DFT.KSO.KSOup(:,1)) .* Hp );
    DE1                 =   sqrt(sum( abs(Hp - Eup(1)*DFT.KSO.KSOup(:,1)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSOup(:,2),true);
    E2                  =   sum( conj(DFT.KSO.KSOup(:,2)) .* Hp );
    DE2                 =   sqrt(sum( abs(Hp - Eup(2)*DFT.KSO.KSOup(:,2)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSOup(:,3),true);
    E3                  =   sum( conj(DFT.KSO.KSOup(:,3)) .* Hp );
    DE3                 =   sqrt(sum( abs(Hp - Eup(3)*DFT.KSO.KSOup(:,3)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSOdw(:,1),false);
    E4                  =   sum( conj(DFT.KSO.KSOdw(:,1)) .* Hp );
    DE4                 =   sqrt(sum( abs(Hp - Edw(1)*DFT.KSO.KSOdw(:,1)).^2 ));

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,DFT.KSO.KSOdw(:,2),false);
    E5                  =   sum( conj(DFT.KSO.KSOdw(:,2)) .* Hp );
    DE5                 =   sqrt(sum( abs(Hp - Edw(2)*DFT.KSO.KSOdw(:,2)).^2 ));

    R_energy            =   max(abs(Eup - [E1;E2;E3])) < 1e-10  && ...
                            max(abs(Edw - [E4;E5   ])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3 DE4 DE5]) < 1e-10;
    
    obj.showResult('computeEigenstates (complex basis)',R_size && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of Kohn-Sham orbitals\n'); end
    if ~R_norm,     fprintf('      - Kohn-Sham orbitals are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Kohn-Sham orbitals are not converged\n'); end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

