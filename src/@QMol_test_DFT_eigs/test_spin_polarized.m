function test_spin_polarized(obj)

    % Initialization
    obj.showSection('Spin restricted DFT model');

    % No symmetry imposed =================================================
    Vext                =   QMol_DFT_Vext('atom', ...
                               {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5), ...
                                QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3), ...
                                QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)});
    Vh                  =   QMol_DFT_Vh_conv;
    Vxc                 =   QMol_DFT_Vx_LDA_exp;

    DFT                 =   QMol_DFT_spinPol(...
                                'xspan',                            -15:.1:20,  ...
                                'selfInteractionCorrection',        'ADSIC',    ...
                                'occupation',                       {[1 .5 .5],[1 1]},    ...
                                'externalPotential',                Vext,       ...
                                'HartreePotential',                 Vh,         ...
                                'exchangeCorrelationPotential',     Vxc);
    DFT.initialize;

    ES                  =   QMol_DFT_eigs;
    ES.initialize(DFT);
    
    KSO                 =   DFT.KSO;
    x                   =   DFT.x(:);
    KSO.KSOup           =  [exp(-x.^2), x.*exp(-(x-3).^2/2), cos(2*x).*exp(-(x-1).^2/1.5)];
    KSO.KSOup           =   DFT.disc.DFT_normalizeOrbital(KSO.KSOup);
    KSO.KSOdw           =  [exp(-(x-1).^2), x.*exp(-(x+1).^2/2)];
    KSO.KSOdw           =   DFT.disc.DFT_normalizeOrbital(KSO.KSOdw);
    
    Vks                 =   DFT.getPotential;
    E                   =   ES.computeEigenstates;
    Eup                 =   E{1};
    Edw                 =   E{2};

    R_size              =   all(size(KSO.KSOup) == [numel(x) 3])    && ...
                            all(size(KSO.KSOdw) == [numel(x) 2]);

    dx                  =   x(2)-x(1);
    R_norm              =   abs( sum(abs(KSO.KSOup(:,1)).^2)*dx - 1) < 1e-10  &&  ... % KSO are normalized
                            abs( sum(abs(KSO.KSOup(:,2)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(KSO.KSOup(:,3)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(KSO.KSOdw(:,1)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(KSO.KSOdw(:,2)).^2)*dx - 1) < 1e-10;

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOup(:,1),true);
    E1                  =   sum( KSO.KSOup(:,1) .* Hp ) * dx;
    DE1                 =   sqrt(sum( abs(Hp - Eup(1)*KSO.KSOup(:,1)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOup(:,2),true);
    E2                  =   sum( KSO.KSOup(:,2) .* Hp ) * dx;
    DE2                 =   sqrt(sum( abs(Hp - Eup(2)*KSO.KSOup(:,2)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOup(:,3),true);
    E3                  =   sum( KSO.KSOup(:,3) .* Hp ) * dx;
    DE3                 =   sqrt(sum( abs(Hp - Eup(3)*KSO.KSOup(:,3)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOdw(:,1),false);
    E4                  =   sum( KSO.KSOdw(:,1) .* Hp ) * dx;
    DE4                 =   sqrt(sum( abs(Hp - Edw(1)*KSO.KSOdw(:,1)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOdw(:,2),false);
    E5                  =   sum( KSO.KSOdw(:,2) .* Hp ) * dx;
    DE5                 =   sqrt(sum( abs(Hp - Edw(2)*KSO.KSOdw(:,2)).^2 ) *dx);

    R_energy            =   max(abs(Eup - [E1;E2;E3])) < 1e-10  && ...
                            max(abs(Edw - [E4;E5   ])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3 DE4 DE5]) < 1e-10;
    
    obj.showResult('computeEigenstates (no symmetry)',R_size && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of Kohn-Sham orbitals\n'); end
    if ~R_norm,     fprintf('      - Kohn-Sham orbitals are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Kohn-Sham orbitals are not converged\n'); end


    Vks                 =   DFT.getPotential;
    Vks.Vup             =   1.2*Vks.Vup;
    Vks.Vdw             =   0.9*Vks.Vdw;

    Ee                  =   ES.computeEigenstates(Vks);
    [E,DE]              =   DFT.getEnergy('orbital',Vks);

    obj.showResult('computeEigenstates (no symmetry, input potential)', ...
        max(abs(Ee{1}-E{1})) < 1e-10 && max(abs(Ee{2}-E{2})) < 1e-10 && ... Energies
        max(abs(DE{1})) < 1e-10 && max(abs(DE{2})) < 1e-10);              % Errors

    % Symmetry imposed ====================================================
    DFT.set('xspan',-20:.1:20,'occupation',{[1 1 1 0],[.5 .5 0]});
    Vext.set('atom',{QMol_Va_softCoulomb('name','(1)','Z',2,'a',.5,'X0',-1.5), ...
                     QMol_Va_softCoulomb('name','(2)','Z',2,'a',.5,'X0', 1.5)});
    DFT.initialize;

    x                   =   DFT.x(:);
    dx                  =   x(2)-x(1);
    rho                 =   DFT.disc.DFT_allocateDensity;
    rho.rhoUp           =   exp(-x.^2 * .5/3^2);
    rho.rhoUp           =   2 * rho.rhoUp / sum(rho.rhoUp) / dx;
    rho.rhoDw           =   exp(-(x-1).^2 * .5/2^2) + exp(-(x+1).^2 * .5/2^2);
    rho.rhoDw           =   2 * rho.rhoDw / sum(rho.rhoDw) / dx;

    ES.set('Symmetry',{'1 Sx + 3 Ax','2 Sx + Ax'});

    ES.initialize(DFT);
    DFT.getPotential(rho,Vks);
    E                   =   ES.computeEigenstates(Vks);
    Eup                 =   E{1};
    Edw                 =   E{2};

    R_size              =   all(size(KSO.KSOup) == [numel(x) 4])    && ...
                            all(size(KSO.KSOdw) == [numel(x) 3]);
    
    dx                  =   x(2)-x(1);
    KSO                 =   DFT.KSO;
    R_norm              =   abs( sum(abs(KSO.KSOup(:,1)).^2)*dx - 1) < 1e-10  &&  ... % KSO are normalized
                            abs( sum(abs(KSO.KSOup(:,2)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(KSO.KSOup(:,3)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(KSO.KSOup(:,4)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(KSO.KSOdw(:,1)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(KSO.KSOdw(:,2)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(KSO.KSOdw(:,3)).^2)*dx - 1) < 1e-10;

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOup(:,1),true);
    E1                  =   sum( KSO.KSOup(:,1) .* Hp ) * dx;
    DE1                 =   sqrt(sum( abs(Hp - Eup(1)*KSO.KSOup(:,1)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOup(:,2),true);
    E2                  =   sum( KSO.KSOup(:,2) .* Hp ) * dx;
    DE2                 =   sqrt(sum( abs(Hp - Eup(2)*KSO.KSOup(:,2)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOup(:,3),true);
    E3                  =   sum( KSO.KSOup(:,3) .* Hp ) * dx;
    DE3                 =   sqrt(sum( abs(Hp - Eup(3)*KSO.KSOup(:,3)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOup(:,4),true);
    E4                  =   sum( KSO.KSOup(:,4) .* Hp ) * dx;
    DE4                 =   sqrt(sum( abs(Hp - Eup(4)*KSO.KSOup(:,4)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOdw(:,1),false);
    E5                  =   sum( KSO.KSOdw(:,1) .* Hp ) * dx;
    DE5                 =   sqrt(sum( abs(Hp - Edw(1)*KSO.KSOdw(:,1)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOdw(:,2),false);
    E6                  =   sum( KSO.KSOdw(:,2) .* Hp ) * dx;
    DE6                 =   sqrt(sum( abs(Hp - Edw(2)*KSO.KSOdw(:,2)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSOdw(:,3),false);
    E7                  =   sum( KSO.KSOdw(:,3) .* Hp ) * dx;
    DE7                 =   sqrt(sum( abs(Hp - Edw(3)*KSO.KSOdw(:,3)).^2 ) *dx);

    R_energy            =   max(abs(Eup - [E1;E2;E3;E4])) < 1e-10  && ...
                            max(abs(Edw - [E5;E6;E7   ])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3 DE4 DE5 DE6 DE7]) < 1e-10;
    
    obj.showResult('computeEigenstates (no symmetry)',R_size && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of Kohn-Sham orbitals\n'); end
    if ~R_norm,     fprintf('      - Kohn-Sham orbitals are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Kohn-Sham orbitals are not converged\n'); end

    
    Vks                 =   DFT.getPotential;
    Vks.Vup             =   1.2*Vks.Vup;
    Vks.Vdw             =   0.9*Vks.Vdw;

    Ee                  =   ES.computeEigenstates(Vks);
    [E,DE]              =   DFT.getEnergy('orbital',Vks);

    obj.showResult('computeEigenstates (imposed symmetry, input potential)', ...
        max(abs(Ee{1}-E{1})) < 1e-10 && max(abs(Ee{2}-E{2})) < 1e-10 && ... Energies
        max(abs(DE{1})) < 1e-10 && max(abs(DE{2})) < 1e-10);              % Errors
end