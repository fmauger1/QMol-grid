function test_spin_restricted(obj)

    % Initialization
    obj.showSection('Spin restricted DFT model');

    % No symmetry imposed =================================================
    Vext                =   QMol_DFT_Vext('atom', ...
                               {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5), ...
                                QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3), ...
                                QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)});
    Vh                  =   QMol_DFT_Vh_conv;
    Vxc                 =   QMol_DFT_Vx_LDA_exp;

    DFT                 =   QMol_DFT_spinRes(...
                                'xspan',                            -15:.1:20,  ...
                                'selfInteractionCorrection',        'ADSIC',    ...
                                'occupation',                       [2 1 1],    ...
                                'externalPotential',                Vext,       ...
                                'HartreePotential',                 Vh,         ...
                                'exchangeCorrelationPotential',     Vxc);
    DFT.initialize;

    ES                  =   QMol_DFT_eigs;
    ES.initialize(DFT);
    
    KSO                 =   DFT.KSO;
    x                   =   DFT.x(:);
    KSO.KSO             =  [exp(-x.^2), x.*exp(-(x-3).^2/2), cos(2*x).*exp(-(x-1).^2/1.5)];
    KSO.KSO             =   DFT.disc.DFT_normalizeOrbital(KSO.KSO);
    
    Vks                 =   DFT.getPotential;
    E                   =   ES.computeEigenstates;

    R_size              =   all(size(KSO.KSO) == [numel(x) 3]);

    dx                  =   x(2)-x(1);
    R_norm              =   abs( sum(abs(KSO.KSO(:,1)).^2)*dx - 1) < 1e-10  &&  ... % KSO are normalized
                            abs( sum(abs(KSO.KSO(:,2)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(KSO.KSO(:,3)).^2)*dx - 1) < 1e-10;

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSO(:,1));
    E1                  =   sum( KSO.KSO(:,1) .* Hp ) * dx;
    DE1                 =   sqrt(sum( abs(Hp - E(1)*KSO.KSO(:,1)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSO(:,2));
    E2                  =   sum( KSO.KSO(:,2) .* Hp ) * dx;
    DE2                 =   sqrt(sum( abs(Hp - E(2)*KSO.KSO(:,2)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSO(:,3));
    E3                  =   sum( KSO.KSO(:,3) .* Hp ) * dx;
    DE3                 =   sqrt(sum( abs(Hp - E(3)*KSO.KSO(:,3)).^2 ) *dx);

    R_energy            =   max(abs(E - [E1;E2;E3])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3]) < 1e-10;
    
    obj.showResult('computeEigenstates (no symmetry)',R_size && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of Kohn-Sham orbitals\n'); end
    if ~R_norm,     fprintf('      - Kohn-Sham orbitals are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Kohn-Sham orbitals are not converged\n'); end


    Vks                 =   DFT.getPotential;
    Vks.V               =   .9*Vks.V;

    Ee                  =   ES.computeEigenstates(Vks);
    [E,DE]              =   DFT.getEnergy('orbital',Vks);

    obj.showResult('computeEigenstates (no symmetry, input potential)', ...
        max(abs(Ee-E)) < 1e-10 && max(abs(DE)) < 1e-10);

    % Symmetry imposed ====================================================
    DFT.set('xspan',-20:.1:20,'occupation',[2 1 1 0]);
    Vext.set('atom',{QMol_Va_softCoulomb('name','(1)','Z',2,'a',.5,'X0',-1.5), ...
                     QMol_Va_softCoulomb('name','(2)','Z',2,'a',.5,'X0', 1.5)});
    DFT.initialize;

    x                   =   DFT.x(:);
    dx                  =   x(2)-x(1);
    rho                 =   DFT.disc.DFT_allocateDensity;
    rho.rho             =   exp(-x.^2 * .5/3^2);
    rho.rho             =   4 * rho.rho / sum(rho.rho) / dx;

    ES.set('Symmetry','1 Sx + 3 Ax');
    ES.initialize(DFT);
    
    DFT.getPotential(rho,Vks);
    E                   =   ES.computeEigenstates(Vks);

    R_size              =   all(size(KSO.KSO) == [numel(x) 4]);
    
    R_norm              =   abs( sum(abs(KSO.KSO(:,1)).^2)*dx - 1) < 1e-10  &&  ... % KSO are normalized
                            abs( sum(abs(KSO.KSO(:,2)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(KSO.KSO(:,3)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(KSO.KSO(:,4)).^2)*dx - 1) < 1e-10;
    
    R_symmetry          =   max(abs(KSO.KSO(:,1)-flip(KSO.KSO(:,1)))) < 1e-10 && ... KSO symmetry
                            max(abs(KSO.KSO(:,2)+flip(KSO.KSO(:,2)))) < 1e-10 && ...
                            max(abs(KSO.KSO(:,3)+flip(KSO.KSO(:,3)))) < 1e-10 && ...
                            max(abs(KSO.KSO(:,4)+flip(KSO.KSO(:,4)))) < 1e-10;

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSO(:,1));
    E1                  =   sum( KSO.KSO(:,1) .* Hp ) * dx;
    DE1                 =   sqrt(sum( abs(Hp - E(1)*KSO.KSO(:,1)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSO(:,2));
    E2                  =   sum( KSO.KSO(:,2) .* Hp ) * dx;
    DE2                 =   sqrt(sum( abs(Hp - E(2)*KSO.KSO(:,2)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSO(:,3));
    E3                  =   sum( KSO.KSO(:,3) .* Hp ) * dx;
    DE3                 =   sqrt(sum( abs(Hp - E(3)*KSO.KSO(:,3)).^2 ) *dx);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(Vks,KSO.KSO(:,4));
    E4                  =   sum( KSO.KSO(:,4) .* Hp ) * dx;
    DE4                 =   sqrt(sum( abs(Hp - E(4)*KSO.KSO(:,4)).^2 ) *dx);

    R_energy            =   max(abs(E - [E1;E2;E3;E4])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3 DE4]) < 1e-10;

    obj.showResult('computeEigenstates (imposed symmetry)',R_size && R_symmetry && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of Kohn-Sham orbitals\n'); end
    if ~R_size,     fprintf('      - Wrong Kohn-Sham-orbital symmetries\n'); end
    if ~R_norm,     fprintf('      - Kohn-Sham orbitals are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Kohn-Sham orbitals are not converged\n'); end

    
    Vks                 =   DFT.getPotential;
    Vks.V               =   .9*Vks.V;

    Ee                  =   ES.computeEigenstates(Vks);
    [E,DE]              =   DFT.getEnergy('orbital',Vks);

    obj.showResult('computeEigenstates (imposed symmetry, input potential)', ...
        max(abs(Ee-E)) < 1e-10 && max(abs(DE)) < 1e-10);
end