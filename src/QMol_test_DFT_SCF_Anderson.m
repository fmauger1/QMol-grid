classdef QMol_test_DFT_SCF_Anderson < QMol_test
%QMol_test_DFT_SCF_Anderson suite of unit tests for QMol_DFT_SCF_Anderson

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
    fprintf('  * QMol_test_DFT_SCF_Anderson\n'); 
    QMol_test_DFT_SCF_Anderson.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test units
    obj.test_spin_restricted_grid;
    obj.test_spin_restricted_basis;
    obj.test_spin_polarized_grid;
    obj.test_spin_polarized_basis;
    
end
end
methods (Access=private)
function test_spin_restricted_grid(obj) %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Initialization
    obj.showSection('Spin restricted DFT model (grid discretization)');

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
    
    % Start from scratch
    SCF                 =   QMol_DFT_SCF_Anderson('display',false);
    SCF.solveSCF(DFT);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (no symmetry, from scratch)',max(abs(DE))<1e-10);

    % Start from density
    SCF.set('mixingMode','potential');
    rho                 =   DFT.getDensity;
    rho.rho             =   rho.rho + flip(rho.rho);
    SCF.solveSCF(DFT,[],rho);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (no symmetry, from density)',max(abs(DE))<1e-10);

    % Start from potential
    SCF.set('convergenceTest','eigenvalues');
    Vks                 =   DFT.getPotential;
    Vks.V               =   .9*Vks.V;
    SCF.solveSCF(DFT,[],Vks);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (no symmetry, from potential)',max(abs(DE))<1e-10);

    % Symmetry imposed ====================================================
    DFT.set('xspan',-20:.1:20,'occupation',[2 1 1 ]);
    Vext.set('atom',{QMol_Va_softCoulomb('name','(1)','Z',2,'a',.5,'X0',-1.5), ...
                     QMol_Va_softCoulomb('name','(2)','Z',2,'a',.5,'X0', 1.5)});
    DFT.initialize;

    ES                  =   QMol_DFT_eigs('Symmetry','Sx + 2 Ax');

    % Start from scratch
    SCF.clear;              SCF.set('display',false);
    SCF.solveSCF(DFT,ES);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (imposed symmetry, from scratch)', ...
        max(abs(DFT.KSO.KSO(:,1)-flip(DFT.KSO.KSO(:,1)))) < 1e-10   && ...  % Symmetry
        max(abs(DFT.KSO.KSO(:,2)+flip(DFT.KSO.KSO(:,2)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSO(:,3)+flip(DFT.KSO.KSO(:,3)))) < 1e-10   && ...
        max(abs(DE))<1e-10);                                                % Converged

    % Start from density
    SCF.set('mixingMode','potential');
    rho                 =   DFT.getDensity;
    rho.rho             =   rho.rho + exp(-DFT.x(:).^2/5);
    SCF.solveSCF(DFT,ES,rho);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (imposed symmetry, from density)', ...
        max(abs(DFT.KSO.KSO(:,1)-flip(DFT.KSO.KSO(:,1)))) < 1e-10   && ...  % Symmetry
        max(abs(DFT.KSO.KSO(:,2)+flip(DFT.KSO.KSO(:,2)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSO(:,3)+flip(DFT.KSO.KSO(:,3)))) < 1e-10   && ...
        max(abs(DE))<1e-10);                                                % Converged

    % Start from potential
    SCF.set('convergenceTest','eigenvalues');
    Vks                 =   DFT.getPotential;
    Vks.V               =   .9*Vks.V;
    SCF.solveSCF(DFT,ES,Vks);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (imposed symmetry, from potential)', ...
        max(abs(DFT.KSO.KSO(:,1)-flip(DFT.KSO.KSO(:,1)))) < 1e-10   && ...  % Symmetry
        max(abs(DFT.KSO.KSO(:,2)+flip(DFT.KSO.KSO(:,2)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSO(:,3)+flip(DFT.KSO.KSO(:,3)))) < 1e-10   && ...
        max(abs(DE))<1e-10);                                                % Converged
end
function test_spin_restricted_basis(obj) %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Initialization
    obj.showSection('Spin restricted DFT model (basis discretization)');

    % Grid and basis set
    x                   =   -15:.1:20;

    A1                  =   QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5);
    A2                  =   QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3);
    A3                  =   QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi);

    AO                  =   @(s,n,x0,x) (x(:)-x0).^n .*exp(-(x(:)-x0).^2 * .5/s^2);
    V                   =  [AO(1.3,0,A1.position,x),AO(1.7,1,A1.position,x),AO(0.8,2,A1.position,x), ...
                            AO(1.3,0,A2.position,x),AO(1.7,1,A2.position,x),AO(0.8,2,A2.position,x), ...
                            AO(1.3,0,A3.position,x),AO(1.7,1,A3.position,x),AO(0.8,2,A3.position,x)];


    disc                =   QMol_disc_basis('x',x,'basis',V);
    disc.initialize;        disc.orthonormalizeBasis;

    % SCF iterations
    Vext                =   QMol_DFT_Vext('atom',{A1,A2,A3});
    Vh                  =   QMol_DFT_Vh_conv;
    Vxc                 =   QMol_DFT_Vx_LDA_exp;

    DFT                 =   QMol_DFT_spinRes(...
                                'discretization',                   disc,       ...
                                'selfInteractionCorrection',        'ADSIC',    ...
                                'occupation',                       [2 1 1],    ...
                                'externalPotential',                Vext,       ...
                                'HartreePotential',                 Vh,         ...
                                'exchangeCorrelationPotential',     Vxc);
    DFT.initialize;
    
    % Start from scratch
    SCF                 =   QMol_DFT_SCF_Anderson('display',false);
    SCF.solveSCF(DFT);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (from scratch)',max(abs(DE))<1e-10);

    % Start from density
    SCF.set('mixingMode','potential');
    rho                 =   DFT.getDensity;
    rho.rho             =   rho.rho + flip(rho.rho);
    SCF.solveSCF(DFT,[],rho);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (from density)',max(abs(DE))<1e-10);

    % Start from potential
    SCF.set('convergenceTest','eigenvalues');
    Vks                 =   DFT.getPotential;
    Vks.V               =   .9*Vks.V;
    SCF.solveSCF(DFT,[],Vks);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (from potential)',max(abs(DE))<1e-10);
end
function test_spin_polarized_grid(obj) %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Initialization
    obj.showSection('Spin polarized DFT model (grid discretization)');

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
    
    % Start from scratch
    SCF                 =   QMol_DFT_SCF_Anderson('display',false);
    SCF.solveSCF(DFT);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (no symmetry, from scratch)', ...
        max(abs(DE{1}))<1e-10   &&   max(abs(DE{2}))<1e-10);

    % Start from density
    SCF.set('mixingMode','potential');
    rho                 =   DFT.getDensity;
    rho.rhoUp           =   rho.rhoUp + 1.2*flip(rho.rhoUp);
    rho.rhoDw           =   rho.rhoDw + 1.8*flip(rho.rhoDw);
    SCF.solveSCF(DFT,[],rho);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (no symmetry, from density)', ...
        max(abs(DE{1}))<1e-10   &&   max(abs(DE{2}))<1e-10);

    % Start from potential
    SCF.set('convergenceTest','eigenvalues');
    Vks                 =   DFT.getPotential;
    Vks.Vup             =   .9*Vks.Vup;
    Vks.Vdw             =   1.1*Vks.Vdw;
    SCF.solveSCF(DFT,[],Vks);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (no symmetry, from potential)', ...
        max(abs(DE{1}))<1e-10   &&   max(abs(DE{2}))<1e-10);

    % Symmetry imposed ====================================================
    DFT.set('xspan',-20:.1:20,'occupation',{[1 .5 .5],[1 1]});
    Vext.set('atom',{QMol_Va_softCoulomb('name','(1)','Z',2,'a',.5,'X0',-1.5), ...
                     QMol_Va_softCoulomb('name','(2)','Z',2,'a',.5,'X0', 1.5)});
    DFT.initialize;

    ES                  =   QMol_DFT_eigs('Symmetry',{'Sx + 2 Ax','Sx + Ax'});

    % Start from scratch
    SCF.clear;              SCF.set('display',false);
    SCF.solveSCF(DFT,ES);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (imposed symmetry, from scratch)', ...
        max(abs(DFT.KSO.KSOup(:,1)-flip(DFT.KSO.KSOup(:,1)))) < 1e-10   && ...  % Symmetry
        max(abs(DFT.KSO.KSOup(:,2)+flip(DFT.KSO.KSOup(:,2)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSOup(:,3)+flip(DFT.KSO.KSOup(:,3)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSOdw(:,1)-flip(DFT.KSO.KSOdw(:,1)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSOdw(:,2)+flip(DFT.KSO.KSOdw(:,2)))) < 1e-10   && ...
        max(abs(DE{1}))<1e-10   &&   max(abs(DE{2}))<1e-10);                    % Converged

    % Start from density
    SCF.set('mixingMode','potential');
    rho                 =   DFT.getDensity;
    rho.rhoUp           =   rho.rhoUp + exp(-DFT.x(:).^2/5);
    rho.rhoDw           =   rho.rhoDw + exp(-DFT.x(:).^2/6);
    SCF.solveSCF(DFT,ES,rho);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (imposed symmetry, from density)', ...
        max(abs(DFT.KSO.KSOup(:,1)-flip(DFT.KSO.KSOup(:,1)))) < 1e-10   && ...  % Symmetry
        max(abs(DFT.KSO.KSOup(:,2)+flip(DFT.KSO.KSOup(:,2)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSOup(:,3)+flip(DFT.KSO.KSOup(:,3)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSOdw(:,1)-flip(DFT.KSO.KSOdw(:,1)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSOdw(:,2)+flip(DFT.KSO.KSOdw(:,2)))) < 1e-10   && ...
        max(abs(DE{1}))<1e-10   &&   max(abs(DE{2}))<1e-10);                    % Converged

    % Start from potential
    SCF.set('convergenceTest','eigenvalues');
    Vks                 =   DFT.getPotential;
    Vks.Vup             =   .8*Vks.Vup;
    Vks.Vdw             =   1.2*Vks.Vdw;
    SCF.solveSCF(DFT,ES,Vks);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (imposed symmetry, from potential)', ...
        max(abs(DFT.KSO.KSOup(:,1)-flip(DFT.KSO.KSOup(:,1)))) < 1e-10   && ...  % Symmetry
        max(abs(DFT.KSO.KSOup(:,2)+flip(DFT.KSO.KSOup(:,2)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSOup(:,3)+flip(DFT.KSO.KSOup(:,3)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSOdw(:,1)-flip(DFT.KSO.KSOdw(:,1)))) < 1e-10   && ...
        max(abs(DFT.KSO.KSOdw(:,2)+flip(DFT.KSO.KSOdw(:,2)))) < 1e-10   && ...
        max(abs(DE{1}))<1e-10   &&   max(abs(DE{2}))<1e-10);                    % Converged
end
function test_spin_polarized_basis(obj) %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Initialization
    obj.showSection('Spin polarized DFT model (basis discretization)');

    % Grid and basis set
    x                   =   -15:.1:20;

    A1                  =   QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5);
    A2                  =   QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3);
    A3                  =   QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi);

    AO                  =   @(s,n,x0,x) (x(:)-x0).^n .*exp(-(x(:)-x0).^2 * .5/s^2);
    V                   =  [AO(1.3,0,A1.position,x),AO(1.7,1,A1.position,x),AO(0.8,2,A1.position,x), ...
                            AO(1.3,0,A2.position,x),AO(1.7,1,A2.position,x),AO(0.8,2,A2.position,x), ...
                            AO(1.3,0,A3.position,x),AO(1.7,1,A3.position,x),AO(0.8,2,A3.position,x)];


    disc                =   QMol_disc_basis('x',x,'basis',V);
    disc.initialize;        disc.orthonormalizeBasis;
    Vext                =   QMol_DFT_Vext('atom',{A1,A2,A3});
    Vh                  =   QMol_DFT_Vh_conv;
    Vxc                 =   QMol_DFT_Vx_LDA_exp;

    DFT                 =   QMol_DFT_spinPol(...
                                'discretization',                   disc,  ...
                                'selfInteractionCorrection',        'ADSIC',    ...
                                'occupation',                       {[1 .5 .5],[1 1]},    ...
                                'externalPotential',                Vext,       ...
                                'HartreePotential',                 Vh,         ...
                                'exchangeCorrelationPotential',     Vxc);
    DFT.initialize;
    
    % Start from scratch
    SCF                 =   QMol_DFT_SCF_Anderson('display',false);
    SCF.solveSCF(DFT);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (from scratch)', ...
        max(abs(DE{1}))<1e-10   &&   max(abs(DE{2}))<1e-10);

    % Start from density
    SCF.set('mixingMode','potential');
    rho                 =   DFT.getDensity;
    rho.rhoUp           =   rho.rhoUp + 1.2*flip(rho.rhoUp);
    rho.rhoDw           =   rho.rhoDw + 1.8*flip(rho.rhoDw);
    SCF.solveSCF(DFT,[],rho);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (from density)', ...
        max(abs(DE{1}))<1e-10   &&   max(abs(DE{2}))<1e-10);

    % Start from potential
    SCF.set('convergenceTest','eigenvalues');
    Vks                 =   DFT.getPotential;
    Vks.Vup             =   .9*Vks.Vup;
    Vks.Vdw             =   1.1*Vks.Vdw;
    SCF.solveSCF(DFT,[],Vks);

    [~,DE]              =   DFT.getEnergy('KSO');
    obj.showResult('solveSCF (from potential)', ...
        max(abs(DE{1}))<1e-10   &&   max(abs(DE{2}))<1e-10);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

