classdef QMol_test_DFT_spinRes < QMol_test
%QMol_test_DFT_spinRes suite of unit tests for QMol_DFT_spinRes

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
    fprintf('  * QMol_test_DFT_spinRes\n'); 
    QMol_test_DFT_spinRes.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test units
    obj.test_initialization;
    obj.test_DFT_density_and_potential;
    obj.test_DFT_energy;
    
end
end
methods (Access=private)
function test_initialization(obj) %========================================
    
    % DFT model
    obj.showSection('Initialization');
    
    X                   =   -20:.1:15;

    Vext                =   QMol_DFT_Vext('externalPotential',@(x)-5./sqrt(x.^2+.7));
    Vh                  =   QMol_DFT_Vh_conv;
    Vx                  =   QMol_DFT_Vx_LDA_exp;
    Vc                  =   QMol_DFT_Vx_LDA_exp('V0',1,'s',8);              % Replace with actual correlation potential

    DFT                 =   QMol_DFT_spinRes( ...
                                'xspan',                        X,                  ...
                                'occupation',                   [2 2 1],    ...
                                'externalPotential',            Vext, ...
                                'HartreePotential',             Vh, ...
                                'exchangeCorrelationPotential', {Vx,Vc}, ...
                                'selfInteractionCorrection',    'ADSIC');
    
    DFT.initialize;

    % discretization 
    obj.showResult('discretization',DFT.discretization == QMol_disc('xspan',X));

    % orbital
    obj.showResult('orbital', ...
        ~DFT.orbital.isSpinPol                                  && ... Spin restricted model
        all(size(DFT.orbital.orbital) == [numel(X) 3])          && ... Orbital is allocated
        isempty(DFT.orbital.orbitalUp)                          && ... Only for spin polarized models
        isempty(DFT.orbital.orbitalDown));

    % Total charge
    obj.showResult('totalCharge',DFT.totalCharge == 5);

    % isInitialized
    obj.showResult('isInitialized',DFT.isInitialized);

    % isSpinPol 
    obj.showResult('isSpinPol',~DFT.isSpinPol);

    % dim
    obj.showResult('dim',DFT.dim == 1);
    
end
function test_DFT_density_and_potential(obj) %=============================

    % DFT model
    obj.showSection('One-body density and Kohn-Sham potential');
    
    X                   =   -20:.1:15;

    Vext                =   QMol_DFT_Vext('externalPotential',@(x)-5./sqrt(x.^2+.7));
    Vh                  =   QMol_DFT_Vh_conv;
    Vx                  =   QMol_DFT_Vx_LDA_exp;
    Vc                  =   QMol_DFT_Vx_LDA_exp('V0',1,'s',8);              % Replace with actual correlation potential

    DFT                 =   QMol_DFT_spinRes( ...
                                'xspan',                        X,                  ...
                                'occupation',                   [2 2 1],    ...
                                'externalPotential',            Vext, ...
                                'HartreePotential',             Vh, ...
                                'exchangeCorrelationPotential', {Vx,Vc}, ...
                                'selfInteractionCorrection',    'ADSIC');
    
    DFT.initialize;

    %getDensity
    dx                  =   X(2)-X(1);      X = X.';
    p                   =  [      exp(-X.^2*.5    ) / sqrt(sum( (      exp(-X.^2*.5    )).^2 )*dx), ...
                            X   .*exp(-X.^2*.5/2^2) / sqrt(sum( (X   .*exp(-X.^2*.5/2^2)).^2 )*dx), ...
                            X.^2.*exp(-X.^2*.5/3^2) / sqrt(sum( (X.^2.*exp(-X.^2*.5/3^2)).^2 )*dx)];
    DFT.KSO.KSO         =   p;

    rho                 =   DFT.disc.DFT_allocateDensity;
    rho.rho             =   2*p(:,1).^2 + 2*p(:,2).^2 + p(:,3).^2;
    obj.showResult('getDensity',DFT.getDensity == rho);

    % getPotential
    Vks                 =   Vext.getPotential;
    Vh.getPotential([],Vks,true);
    Vx.getPotential([],Vks,true);
    Vc.getPotential([],Vks,true);
    V = DFT.getPotential; 
    obj.showResult('getPotential (DFT density)', ...
        max(abs(Vks.V-V.V)) < 1e-10);
    
    rho.rho             =   1.2*rho.rho;
    Vks                 =   Vext.getPotential;
    Vh.getPotential(rho,Vks,true);
    Vx.getPotential(rho,Vks,true);
    Vc.getPotential(rho,Vks,true);
    V = DFT.getPotential(rho); 
    obj.showResult('getPotential (rho supplied)', ...
        max(abs(Vks.V-V.V)) < 1e-10);

    Vks                 =   Vext.getPotential;
    Vh.getPotential([],Vks,true);
    Vx.getPotential([],Vks,true);
    Vc.getPotential([],Vks,true);
    DFT.getPotential([],V); 
    obj.showResult('getPotential (overwrite)', ...
        max(abs(Vks.V-V.V)) < 1e-10);

    %getPotentialGradient
    DVks                =   Vext.getPotentialDerivative(1);
    Vh.getPotentialDerivative(1,[],DVks,true);
    Vx.getPotentialDerivative(1,[],DVks,true);
    Vc.getPotentialDerivative(1,[],DVks,true);
    DV = DFT.getPotentialGradient; 
    obj.showResult('getPotentialGradient (DFT density)', ...
        max(abs(DVks.DV-DV.DV)) < 1e-10);
    
    rho.rho             =   1.2*rho.rho;
    DVks                =   Vext.getPotentialDerivative(1);
    Vh.getPotentialDerivative(1,rho,DVks,true);
    Vx.getPotentialDerivative(1,rho,DVks,true);
    Vc.getPotentialDerivative(1,rho,DVks,true);
    DV = DFT.getPotentialGradient(rho); 
    obj.showResult('getPotentialGradient (rho supplied)', ...
        max(abs(DVks.DV-DV.DV)) < 1e-10);

    DVks                =   Vext.getPotentialDerivative(1);
    Vh.getPotentialDerivative(1,[],DVks,true);
    Vx.getPotentialDerivative(1,[],DVks,true);
    Vc.getPotentialDerivative(1,[],DVks,true);
    DFT.getPotentialGradient([],DV); 
    obj.showResult('getPotentialGradient (overwrite)', ...
        max(abs(DVks.DV-DV.DV)) < 1e-10);

end
function test_DFT_energy(obj) %============================================

    % DFT model
    obj.showSection('Energy components');
    
    X                   =   -20:.1:15;

    Vext                =   QMol_DFT_Vext('externalPotential',@(x)-5./sqrt(x.^2+.7));
    Vh                  =   QMol_DFT_Vh_conv;
    Vx                  =   QMol_DFT_Vx_LDA_exp;
    Vc                  =   QMol_DFT_Vx_LDA_exp('V0',1,'s',8);              % Replace with actual correlation potential

    DFT                 =   QMol_DFT_spinRes( ...
                                'xspan',                        X,                  ...
                                'occupation',                   [2 2 1],    ...
                                'externalPotential',            Vext, ...
                                'HartreePotential',             Vh, ...
                                'exchangeCorrelationPotential', {Vx,Vc}, ...
                                'selfInteractionCorrection',    'ADSIC');
    
    DFT.initialize;

    dx                  =   X(2)-X(1);      X = X.';
    p                   =  [      exp(-X.^2*.5    ) / sqrt(sum( (      exp(-X.^2*.5    )).^2 )*dx), ...
                            X   .*exp(-X.^2*.5/2^2) / sqrt(sum( (X   .*exp(-X.^2*.5/2^2)).^2 )*dx), ...
                            X.^2.*exp(-X.^2*.5/3^2) / sqrt(sum( (X.^2.*exp(-X.^2*.5/3^2)).^2 )*dx)];
    DFT.KSO.KSO         =   p;

    % DFT energy
    [Etot,Ekin,Eext,Eh,Exc] = DFT.getEnergy('DFT');
    rho                 =   DFT.getDensity;

    obj.showResult('getEnergy (DFT -- total)',Etot == Ekin+Eext+Eh+Exc);

    E                   =   DFT.disc.DFT_energyKinetic(DFT.occ,DFT.KSO);
    obj.showResult('getEnergy (DFT -- kinetic)', E == Ekin);

    E                   =   Vext.getEnergy(rho);
    obj.showResult('getEnergy (DFT -- external)', E == Eext);

    E                   =   Vh.getEnergy(rho);
    obj.showResult('getEnergy (DFT -- Hartree)', E == Eh);

    E                   =   Vx.getEnergy(rho) + Vc.getEnergy(rho);
    obj.showResult('getEnergy (DFT -- Hartree)', E == Exc);

    % Orbital energy
    Vks                 =   DFT.getPotential;
    [E_1,DE_1]          =   DFT.getEnergy('orbital');
    [E_2,DE_2]          =   DFT.disc.DFT_energyOrbital(Vks,DFT.KSO);
    obj.showResult('getEnergy (orbital)', ...
        all(E_1 == E_2)   &&  all(DE_1 == DE_2));

    rho.rho             =   1.2*rho.rho;
    DFT.getPotential(rho,Vks);
    [E_1,DE_1]          =   DFT.getEnergy('orbital',rho);
    [E_2,DE_2]          =   DFT.disc.DFT_energyOrbital(Vks,DFT.KSO);
    obj.showResult('getEnergy (orbital -- rho supplied)', ...
        all(E_1 == E_2)   &&  all(DE_1 == DE_2));

    Vks.V               =   1.2*Vks.V;
    [E_1,DE_1]          =   DFT.getEnergy('orbital',Vks);
    [E_2,DE_2]          =   DFT.disc.DFT_energyOrbital(Vks,DFT.KSO);
    obj.showResult('getEnergy (orbital -- Vks supplied)', ...
        all(E_1 == E_2)   &&  all(DE_1 == DE_2));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

