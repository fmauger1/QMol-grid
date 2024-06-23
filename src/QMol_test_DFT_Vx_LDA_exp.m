classdef QMol_test_DFT_Vx_LDA_exp < QMol_test
%QMol_test_DFT_Vx_LDA_exp unit test for the QMol_DFT_Vx_LDA_exp class

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
    fprintf('  * QMol_test_DFT_Vx_LDA_exp\n'); 
    QMol_test_DFT_Vx_LDA_exp.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test units
    obj.runSpinRestricted;
    obj.runSpinPolarized;
    
end
end
methods (Access=private)
function runSpinRestricted(obj) %==========================================
%runSpinRestricted unit tests for spin-restricted DFT
    
    % Initialization
    obj.showSection('Spin-restricted DFT');

    A                   =   0.8;
    s                   =   3;
    Vx                  =   QMol_DFT_Vx_LDA_exp(                            ...
                                'potentialHeight',       A,                  ...
                                'potentialWidth',       s);
    DFT                 =   QMol_DFT_spinRes(...
                                'xspan',                -15:.3:20,          ...
                                'occupation',           [5 5 1.3],          ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vx);
    DFT.initialize;

    rho                 =   DFT.disc.DFT_allocateDensity;
    rho.set('rho',10.^linspace(-5,3,numel(DFT.x)).');
    rho_D               =   DFT.disc.DFT_allocateDensity;
    rho_D.set('rho',exp(-.5*(DFT.x(:)-3).^2));

    funV                =   @(r)  A  *exp(-r/s);
    funDV               =   @(r) -A/s*exp(-r/s);

    % Potential -- no SIC ------------------
    DFT.set('SIC','none');  DFT.initialize;

    V_th                =   obj.DFT_1D_LDA_eps_x(funV,rho.rho) + rho.rho.*obj.DFT_1D_LDA_d_eps_x(funDV,rho.rho);
    V                   =   Vx.getPotential(rho);

    R                   =   isa(V,'QMol_DFT_Vks')                           && ... potential object
                            all(abs( (V.potential - V_th)./V_th ) < 1e-6);       % potential discretization
    obj.showResult('getPotential (no SIC and new)' ,R);
    
    % ADSIC
    DFT.set('SIC','averageDensity'); DFT.initialize;

    N                   =   sum(DFT.occupation);
    V_th                =   2*V_th - obj.DFT_1D_LDA_eps_x(funV,rho.rho*2/N) - 2/N*rho.rho.*obj.DFT_1D_LDA_d_eps_x(funDV,rho.rho*2/N);
    Vx.getPotential(rho,V,true);
    R                   =   isa(V,'QMol_DFT_Vks')                           && ... potential object
                            all(abs( (V.potential - V_th)./V_th ) < 1e-6);       % potential discretization
    obj.showResult('getPotential (ADSIC and add to existing potential)' ,R);

    % Potential derivative -- no SIC ------
    DFT.set('SIC','none');           DFT.initialize;
    V_D                 =   Vx.getPotential(rho_D);
    DV_D                =   Vx.getPotentialDerivative(1,rho_D);
    R                   =   isa(DV_D,'QMol_DFT_Vks_grad') && isa(V_D,'QMol_DFT_Vks') && ... potential objects
                            all(abs(DV_D.potentialGradient - ifft(DFT.disc.D.*fft(V_D.potential))) < 1e-10);
    obj.showResult('getPotentialDerivative (no SIC and new)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity');    DFT.initialize;
    V_D                 =   Vx.getPotential(rho_D,V_D,true);
    DV_D                =   Vx.getPotentialDerivative(1,rho_D,DV_D,true);
    R                   =   isa(DV_D,'QMol_DFT_Vks_grad') && isa(V_D,'QMol_DFT_Vks') && ... potential objects
                            all(abs(DV_D.potentialGradient - ifft(DFT.disc.D.*fft(V_D.potential))) < 1e-10);
    obj.showResult('getPotentialDerivative (ADSIC and add to exist. pot. grad.)' ,R);

    % Energy -- no SIC ---------------------
    DFT.set('SIC','none');              DFT.initialize;
    E                   =   obj.DFT_1D_LDA_eps_x(funV,rho_D.rho);
    R                   =   Vx.getEnergy(rho_D) - sum(rho_D.rho.*E)*DFT.disc.dx;
    obj.showResult('getEnergy (no SIC)' ,R < 1e-9);

    % ADSIC
    DFT.set('SIC','averageDensity');    DFT.initialize;
    N                   =   sum(DFT.occupation);
    E                   =   obj.DFT_1D_LDA_eps_x(funV,rho_D.rho) - obj.DFT_1D_LDA_eps_x(funV,2*rho_D.rho/N);
    R                   =   Vx.getEnergy(rho_D) - sum(rho_D.rho.*E)*DFT.disc.dx;
    obj.showResult('getEnergy (no ADSIC)' ,R < 1e-9);

end
function runSpinPolarized(obj) %===========================================
%runSpinRestricted unit tests for spin-restricted DFT
    
    % Initialization
    obj.showSection('Spin-polarized DFT');

    A                   =   0.8;
    s                   =   3;
    Vx                  =   QMol_DFT_Vx_LDA_exp(                            ...
                                'potentialHeight',       A,                  ...
                                'potentialWidth',       s);
    DFT                 =   QMol_DFT_spinPol(...
                                'xspan',                -15:.3:20,          ...
                                'occupation',           {[5 2.4],[2 2 .1]}, ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vx);
    DFT.initialize;
    rho                 =   DFT.disc.DFT_allocateDensity;
    rho.set('rhoUp',10.^linspace(-5.01,2.97,numel(DFT.x)).');
    rho.set('rhoDw',10.^linspace(-4.52,2.87,numel(DFT.x)).');

    rho_D               =   DFT.disc.DFT_allocateDensity;
    rho_D.set('rhoUp',exp(-.5*(DFT.x(:)-3).^2),'rhoDw',exp(-.4*(DFT.x(:)-2).^2));

    funV                =   @(r)  A  *exp(-r/s);
    funDV               =   @(r) -A/s*exp(-r/s);

    % Poptential -- no SIC -----------------
    DFT.set('SIC','none');  DFT.initialize;

    V_th_up             =   obj.DFT_1D_LDA_eps_x(funV,2*rho.rhoUp) + 2*rho.rhoUp.*obj.DFT_1D_LDA_d_eps_x(funDV,2*rho.rhoUp);
    V_th_dw             =   obj.DFT_1D_LDA_eps_x(funV,2*rho.rhoDw) + 2*rho.rhoDw.*obj.DFT_1D_LDA_d_eps_x(funDV,2*rho.rhoDw);
    V                   =   Vx.getPotential(rho);
    R                   =   isa(V,'QMol_DFT_Vks')                                       && ... potential object
                            all(abs( (V.potentialUp   - V_th_up)./V_th_up ) < 1e-5)     && ... potential values
                            all(abs( (V.potentialDown - V_th_dw)./V_th_dw ) < 1e-5);
    obj.showResult('getPotential (no SIC and new)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity'); DFT.initialize;

    N                   =   [sum(DFT.occupation{1}), sum(DFT.occupation{2})];

    V_th_up             =   2*V_th_up - obj.DFT_1D_LDA_eps_x(funV,2*rho.rhoUp/N(1)) - 2*rho.rhoUp.*obj.DFT_1D_LDA_d_eps_x(funDV,2*rho.rhoUp/N(1))/N(1);
    V_th_dw             =   2*V_th_dw - obj.DFT_1D_LDA_eps_x(funV,2*rho.rhoDw/N(2)) - 2*rho.rhoDw.*obj.DFT_1D_LDA_d_eps_x(funDV,2*rho.rhoDw/N(2))/N(2);
    Vx.getPotential(rho,V,true);
    R                   =   isa(V,'QMol_DFT_Vks')                                       && ... potential object
                            all(abs( (V.potentialUp   - V_th_up)./V_th_up ) < 2e-5) && ... potential values
                            all(abs( (V.potentialDown - V_th_dw)./V_th_dw ) < 2e-5);
    obj.showResult('getPotential (ADSIC and add to existing potential)' ,R);

    % Potential derivative -- no SIC ------
    DFT.set('SIC','none');              DFT.initialize;
    V_D                 =   Vx.getPotential(rho_D);
    DV_D                =   Vx.getPotentialDerivative(1,rho_D);
    R                   =   isa(DV_D,'QMol_DFT_Vks_grad') && isa(V_D,'QMol_DFT_Vks')                               && ... potential objects
                            all(abs(DV_D.potentialGradientUp - ifft(DFT.disc.D.*fft(V_D.potentialUp)))    < 1e-10) && ...
                            all(abs(DV_D.potentialGradientDown - ifft(DFT.disc.D.*fft(V_D.potentialDown))) < 1e-10);
    obj.showResult('getPotentialDerivative (no SIC and new)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity');    DFT.initialize;
    V_D                 =   Vx.getPotential(rho_D,V_D,true);
    DV_D                =   Vx.getPotentialDerivative(1,rho_D,DV_D,true);
    R                   =   isa(DV_D,'QMol_DFT_Vks_grad') && isa(V_D,'QMol_DFT_Vks')                               && ... potential objects
                            all(abs(DV_D.potentialGradientUp - ifft(DFT.disc.D.*fft(V_D.potentialUp)))    < 1e-10) && ...
                            all(abs(DV_D.potentialGradientDown - ifft(DFT.disc.D.*fft(V_D.potentialDown))) < 1e-10);
    obj.showResult('getPotentialDerivative (ADSIC and add to exist. pot. grad.)' ,R);

    % Energy -- no SIC ---------------------
    DFT.set('SIC','none');              DFT.initialize;
    Eup                 =   obj.DFT_1D_LDA_eps_x(funV,2*rho_D.rhoUp);
    Edw                 =   obj.DFT_1D_LDA_eps_x(funV,2*rho_D.rhoDw);
    R                   =   Vx.getEnergy(rho_D) - sum(rho_D.rhoUp.*Eup+rho_D.rhoDw.*Edw)*DFT.disc.dx;
    obj.showResult('getEnergy (no SIC)' ,R < 1e-9);

    % ADSIC
    DFT.set('SIC','averageDensity');    DFT.initialize;
    N                   =   [sum(DFT.occupation{1}), sum(DFT.occupation{2})];
    Eup                 =   obj.DFT_1D_LDA_eps_x(funV,2*rho_D.rhoUp) - obj.DFT_1D_LDA_eps_x(funV,2*rho_D.rhoUp/N(1));
    Edw                 =   obj.DFT_1D_LDA_eps_x(funV,2*rho_D.rhoDw) - obj.DFT_1D_LDA_eps_x(funV,2*rho_D.rhoDw/N(2));
    R                   =   Vx.getEnergy(rho_D) - sum(rho_D.rhoUp.*Eup+rho_D.rhoDw.*Edw)*DFT.disc.dx;
    obj.showResult('getEnergy (no ADSIC)' ,R < 1e-9);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

