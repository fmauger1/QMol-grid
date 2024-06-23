classdef QMol_test_DFT_Vx_LDA_soft < QMol_test
%QMol_test_DFT_Vx_LDA_soft unit test for the QMol_DFT_Vx_LDA_soft class

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
    fprintf('  * QMol_test_DFT_Vx_LDA_soft\n'); 
    QMol_test_DFT_Vx_LDA_soft.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test units
    obj.runFunctionalFit;
    obj.runSpinRestricted;
    obj.runSpinPolarized;
    
end
end
methods (Access=private)
function runFunctionalFit(obj) %===========================================
%runFunctionalFit unit test for the quality of the functional fit (as
%   compared to the analytical integral formula)
    
    % Initialization 
    obj.showSection('Functional fit');

    Z                   =   .8;
    a                   =   1.3;

    % Spin restricted ----------------------
    Vx                  =   QMol_DFT_Vx_LDA_soft('Z',Z,'a',a);
    DFT                 =   QMol_DFT_spinRes(...
                                'xspan',                -15:.3:20,          ...
                                'occupation',           [5.2 3.4],          ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vx);
    DFT.initialize;
    rho                 =   DFT.disc.DFT_allocateDensity;
    rho.set('rho',10.^linspace(-2.99,2.02,numel(DFT.x)).');

    funV                =   @(r)  Z  ./sqrt(r.^2+a^2);
    funDV               =   @(r) -Z*r.*(r.^2+a^2).^-1.5;

    % No SIC
    DFT.set('SIC','none');  DFT.initialize;

    V_th                =   obj.DFT_1D_LDA_eps_x(funV,rho.rho) + rho.rho.*obj.DFT_1D_LDA_d_eps_x(funDV,rho.rho);
    V                   =   Vx.getPotential(rho);
    R                   =   isa(V,'QMol_DFT_Vks')                           && ... potential object
                            max(abs( (V.potential - V_th)./V_th )) < 2e-3;       % potential discretization
    obj.showResult('spin restricted (no SIC)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity'); DFT.initialize;

    N                   =   sum(DFT.occupation);
    V_th                =   2*V_th - obj.DFT_1D_LDA_eps_x(funV,rho.rho*2/N) - 2/N*rho.rho.*obj.DFT_1D_LDA_d_eps_x(funDV,rho.rho*2/N);
    Vx.getPotential(rho,V,true);
    R                   =   isa(V,'QMol_DFT_Vks')                           && ... potential object
                            max(abs( (V.potential - V_th)./V_th )) < 4e-3;       % potential discretization
    obj.showResult('spin restricted (no ADSIC)' ,R);

    clear DFT Vx

    % Spin-polarized -----------------------
    Vx                  =   QMol_DFT_Vx_LDA_soft('Z',Z,'a',a);
    DFT                 =   QMol_DFT_spinPol(...
                                'xspan',                -15:.3:20,          ...
                                'occupation',           {[5 2.4],[2 2 .1]}, ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vx);
    DFT.initialize;
    rho                 =   DFT.disc.DFT_allocateDensity;
    rho.set('rhoUp',10.^linspace(-3.01,2.03,numel(DFT.x)).');
    rho.set('rhoDw',10.^linspace(-2.52,1.87,numel(DFT.x)).');
    
    % No SIC
    DFT.set('SIC','none');  DFT.initialize;

    V_th_up             =   obj.DFT_1D_LDA_eps_x(funV,2*rho.rhoUp) + 2*rho.rhoUp.*obj.DFT_1D_LDA_d_eps_x(funDV,2*rho.rhoUp);
    V_th_dw             =   obj.DFT_1D_LDA_eps_x(funV,2*rho.rhoDw) + 2*rho.rhoDw.*obj.DFT_1D_LDA_d_eps_x(funDV,2*rho.rhoDw);
    V                   =   Vx.getPotential(rho);
    R                   =   isa(V,'QMol_DFT_Vks')                                       && ... potential object
                            max(abs( (V.potentialUp   - V_th_up)./V_th_up )) < 2e-3 && ... potential values
                            max(abs( (V.potentialDown - V_th_dw)./V_th_dw )) < 2e-3;
    obj.showResult('spin polarized (no SIC)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity'); DFT.initialize;

    N                   =   [sum(DFT.occupation{1}), sum(DFT.occupation{2})];
    V_th_up             =   2*V_th_up - obj.DFT_1D_LDA_eps_x(funV,2*rho.rhoUp/N(1)) - 2*rho.rhoUp.*obj.DFT_1D_LDA_d_eps_x(funDV,2*rho.rhoUp/N(1))/N(1);
    V_th_dw             =   2*V_th_dw - obj.DFT_1D_LDA_eps_x(funV,2*rho.rhoDw/N(2)) - 2*rho.rhoDw.*obj.DFT_1D_LDA_d_eps_x(funDV,2*rho.rhoDw/N(2))/N(2);
    Vx.getPotential(rho,V,true);
    R                   =   isa(V,'QMol_DFT_Vks')                                       && ... potential object
                            max(abs( (V.potentialUp   - V_th_up)./V_th_up )) < 4e-3 && ... potential values
                            max(abs( (V.potentialDown - V_th_dw)./V_th_dw )) < 4e-3;
    obj.showResult('spin polarized (ADSIC)' ,R);

end
function runSpinRestricted(obj) %==========================================
%runSpinRestricted unit tests for spin-restricted DFT
    
    % Initialization
    obj.showSection('Spin-restricted DFT');

    Z                   =   .8;
    a                   =   1.3;

    alp_fit             =   10.18001817;
    bet_fit             =   5.502143989;
    gam_fit             =   14.64700068;
    m_fit               =   2.301803657;

    Vx                  =   QMol_DFT_Vx_LDA_soft('Z',Z,'a',a);
    DFT                 =   QMol_DFT_spinRes(...
                                'xspan',                -15:.3:20,          ...
                                'occupation',           [5.2 3.4],          ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vx);
    DFT.initialize;
    rho                 =   DFT.disc.DFT_allocateDensity;
    rho.set('rho',10.^linspace(-2.99,2.02,numel(DFT.x)).');

    rho_D               =   DFT.disc.DFT_allocateDensity;
    rho_D.set('rho',exp(-.5*DFT.x(:).^2));

    dx                  =   DFT.x(2)-DFT.x(1);

    % Potential -- no SIC ------------------
    DFT.set('SIC','none');  DFT.initialize;

    V_th                =   Z*obj.V_fit(  alp_fit,bet_fit,gam_fit,m_fit,.5/a./rho.rho)/a;
    V                   =   Vx.getPotential(rho);
    R                   =   isa(V,'QMol_DFT_Vks')                           && ... potential object
                            max(abs( (V.potential - V_th)./V_th )) < 1e-10;      % potential discretization
    obj.showResult('getPotential (no SIC and new)' ,R);
    
    % ADSIC
    DFT.set('SIC','averageDensity'); DFT.initialize;

    N                   =   sum(DFT.occupation);
    V_th                =   V_th ...
                           +Z*obj.V_fit(  alp_fit,bet_fit,gam_fit,m_fit,.50  /a./rho.rho)/a ...
                           -Z*obj.V_fit(  alp_fit,bet_fit,gam_fit,m_fit,.25*N/a./rho.rho)/a;
    Vx.getPotential(rho,V,true);
    R                   =   isa(V,'QMol_DFT_Vks')                           && ... potential object
                            max(abs( (V.potential - V_th)./V_th )) < 1e-10;       % potential discretization
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

    E_th                =   Z*obj.eps_fit(alp_fit,bet_fit,gam_fit,m_fit,.5/a./rho_D.rho)/a;
    E                   =   Vx.getEnergy(rho_D);
    R                   =   abs(E - sum(rho_D.rho.*E_th)*dx ) < 1e-10;
    obj.showResult('getEnergy (no SIC)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity');    DFT.initialize;

    N                   =   sum(DFT.occupation);
    E_th                =   Z*obj.eps_fit(alp_fit,bet_fit,gam_fit,m_fit,.50  /a./rho_D.rho)/a ...
                           -Z*obj.eps_fit(alp_fit,bet_fit,gam_fit,m_fit,.25*N/a./rho_D.rho)/a;
    E                   =   Vx.getEnergy(rho_D);
    R                   =   abs(E - sum(rho_D.rho.*E_th)*dx ) < 1e-10;
    obj.showResult('getEnergy (ADSIC)' ,R);

end
function runSpinPolarized(obj) %===========================================
%runSpinRestricted unit tests for spin-restricted DFT
    
    % Initialization
    obj.showSection('Spin-polarized DFT');

    Z                   =   .8;
    a                   =   1.3;

    alp_fit             =   10.18001817;
    bet_fit             =   5.502143989;
    gam_fit             =   14.64700068;
    m_fit               =   2.301803657;

    Vx                  =   QMol_DFT_Vx_LDA_soft('Z',Z,'a',a);
    DFT                 =   QMol_DFT_spinPol(...
                                'xspan',                -15:.3:20,          ...
                                'occupation',           {[5 2.4],[2 2 .1]}, ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vx);
    DFT.initialize;
    rho                 =   DFT.disc.DFT_allocateDensity;
    rho.set('rhoUp',10.^linspace(-3.01,2.03,numel(DFT.x)).');
    rho.set('rhoDw',10.^linspace(-2.52,1.87,numel(DFT.x)).');

    rho_D               =   DFT.disc.DFT_allocateDensity;
    rho_D.set('rhoUp',exp(-.5*DFT.x(:).^2),'rhoDw',exp(-.4*DFT.x(:).^2));

    dx                  =   DFT.x(2)-DFT.x(1);

    % Poptential -- no SIC -----------------
    DFT.set('SIC','none');  DFT.initialize;

    V_th                =  [Z*obj.V_fit(  alp_fit,bet_fit,gam_fit,m_fit,.25/a./rho.rhoUp)/a, ...
                            Z*obj.V_fit(  alp_fit,bet_fit,gam_fit,m_fit,.25/a./rho.rhoDw)/a];

    V                   =   Vx.getPotential(rho);
    R                   =   isa(V,'QMol_DFT_Vks')                                       && ... potential object
                            max(abs( (V.potentialUp   - V_th(:,1))./V_th(:,1) )) < 1e-10 && ... potential values
                            max(abs( (V.potentialDown - V_th(:,2))./V_th(:,2) )) < 1e-10;
    obj.showResult('getPotential (no SIC and new)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity'); DFT.initialize;

    N                   =   [sum(DFT.occupation{1}), sum(DFT.occupation{2})];
    V_th                =   V_th + ...
                           [Z*obj.V_fit(alp_fit,bet_fit,gam_fit,m_fit,.25     /a./rho.rhoUp)/a, ...
                            Z*obj.V_fit(alp_fit,bet_fit,gam_fit,m_fit,.25     /a./rho.rhoDw)/a] ...
                          -[Z*obj.V_fit(alp_fit,bet_fit,gam_fit,m_fit,.25*N(1)/a./rho.rhoUp)/a, ...
                            Z*obj.V_fit(alp_fit,bet_fit,gam_fit,m_fit,.25*N(2)/a./rho.rhoDw)/a];

    Vx.getPotential(rho,V,true);
    R                   =   isa(V,'QMol_DFT_Vks')                                       && ... potential object
                            max(abs( (V.potentialUp   - V_th(:,1))./V_th(:,1) )) < 1e-10 && ... potential values
                            max(abs( (V.potentialDown - V_th(:,2))./V_th(:,2) )) < 1e-10;
    obj.showResult('getPotential (ADSIC and add to existing potential)' ,R);

    % Potential derivative -- no SIC ------
    DFT.set('SIC','none');              DFT.initialize;

    V_D                 =   Vx.getPotential(rho_D);
    DV_D                =   Vx.getPotentialDerivative(1,rho_D);
    R                   =   isa(DV_D,'QMol_DFT_Vks_grad') && isa(V_D,'QMol_DFT_Vks')                                && ... potential objects
                            max(abs(DV_D.potentialGradientUp   - ifft(DFT.disc.D.*fft(V_D.potentialUp)))) < 1e-10   && ...
                            max(abs(DV_D.potentialGradientDown - ifft(DFT.disc.D.*fft(V_D.potentialDown)))) < 1e-10;
    obj.showResult('getPotentialDerivative (no SIC and new)' ,R);  

    % ADSIC
    DFT.set('SIC','averageDensity');    DFT.initialize; 

    V_D                 =   Vx.getPotential(rho_D,V_D,true);
    DV_D                =   Vx.getPotentialDerivative(1,rho_D,DV_D,true);
    R                   =   isa(DV_D,'QMol_DFT_Vks_grad') && isa(V_D,'QMol_DFT_Vks')                                && ... potential objects
                            all(abs(DV_D.potentialGradientUp - ifft(DFT.disc.D.*fft(V_D.potentialUp)))    < 1e-10)  && ...
                            all(abs(DV_D.potentialGradientDown - ifft(DFT.disc.D.*fft(V_D.potentialDown))) < 1e-10);
    obj.showResult('getPotentialDerivative (ADSIC and add to existing potential)' ,R);

    % Energy -- no SIC ---------------------
    DFT.set('SIC','none');              DFT.initialize;

    E_th                =  [Z*obj.eps_fit(alp_fit,bet_fit,gam_fit,m_fit,.25/a./rho_D.rhoUp)/a, ...
                            Z*obj.eps_fit(alp_fit,bet_fit,gam_fit,m_fit,.25/a./rho_D.rhoDw)/a];
    E                   =   Vx.getEnergy(rho_D);
    R                   =   abs(E - sum(rho_D.rhoUp.*E_th(:,1) + rho_D.rhoDw.*E_th(:,2))*dx ) < 1e-10;
    obj.showResult('getEnergy (no SIC)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity');    DFT.initialize;

    N                   =   [sum(DFT.occupation{1}), sum(DFT.occupation{2})];
    
    E_th                =  [Z*obj.eps_fit(alp_fit,bet_fit,gam_fit,m_fit,.25     /a./rho_D.rhoUp)/a, ...
                            Z*obj.eps_fit(alp_fit,bet_fit,gam_fit,m_fit,.25     /a./rho_D.rhoDw)/a] ...
                          -[Z*obj.eps_fit(alp_fit,bet_fit,gam_fit,m_fit,.25*N(1)/a./rho_D.rhoUp)/a, ...
                            Z*obj.eps_fit(alp_fit,bet_fit,gam_fit,m_fit,.25*N(2)/a./rho_D.rhoDw)/a];
    E                   =   Vx.getEnergy(rho_D);
    R                   =   abs(E - sum(rho_D.rhoUp.*E_th(:,1) + rho_D.rhoDw.*E_th(:,2))*dx ) < 1e-10;
    obj.showResult('getEnergy (ADSIC)' ,R);

end
end
%% LDA exchange functional components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=private)
function E = eps_fit(~,alp,bet,gam,m,r)
%eps_fit fitted exchange energy per particle
    E                   =   -.5 * (1+alp*r)./(bet*r+2*alp*m*r.^2) .* log(1 + bet*r + gam*r.^m);
end
function DE = D_eps_fit(~,alp,bet,gam,m,r)
%D_eps_fit derivative of the fitted exchange energy per particle
    DE                  =   -.5*( alp./(bet*r+2*alp*m*r.^2) - (1+alp*r).*(bet+4*alp*m*r)./(bet*r+2*alp*m*r.^2).^2 ) .* log(1 + bet*r + gam*r.^m) ...
                            -.5 * (1+alp*r)./(bet*r+2*alp*m*r.^2) .* (bet + m*gam*r.^(m-1)) ./ (1 + bet*r + gam*r.^m);
end
function V = V_fit(obj,alp,bet,gam,m,r)
%D_eps_fit exchange potential associated with the fitted energy per particle
    V                   =   obj.eps_fit(alp,bet,gam,m,r) - r .* obj.D_eps_fit(alp,bet,gam,m,r);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

