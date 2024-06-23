classdef QMol_test_DFT_Vc_LDA_soft < QMol_test
%QMol_test_DFT_Vc_LDA_soft unit test for the QMol_DFT_Vc_LDA_soft class

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
    fprintf('  * QMol_test_DFT_Vc_LDA_soft\n'); 
    QMol_test_DFT_Vc_LDA_soft.version;
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
    
    A0  =   18.40;  C0  =   7.501;  D0  =   0.10185;    E0  =   0.012827;
    A1  =   5.24;   C1  =   1.568;  D1  =   0.1286;     E1  =   0.00320;
    
    alpha0  =   1.511;  beta0   =   0.258;      m0  =   4.424;
    alpha1  =   0.0538; beta1   =   1.56e-5;    m1  =   2.958;
    
    eps_c               =   @(A,C,D,E,alp,bet,m,r) -.5 * (r+E*r.^2)./(A+C*r.^2+D*r.^3) .* log(1+alp*r+bet*r.^m);
    D_eps_c             =   @(A,C,D,E,alp,bet,m,r) -.5 * (A+2*A*E*r-C*r.^2-2*D*r.^3-D*E*r.^4)./(A+C*r.^2+D*r.^3).^2 .* log(1+alp*r+bet*r.^m) ...
                                                   -.5 * (r+E*r.^2)./(A+C*r.^2+D*r.^3) .* (alp+bet*m*r.^(m-1))./(1+alp*r+bet*r.^m);
                                               
    r_s                 =   @(r1,r2) .5 ./ (r1+r2);
    xi                  =   @(r1,r2) (r1-r2) ./ (r1+r2);
    eps_c_tot           =   @(r1,r2) eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ...
                                    +xi(r1,r2).^2 .*( ...
                                        eps_c(A1,C1,D1,E1,alpha1,beta1,m1,r_s(r1,r2)) ...
                                       -eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) );
    V_c_up              =   @(r1,r2) eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ...
                                    +xi(r1,r2).*(2-xi(r1,r2)) .* (...
                                     	eps_c(A1,C1,D1,E1,alpha1,beta1,m1,r_s(r1,r2)) ...
                                       -eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ) ...
                                    -r_s(r1,r2) .* (...
                                        D_eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ...
                                       +xi(r1,r2).^2 .* (...
                                            D_eps_c(A1,C1,D1,E1,alpha1,beta1,m1,r_s(r1,r2)) ...
                                           -D_eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ) );

    Vc                  =   QMol_DFT_Vc_LDA_soft;
    DFT                 =   QMol_DFT_spinRes(...
                                'xspan',                -15:.3:20,          ...
                                'occupation',           [5.2 3.4],          ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vc);
    DFT.initialize;
    rho                 =   DFT.disc.DFT_allocateDensity;
    rho.set('rho',10.^linspace(-2.99,2.02,numel(DFT.x)).');

    rho_D               =   DFT.disc.DFT_allocateDensity;
    rho_D.set('rho',exp(-.5*DFT.x(:).^2));

    dx                  =   DFT.x(2)-DFT.x(1);

    % Energy per particle derivative -------
    rs                  =   .5 ./ rho.rho;
    dr                  =   1e-4;
    
    DE                  =   (eps_c(A0,C0,D0,E0,alpha0,beta0,m0,rs*(1+.5*dr))-eps_c(A0,C0,D0,E0,alpha0,beta0,m0,rs*(1-.5*dr)))./rs/dr;
    R                   =   max((DE-D_eps_c(A0,C0,D0,E0,alpha0,beta0,m0,rs))./DE) < 1e-6;
    if ~R,  obj.showResult('run-time error: test results may not be accurate' ,R); end

    % Potential -- no SIC ------------------
    DFT.set('SIC','none');  DFT.initialize;

    V_th                =   V_c_up(.5*rho.rho,.5*rho.rho);
    V                   =   Vc.getPotential(rho);

    R                   =   isa(V,'QMol_DFT_Vks')                           && ... potential object
                            max(abs( (V.potential - V_th)./V_th )) < 1e-10;      % potential discretization
    obj.showResult('getPotential (no SIC and new)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity'); DFT.initialize;

    N                   =   sum(DFT.occupation);
    V_th                =   V_th + V_c_up(.5*rho.rho,.5*rho.rho) - V_c_up(rho.rho/N,0);

    Vc.getPotential(rho,V,true);
    R                   =   isa(V,'QMol_DFT_Vks')                           && ... potential object
                            max(abs( (V.potential - V_th)./V_th )) < 1e-10;       % potential discretization
    obj.showResult('getPotential (ADSIC and add to existing potential)' ,R);

    % Potential derivative -- no SIC ------
    DFT.set('SIC','none');           DFT.initialize;

    V_D                 =   Vc.getPotential(rho_D);
    DV_D                =   Vc.getPotentialDerivative(1,rho_D);
    R                   =   isa(DV_D,'QMol_DFT_Vks_grad') && isa(V_D,'QMol_DFT_Vks') && ... potential objects
                            max(abs(DV_D.potentialGradient - ifft(DFT.disc.D.*fft(V_D.potential)))) < 1e-10;
    obj.showResult('getPotentialDerivative (no SIC and new)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity');    DFT.initialize;

    V_D                 =   Vc.getPotential(rho_D,V_D,true);
    DV_D                =   Vc.getPotentialDerivative(1,rho_D,DV_D,true);
    R                   =   isa(DV_D,'QMol_DFT_Vks_grad') && isa(V_D,'QMol_DFT_Vks') && ... potential objects
                            max(abs(DV_D.potentialGradient - ifft(DFT.disc.D.*fft(V_D.potential)))) < 1e-10;
    obj.showResult('getPotentialDerivative (ADSIC and add to existing potential)' ,R);

    % Energy -- no SIC ---------------------
    DFT.set('SIC','none');              DFT.initialize;

    E_th                =   eps_c_tot(.5*rho_D.rho,.5*rho_D.rho);
    E                   =   Vc.getEnergy(rho_D);
    R                   =   abs(E - sum(rho_D.rho.*E_th,'omitnan')*dx ) < 1e-10;
    obj.showResult('getEnergy (no SIC)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity');    DFT.initialize;

    N                   =   sum(DFT.occupation);
    E_th1               =   eps_c_tot(.5*rho_D.rho,.5*rho_D.rho);
    E_th2               =   eps_c_tot(rho_D.rho/N,0);

    E                   =   Vc.getEnergy(rho_D);
    R                   =   abs(E - sum(rho_D.rho.*(E_th1-E_th2),'omitnan')*dx ) < 1e-10;
    obj.showResult('getEnergy (ADSIC)' ,R);

end
function runSpinPolarized(obj) %===========================================
%runSpinRestricted unit tests for spin-restricted DFT
    
    % Initialization
    obj.showSection('Spin-polarized DFT');
    
    A0  =   18.40;  C0  =   7.501;  D0  =   0.10185;    E0  =   0.012827;
    A1  =   5.24;   C1  =   1.568;  D1  =   0.1286;     E1  =   0.00320;
    
    alpha0  =   1.511;  beta0   =   0.258;      m0  =   4.424;
    alpha1  =   0.0538; beta1   =   1.56e-5;    m1  =   2.958;
    
    eps_c               =   @(A,C,D,E,alp,bet,m,r) -.5 * (r+E*r.^2)./(A+C*r.^2+D*r.^3) .* log(1+alp*r+bet*r.^m);
    D_eps_c             =   @(A,C,D,E,alp,bet,m,r) -.5 * (A+2*A*E*r-C*r.^2-2*D*r.^3-D*E*r.^4)./(A+C*r.^2+D*r.^3).^2 .* log(1+alp*r+bet*r.^m) ...
                                                   -.5 * (r+E*r.^2)./(A+C*r.^2+D*r.^3) .* (alp+bet*m*r.^(m-1))./(1+alp*r+bet*r.^m);
                                               
    r_s                 =   @(r1,r2) .5 ./ (r1+r2);
    xi                  =   @(r1,r2) (r1-r2) ./ (r1+r2);
    eps_c_tot           =   @(r1,r2) eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ...
                                    +xi(r1,r2).^2 .*( ...
                                        eps_c(A1,C1,D1,E1,alpha1,beta1,m1,r_s(r1,r2)) ...
                                       -eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) );
    V_c_up              =   @(r1,r2) eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ...
                                    +xi(r1,r2).*(2-xi(r1,r2)) .* (...
                                     	eps_c(A1,C1,D1,E1,alpha1,beta1,m1,r_s(r1,r2)) ...
                                       -eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ) ...
                                    -r_s(r1,r2) .* (...
                                        D_eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ...
                                       +xi(r1,r2).^2 .* (...
                                            D_eps_c(A1,C1,D1,E1,alpha1,beta1,m1,r_s(r1,r2)) ...
                                           -D_eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ) );
    V_c_down            =   @(r1,r2) eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ...
                                    -xi(r1,r2).*(2+xi(r1,r2)) .* (...
                                     	eps_c(A1,C1,D1,E1,alpha1,beta1,m1,r_s(r1,r2)) ...
                                       -eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ) ...
                                    -r_s(r1,r2) .* (...
                                        D_eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ...
                                       +xi(r1,r2).^2 .* (...
                                            D_eps_c(A1,C1,D1,E1,alpha1,beta1,m1,r_s(r1,r2)) ...
                                           -D_eps_c(A0,C0,D0,E0,alpha0,beta0,m0,r_s(r1,r2)) ) );

    Vc                  =   QMol_DFT_Vc_LDA_soft;
    DFT                 =   QMol_DFT_spinPol(...
                                'xspan',                -15:.3:20,          ...
                                'occupation',           {[5 2.4],[2 2 .1]}, ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vc);
    DFT.initialize;
    rho                 =   DFT.disc.DFT_allocateDensity;
    rho.set('rhoUp',10.^linspace(-3.01,2.03,numel(DFT.x)).');
    rho.set('rhoDw',10.^linspace(-2.52,1.87,numel(DFT.x)).');

    rho_D               =   DFT.disc.DFT_allocateDensity;
    rho_D.set('rhoUp',exp(-.5*DFT.x(:).^2),'rhoDw',exp(-.4*DFT.x(:).^2));

    dx                  =   DFT.x(2)-DFT.x(1);

    % Energy per particle derivative -------
    rs                  =   .5 ./ rho.rhoUp;
    dr                  =   1e-4;
    
    DE                  =   (eps_c(A0,C0,D0,E0,alpha0,beta0,m0,rs*(1+.5*dr))-eps_c(A0,C0,D0,E0,alpha0,beta0,m0,rs*(1-.5*dr)))./rs/dr;
    R                   =   max((DE-D_eps_c(A0,C0,D0,E0,alpha0,beta0,m0,rs))./DE) < 1e-6;
    if ~R,  obj.showResult('run-time error: test results may not be accurate' ,R); end


    % Poptential -- no SIC -----------------
    DFT.set('SIC','none');  DFT.initialize;

    V_th                =  [V_c_up(rho.rhoUp,rho.rhoDw), V_c_down(rho.rhoUp,rho.rhoDw)];

    V                   =   Vc.getPotential(rho);
    R                   =   isa(V,'QMol_DFT_Vks')                                       && ... potential object
                            max(abs( (V.potentialUp   - V_th(:,1))./V_th(:,1) )) < 1e-10 && ... potential values
                            max(abs( (V.potentialDown - V_th(:,2))./V_th(:,2) )) < 1e-10;
    obj.showResult('getPotential (no SIC and new)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity'); DFT.initialize;

    N                   =   [sum(DFT.occupation{1}), sum(DFT.occupation{2})];
    V_th                =   V_th + ...
                           [V_c_up(  rho.rhoUp,rho.rhoDw)-V_c_up(  rho.rhoUp/N(1),0), ...
                            V_c_down(rho.rhoUp,rho.rhoDw)-V_c_down(0,rho.rhoDw/N(2))];

    Vc.getPotential(rho,V,true);
    R                   =   isa(V,'QMol_DFT_Vks')                                       && ... potential object
                            max(abs( (V.potentialUp   - V_th(:,1))./V_th(:,1) )) < 1e-10 && ... potential values
                            max(abs( (V.potentialDown - V_th(:,2))./V_th(:,2) )) < 1e-10;
    obj.showResult('getPotential (ADSIC and add to existing potential)' ,R);

    % Potential derivative -- no SIC ------
    DFT.set('SIC','none');              DFT.initialize;

    V_D                 =   Vc.getPotential(rho_D);
    DV_D                =   Vc.getPotentialDerivative(1,rho_D);
    R                   =   isa(DV_D,'QMol_DFT_Vks_grad') && isa(V_D,'QMol_DFT_Vks')                               && ... potential objects
                            max(abs(DV_D.potentialGradientUp   - ifft(DFT.disc.D.*fft(V_D.potentialUp)))) < 1e-10  && ...
                            max(abs(DV_D.potentialGradientDown - ifft(DFT.disc.D.*fft(V_D.potentialDown)))) < 1e-10;
    obj.showResult('getPotentialDerivative (no SIC and new)' ,R);
    
    % ADSIC
    DFT.set('SIC','averageDensity');    DFT.initialize;

    V_D                 =   Vc.getPotential(rho_D,V_D,true);
    DV_D                =   Vc.getPotentialDerivative(1,rho_D,DV_D,true);
    R                   =   isa(DV_D,'QMol_DFT_Vks_grad') && isa(V_D,'QMol_DFT_Vks')                               && ... potential objects
                            all(abs(DV_D.potentialGradientUp - ifft(DFT.disc.D.*fft(V_D.potentialUp)))    < 1e-10) && ...
                            all(abs(DV_D.potentialGradientDown - ifft(DFT.disc.D.*fft(V_D.potentialDown))) < 1e-10);
    obj.showResult('getPotentialDerivative (ADSIC and add to exist. pot. grad.)' ,R);

    % Energy -- no SIC ---------------------
    DFT.set('SIC','none');              DFT.initialize;

    E_th                =   eps_c_tot(rho_D.rhoUp,rho_D.rhoDw);
    E                   =   Vc.getEnergy(rho_D);
    R                   =   abs(E - sum( (rho_D.rhoUp+rho_D.rhoDw).*E_th,'omitnan' )*dx ) < 1e-10;
    obj.showResult('getEnergy (no SIC)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity');    DFT.initialize;

    E_th                =   eps_c_tot(rho_D.rhoUp,rho_D.rhoDw);
    E_u                 =   eps_c_tot(rho_D.rhoUp/N(1),0);
    E_d                 =   eps_c_tot(0,rho_D.rhoDw/N(2));

    E                   =   Vc.getEnergy(rho_D);
    R                   =   abs(E - sum( (rho_D.rhoUp+rho_D.rhoDw).*E_th - rho_D.rhoUp.*E_u - rho_D.rhoDw.*E_d,'omitnan')*dx ) < 1e-10;
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

