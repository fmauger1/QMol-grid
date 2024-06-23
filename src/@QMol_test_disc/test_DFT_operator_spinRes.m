function test_DFT_operator_spinRes(obj)
%test_DFT_operator_spinRes unit tests for DFT operators

    % Initialization
    obj.showSection('DFT operators and functionals (spin restricted)');

    Vext                =   QMol_DFT_Vext('atom', ...
                               {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5), ...
                                QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3), ...
                                QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)});
    Vh                  =   QMol_DFT_Vh_conv;

    DFT                 =   QMol_DFT_spinRes(...
                                'xspan',                            -15:.1:20,  ...
                                'occupation',                       [2 1 1],    ...
                                'externalPotential',                Vext,       ...
                                'HartreePotential',                 Vh,         ...
                                'exchangeCorrelationPotential',     []);
    DFT.initialize;

    A                   =   .1;

    % DFT_energyKinetic ===================================================
    x                   =   DFT.x(:);   dx  =   x(2)-x(1);
    s                   =   2;
    p                   =  [exp(-x.^2*.25/s^2),x.*exp(-x.^2*.25/s^2),x.^2.*exp(-x.^2*.25/s^2)];
    p(:,1)              =   p(:,1) / sqrt(sum(p(:,1).^2)*dx);
    p(:,2)              =   p(:,2) / sqrt(sum(p(:,2).^2)*dx);
    p(:,3)              =   p(:,3) / sqrt(sum(p(:,3).^2)*dx);

    % Length gauge
    DFT.KSO.KSO         =   p;
    E                   =   DFT.disc.DFT_energyKinetic(DFT.occ,DFT.KSO);
    
    E_test              =  [sum( p(:,1) .* real(ifft(DFT.disc.T.*fft(p(:,1)))) ), ...
                            sum( p(:,2) .* real(ifft(DFT.disc.T.*fft(p(:,2)))) ), ...
                            sum( p(:,3) .* real(ifft(DFT.disc.T.*fft(p(:,3)))) )];
    E_test              =   sum(DFT.occ.*E_test) * dx;
    obj.showResult('DFT_energyKinetic (real orbitals)',isreal(E) && abs(E-E_test) < 1e-10);

    DFT.KSO.KSO(:,1)    =   DFT.KSO.KSO(:,1) * exp(1i*.3);
    DFT.KSO.KSO(:,2)    =   DFT.KSO.KSO(:,2) * exp(1i*.1);
    DFT.KSO.KSO(:,3)    =   DFT.KSO.KSO(:,3) * exp(1i*.5);
    E                   =   DFT.disc.DFT_energyKinetic(DFT.occ,DFT.KSO);
    obj.showResult('DFT_energyKinetic (complex orbitals)',isreal(E) && abs(E-E_test) < 1e-10);

    % Velocity gauge
    DFT.disc.setTv(A);

    DFT.KSO.KSO         =   p;
    E                   =   DFT.disc.DFT_energyKinetic(DFT.occ,DFT.KSO);
    
    E_test              =   E_test + .5*A^2 * sum(DFT.occ);
    obj.showResult('DFT_energyKinetic (velocity gauge, real orbitals)',isreal(E) && abs(E-E_test) < 1e-10);

    DFT.KSO.KSO(:,1)    =   DFT.KSO.KSO(:,1) * exp(1i*.3);
    DFT.KSO.KSO(:,2)    =   DFT.KSO.KSO(:,2) * exp(1i*.1);
    DFT.KSO.KSO(:,3)    =   DFT.KSO.KSO(:,3) * exp(1i*.5);
    E                   =   DFT.disc.DFT_energyKinetic(DFT.occ,DFT.KSO);
    obj.showResult('DFT_energyKinetic (velocity gauge, complex orbitals)',isreal(E) && abs(E-E_test) < 1e-10);

    DFT.disc.setTv([]);

    % DFT_energyOrbital ===================================================
    x                   =   DFT.x(:);   dx  =   x(2)-x(1);
    s                   =   2;
    p                   =  [exp(-x.^2*.25/s^2),x.*exp(-x.^2*.25/s^2),x.^2.*exp(-x.^2*.25/s^2)];
    p(:,1)              =   p(:,1) / sqrt(sum(p(:,1).^2)*dx);
    p(:,2)              =   p(:,2) / sqrt(sum(p(:,2).^2)*dx);
    p(:,3)              =   p(:,3) / sqrt(sum(p(:,3).^2)*dx);

    V                   =   Vext.getPotential;

    % Length gauge
    DFT.KSO.KSO         =   p;
    [E,err]             =   DFT.disc.DFT_energyOrbital(V,DFT.KSO);
    
    Hp                  =  [real(ifft(DFT.disc.T.*fft(p(:,1)))) + V.V.*p(:,1), ...
                            real(ifft(DFT.disc.T.*fft(p(:,2)))) + V.V.*p(:,2), ...
                            real(ifft(DFT.disc.T.*fft(p(:,3)))) + V.V.*p(:,3)];
    E_test              =   [sum(Hp(:,1).*p(:,1)); sum(Hp(:,2).*p(:,2)); sum(Hp(:,3).*p(:,3))]*dx;
    err_test            =   sqrt([sum(abs(Hp(:,1)-E_test(1)*p(:,1)).^2); ...
                                  sum(abs(Hp(:,2)-E_test(2)*p(:,2)).^2); ...
                                  sum(abs(Hp(:,3)-E_test(3)*p(:,3)).^2)]*dx);
    obj.showResult('DFT_energyOrbital (real orbitals)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));
    
    DFT.KSO.KSO(:,1)    =   DFT.KSO.KSO(:,1) * exp(1i*.3);
    DFT.KSO.KSO(:,2)    =   DFT.KSO.KSO(:,2) * exp(1i*.1);
    DFT.KSO.KSO(:,3)    =   DFT.KSO.KSO(:,3) * exp(1i*.5);
    [E,err]             =   DFT.disc.DFT_energyOrbital(V,DFT.KSO);
    
    obj.showResult('DFT_energyOrbital (complex orbitals)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    % Velocity gauge
    DFT.disc.setTv(A);

    DFT.KSO.KSO         =   p;
    [E,err]             =   DFT.disc.DFT_energyOrbital(V,DFT.KSO);
    
    Hpv                 =  [ifft(DFT.disc.Tv.*fft(p(:,1))) + V.V.*p(:,1), ...
                            ifft(DFT.disc.Tv.*fft(p(:,2))) + V.V.*p(:,2), ...
                            ifft(DFT.disc.Tv.*fft(p(:,3))) + V.V.*p(:,3)];
    E_test              =   E_test + .5*A^2;
    err_test            =   sqrt([sum(abs(Hpv(:,1)-E_test(1)*p(:,1)).^2); ...
                                  sum(abs(Hpv(:,2)-E_test(2)*p(:,2)).^2); ...
                                  sum(abs(Hpv(:,3)-E_test(3)*p(:,3)).^2)]*dx);
    obj.showResult('DFT_energyOrbital (velocity gauge, real orbitals)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));
    
    DFT.KSO.KSO(:,1)    =   DFT.KSO.KSO(:,1) * exp(1i*.3);
    DFT.KSO.KSO(:,2)    =   DFT.KSO.KSO(:,2) * exp(1i*.1);
    DFT.KSO.KSO(:,3)    =   DFT.KSO.KSO(:,3) * exp(1i*.5);
    [E,err]             =   DFT.disc.DFT_energyOrbital(V,DFT.KSO);
    
    obj.showResult('DFT_energyOrbital (velocity gauge, complex orbitals)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    DFT.disc.setTv([]);

    % DFT_operatorHamiltonian =============================================
    x0                  =   5;
    s                   =   2;
    p                   =   exp(-(DFT.x(:)-x0).^2 * .5/s^2);
    D2p                 =  -.5 * ((DFT.x(:)-x0).^2/s^4 - 1/s^2) .* p;

    A                   =   .1;
    D2pv                =   D2p + .5*A^2*p + A*1i/s^2 * (DFT.x(:)-x0) .*p;

    V                   =   Vext.getPotential;
    Vp                  =   V.V .* p;

    % Length gauge
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p);
    R                   =   max(abs(Hp-D2p-Vp)) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (no symmetry)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,1);
    R                   =   max(abs(Hp-.5*(D2p+Vp+flip(D2p+Vp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (symmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,-1);
    R                   =   max(abs(Hp-.5*(D2p+Vp-flip(D2p+Vp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (antisymmetrized)' ,R);

    % Velocity gauge
    DFT.disc.setTv(A);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p);
    R                   =   max(abs(Hp-D2pv-Vp)) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (velocity gauge, no symmetry)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,1);
    R                   =   max(abs(Hp-.5*(D2pv+Vp+flip(D2pv+Vp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (velocity gauge, symmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,-1);
    R                   =   max(abs(Hp-.5*(D2pv+Vp-flip(D2pv+Vp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (velocity gauge, antisymmetrized)' ,R);

    DFT.disc.setTv([]);

    % DFT_operatorKinetic =================================================
    x0                  =   5;
    s                   =   2;
    p                   =   exp(-(DFT.x(:)-x0).^2 * .5/s^2);
    D2p                 =  -.5 * ((DFT.x(:)-x0).^2/s^4 - 1/s^2) .* p;

    A                   =   .1;
    D2pv                =   D2p + .5*A^2*p + A*1i/s^2 * (DFT.x(:)-x0) .*p;

    % Length gauge
    Hp                  =   DFT.disc.DFT_operatorKinetic(p);
    R                   =   max(abs(Hp-D2p)) < 1e-10;
    obj.showResult('DFT_operatorKinetic (no symmetry)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorKinetic(p,1);
    R                   =   max(abs(Hp-.5*(D2p+flip(D2p)))) < 1e-10;
    obj.showResult('DFT_operatorKinetic (symmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorKinetic(p,-1);
    R                   =   max(abs(Hp-.5*(D2p-flip(D2p)))) < 1e-10;
    obj.showResult('DFT_operatorKinetic (antisymmetrized)' ,R);

    % Velocity gauge
    DFT.disc.setTv(A);

    Hp                  =   DFT.disc.DFT_operatorKinetic(p);
    R                   =   max(abs(Hp-D2pv)) < 1e-10;
    obj.showResult('DFT_operatorKinetic (velocity gauge, no symmetry)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorKinetic(p,1);
    R                   =   max(abs(Hp-.5*(D2pv+flip(D2pv)))) < 1e-10;
    obj.showResult('DFT_operatorKinetic (velocity gauge, symmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorKinetic(p,-1);
    R                   =   max(abs(Hp-.5*(D2pv-flip(D2pv)))) < 1e-10;
    obj.showResult('DFT_operatorKinetic (velocity gauge, antisymmetrized)' ,R);

    DFT.disc.setTv([]);

    % DFT_operatorPotential ===============================================
    x0                  =   5;
    s                   =   2;
    p                   =   exp(-(DFT.x(:)-x0).^2 * .5/s^2);

    V                   =   Vext.getPotential;
    Vp                  =   V.V .* p;

    % (same for length and velocity gauges)
    Hp                  =   DFT.disc.DFT_operatorPotential(V,p);
    R                   =   max(abs(Hp-Vp)) < 1e-10;
    obj.showResult('DFT_operatorPotential (no symmetry)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorPotential(V,p,1);
    R                   =   max(abs(Hp-.5*(Vp+flip(Vp)))) < 1e-10;
    obj.showResult('DFT_operatorPotential (symmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorPotential(V,p,-1);
    R                   =   max(abs(Hp-.5*(Vp-flip(Vp)))) < 1e-10;
    obj.showResult('DFT_operatorPotential (antisymmetrized)' ,R);

end