function test_SE_operator(obj)
%test_SE_operator unit tests for Schrodinger-equation operators

    % Initialization
    obj.showSection('Schrodinger-equation operators and functionals');

    Vse                 =   QMol_SE_V('atom', ...
                               {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5), ...
                                QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3), ...
                                QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)});

    SE                  =   QMol_SE(...
                                'xspan',                       -15:.1:20,  ...
                                'potential',                    Vse);
    SE.initialize;

    A                   =   .1;

    % SE_energyKinetic ====================================================
    x                   =   SE.x(:);   dx  =   x(2)-x(1);
    s                   =   2;
    p                   =  [exp(-x.^2*.25/s^2),x.*exp(-x.^2*.25/s^2),x.^2.*exp(-x.^2*.25/s^2)];
    p(:,1)              =   p(:,1) / sqrt(sum(p(:,1).^2)*dx);
    p(:,2)              =   p(:,2) / sqrt(sum(p(:,2).^2)*dx);
    p(:,3)              =   p(:,3) / sqrt(sum(p(:,3).^2)*dx);

    % Length gauge
    SE.wfcn.wfcn        =   p;
    E                   =   SE.disc.SE_energyKinetic(SE.wfcn);
    
    E_test              =  [sum( p(:,1) .* real(ifft(SE.disc.T.*fft(p(:,1)))) ), ...
                            sum( p(:,2) .* real(ifft(SE.disc.T.*fft(p(:,2)))) ), ...
                            sum( p(:,3) .* real(ifft(SE.disc.T.*fft(p(:,3)))) )];
    E_test              =   sum(E_test) * dx;
    obj.showResult('SE_energyKinetic (real orbitals)',isreal(E) && abs(E-E_test) < 1e-10);

    SE.wfcn.wfcn(:,1)   =   SE.wfcn.wfcn(:,1) * exp(1i*.3);
    SE.wfcn.wfcn(:,2)   =   SE.wfcn.wfcn(:,2) * exp(1i*.1);
    SE.wfcn.wfcn(:,3)   =   SE.wfcn.wfcn(:,3) * exp(1i*.5);
    E                   =   SE.disc.SE_energyKinetic(SE.wfcn);
    obj.showResult('SE_energyKinetic (complex orbitals)',isreal(E) && abs(E-E_test) < 1e-10);

    % Velocity gauge
    SE.disc.setTv(A);

    SE.wfcn.wfcn        =   p;
    E                   =   SE.disc.SE_energyKinetic(SE.wfcn);
    
    E_test              =   E_test + .5*A^2 * SE.numberWaveFunction;
    obj.showResult('SE_energyKinetic (velocity gauge, real orbitals)',isreal(E) && abs(E-E_test) < 1e-10);

    SE.wfcn.wfcn(:,1)   =   SE.wfcn.wfcn(:,1) * exp(1i*.3);
    SE.wfcn.wfcn(:,2)   =   SE.wfcn.wfcn(:,2) * exp(1i*.1);
    SE.wfcn.wfcn(:,3)   =   SE.wfcn.wfcn(:,3) * exp(1i*.5);
    E                   =   SE.disc.SE_energyKinetic(SE.wfcn);
    obj.showResult('SE_energyKinetic (velocity gauge, complex orbitals)',isreal(E) && abs(E-E_test) < 1e-10);

    SE.disc.setTv([]);

    % SE_energyWaveFunction ===============================================
    x                   =   SE.x(:);   dx  =   x(2)-x(1);
    s                   =   2;
    p                   =  [exp(-x.^2*.25/s^2),x.*exp(-x.^2*.25/s^2),x.^2.*exp(-x.^2*.25/s^2)];
    p(:,1)              =   p(:,1) / sqrt(sum(p(:,1).^2)*dx);
    p(:,2)              =   p(:,2) / sqrt(sum(p(:,2).^2)*dx);
    p(:,3)              =   p(:,3) / sqrt(sum(p(:,3).^2)*dx);

    % Length gauge
    SE.wfcn.wfcn        =   p;
    [E,err]             =   SE.disc.SE_energyWaveFunction(Vse,SE.wfcn);
    
    Hp                  =  [real(ifft(SE.disc.T.*fft(p(:,1)))) + Vse.V.*p(:,1), ...
                            real(ifft(SE.disc.T.*fft(p(:,2)))) + Vse.V.*p(:,2), ...
                            real(ifft(SE.disc.T.*fft(p(:,3)))) + Vse.V.*p(:,3)];
    E_test              =   [sum(Hp(:,1).*p(:,1)); sum(Hp(:,2).*p(:,2)); sum(Hp(:,3).*p(:,3))]*dx;
    err_test            =   sqrt([sum(abs(Hp(:,1)-E_test(1)*p(:,1)).^2); ...
                                  sum(abs(Hp(:,2)-E_test(2)*p(:,2)).^2); ...
                                  sum(abs(Hp(:,3)-E_test(3)*p(:,3)).^2)]*dx);
    obj.showResult('SE_energyWaveFunction (real orbitals)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));
    
    SE.wfcn.wfcn(:,1)   =   SE.wfcn.wfcn(:,1) * exp(1i*.3);
    SE.wfcn.wfcn(:,2)   =   SE.wfcn.wfcn(:,2) * exp(1i*.1);
    SE.wfcn.wfcn(:,3)   =   SE.wfcn.wfcn(:,3) * exp(1i*.5);
    [E,err]             =   SE.disc.SE_energyWaveFunction(Vse,SE.wfcn);
    
    obj.showResult('SE_energyWaveFunction (complex orbitals)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    % Velocity gauge
    SE.disc.setTv(A);

    SE.wfcn.wfcn        =   p;
    [E,err]             =   SE.disc.SE_energyWaveFunction(Vse,SE.wfcn);
    
    Hpv                 =  [ifft(SE.disc.Tv.*fft(p(:,1))) + Vse.V.*p(:,1), ...
                            ifft(SE.disc.Tv.*fft(p(:,2))) + Vse.V.*p(:,2), ...
                            ifft(SE.disc.Tv.*fft(p(:,3))) + Vse.V.*p(:,3)];
    E_test              =   E_test + .5*A^2;
    err_test            =   sqrt([sum(abs(Hpv(:,1)-E_test(1)*p(:,1)).^2); ...
                                  sum(abs(Hpv(:,2)-E_test(2)*p(:,2)).^2); ...
                                  sum(abs(Hpv(:,3)-E_test(3)*p(:,3)).^2)]*dx);
    obj.showResult('SE_energyWaveFunction (velocity gauge, real orbitals)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));
    
    SE.wfcn.wfcn(:,1)  =   SE.wfcn.wfcn(:,1) * exp(1i*.3);
    SE.wfcn.wfcn(:,2)  =   SE.wfcn.wfcn(:,2) * exp(1i*.1);
    SE.wfcn.wfcn(:,3)  =   SE.wfcn.wfcn(:,3) * exp(1i*.5);
    [E,err]            =   SE.disc.SE_energyWaveFunction(Vse,SE.wfcn);
    
    obj.showResult('SE_energyWaveFunction (velocity gauge, complex orbitals)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    SE.disc.setTv([]);

    % SE_operatorHamiltonian ==============================================
    x0                  =   5;
    s                   =   2;
    p                   =   exp(-(SE.x(:)-x0).^2 * .5/s^2);
    D2p                 =  -.5 * ((SE.x(:)-x0).^2/s^4 - 1/s^2) .* p;

    A                   =   .1;
    D2pv                =   D2p + .5*A^2*p + A*1i/s^2 * (SE.x(:)-x0) .*p;

    Vp                  =   Vse.V .* p;

    % Length gauge
    Hp                  =   SE.disc.SE_operatorHamiltonian(Vse,p);
    R                   =   max(abs(Hp-D2p-Vp)) < 1e-10;
    obj.showResult('SE_operatorHamiltonian (no symmetry)' ,R);
    
    Hp                  =   SE.disc.SE_operatorHamiltonian(Vse,p,1);
    R                   =   max(abs(Hp-.5*(D2p+Vp+flip(D2p+Vp)))) < 1e-10;
    obj.showResult('SE_operatorHamiltonian (symmetrized)' ,R);
    
    Hp                  =   SE.disc.SE_operatorHamiltonian(Vse,p,-1);
    R                   =   max(abs(Hp-.5*(D2p+Vp-flip(D2p+Vp)))) < 1e-10;
    obj.showResult('SE_operatorHamiltonian (antisymmetrized)' ,R);

    % Velocity gauge
    SE.disc.setTv(A);

    Hp                  =   SE.disc.SE_operatorHamiltonian(Vse,p);
    R                   =   max(abs(Hp-D2pv-Vp)) < 1e-10;
    obj.showResult('SE_operatorHamiltonian (velocity gauge, no symmetry)' ,R);
    
    Hp                  =   SE.disc.SE_operatorHamiltonian(Vse,p,1);
    R                   =   max(abs(Hp-.5*(D2pv+Vp+flip(D2pv+Vp)))) < 1e-10;
    obj.showResult('SE_operatorHamiltonian (velocity gauge, symmetrized)' ,R);
    
    Hp                  =   SE.disc.SE_operatorHamiltonian(Vse,p,-1);
    R                   =   max(abs(Hp-.5*(D2pv+Vp-flip(D2pv+Vp)))) < 1e-10;
    obj.showResult('SE_operatorHamiltonian (velocity gauge, antisymmetrized)' ,R);

    SE.disc.setTv([]);

    % DFT_operatorKinetic =================================================
    x0                  =   5;
    s                   =   2;
    p                   =   exp(-(SE.x(:)-x0).^2 * .5/s^2);
    D2p                 =  -.5 * ((SE.x(:)-x0).^2/s^4 - 1/s^2) .* p;

    A                   =   .1;
    D2pv                =   D2p + .5*A^2*p + A*1i/s^2 * (SE.x(:)-x0) .*p;

    % Length gauge
    Hp                  =   SE.disc.SE_operatorKinetic(p);
    R                   =   max(abs(Hp-D2p)) < 1e-10;
    obj.showResult('SE_operatorKinetic (no symmetry)' ,R);
    
    Hp                  =   SE.disc.SE_operatorKinetic(p,1);
    R                   =   max(abs(Hp-.5*(D2p+flip(D2p)))) < 1e-10;
    obj.showResult('SE_operatorKinetic (symmetrized)' ,R);
    
    Hp                  =   SE.disc.SE_operatorKinetic(p,-1);
    R                   =   max(abs(Hp-.5*(D2p-flip(D2p)))) < 1e-10;
    obj.showResult('SE_operatorKinetic (antisymmetrized)' ,R);

    % Velocity gauge
    SE.disc.setTv(A);

    Hp                  =   SE.disc.SE_operatorKinetic(p);
    R                   =   max(abs(Hp-D2pv)) < 1e-10;
    obj.showResult('SE_operatorKinetic (velocity gauge, no symmetry)' ,R);
    
    Hp                  =   SE.disc.SE_operatorKinetic(p,1);
    R                   =   max(abs(Hp-.5*(D2pv+flip(D2pv)))) < 1e-10;
    obj.showResult('SE_operatorKinetic (velocity gauge, symmetrized)' ,R);
    
    Hp                  =   SE.disc.SE_operatorKinetic(p,-1);
    R                   =   max(abs(Hp-.5*(D2pv-flip(D2pv)))) < 1e-10;
    obj.showResult('SE_operatorKinetic (velocity gauge, antisymmetrized)' ,R);

    SE.disc.setTv([]);

    % DFT_operatorPotential ===============================================
    x0                  =   5;
    s                   =   2;
    p                   =   exp(-(SE.x(:)-x0).^2 * .5/s^2);

    Vp                  =   Vse.V .* p;

    % (same for length and velocity gauges)
    Hp                  =   SE.disc.SE_operatorPotential(Vse,p);
    R                   =   max(abs(Hp-Vp)) < 1e-10;
    obj.showResult('SE_operatorPotential (no symmetry)' ,R);
    
    Hp                  =   SE.disc.SE_operatorPotential(Vse,p,1);
    R                   =   max(abs(Hp-.5*(Vp+flip(Vp)))) < 1e-10;
    obj.showResult('SE_operatorPotential (symmetrized)' ,R);
    
    Hp                  =   SE.disc.SE_operatorPotential(Vse,p,-1);
    R                   =   max(abs(Hp-.5*(Vp-flip(Vp)))) < 1e-10;
    obj.showResult('SE_operatorPotential (antisymmetrized)' ,R);

end