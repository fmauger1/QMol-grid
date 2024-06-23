function test_DFT_operator_spinPol(obj)
%test_operator unit tests for the operator overload in the class

    % Initialization
    obj.showSection('DFT operators and functionals (spin polarized)');

    Vext                =   QMol_DFT_Vext('atom', ...
                               {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5), ...
                                QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3), ...
                                QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)});
    Vh                  =   QMol_DFT_Vh_conv;

    DFT                 =   QMol_DFT_spinPol(...
                                'xspan',                            -15:.1:20,  ...
                                'occupation',                       {[1 .5 .5],[1 .1]},    ...
                                'externalPotential',                Vext,       ...
                                'HartreePotential',                 Vh,         ...
                                'exchangeCorrelationPotential',     []);
    DFT.initialize;

    A                   =   .1;

    % DFT_energyKinetic ===================================================
    x                   =   DFT.x(:);   dx  =   x(2)-x(1);
    s                   =  [2 1.5];
    p                   =  [exp(-x.^2*.25/s(1)^2),x.*exp(-x.^2*.25/s(1)^2),x.^2.*exp(-x.^2*.25/s(1)^2),...
                            exp(-x.^2*.25/s(2)^2),x.*exp(-x.^2*.25/s(2)^2),];
    p(:,1)              =   p(:,1) / sqrt(sum(p(:,1).^2)*dx);
    p(:,2)              =   p(:,2) / sqrt(sum(p(:,2).^2)*dx);
    p(:,3)              =   p(:,3) / sqrt(sum(p(:,3).^2)*dx);
    p(:,4)              =   p(:,4) / sqrt(sum(p(:,4).^2)*dx);
    p(:,5)              =   p(:,5) / sqrt(sum(p(:,5).^2)*dx);

    % Length gauge
    x                   =   DFT.x(:);   dx  =   x(2)-x(1);
    s                   =  [2 1.5];
    p                   =  [exp(-x.^2*.25/s(1)^2),x.*exp(-x.^2*.25/s(1)^2),x.^2.*exp(-x.^2*.25/s(1)^2),...
                            exp(-x.^2*.25/s(2)^2),x.*exp(-x.^2*.25/s(2)^2),];
    p(:,1)              =   p(:,1) / sqrt(sum(p(:,1).^2)*dx);
    p(:,2)              =   p(:,2) / sqrt(sum(p(:,2).^2)*dx);
    p(:,3)              =   p(:,3) / sqrt(sum(p(:,3).^2)*dx);
    p(:,4)              =   p(:,4) / sqrt(sum(p(:,4).^2)*dx);
    p(:,5)              =   p(:,5) / sqrt(sum(p(:,5).^2)*dx);

    DFT.KSO.KSOup       =   p(:,1:3);
    DFT.KSO.KSOdw       =   p(:,4:5);
    E                   =   DFT.disc.DFT_energyKinetic(DFT.occ,DFT.KSO);
    
    E_test              =  [sum( p(:,1) .* real(ifft(DFT.disc.T.*fft(p(:,1)))) ), ...
                            sum( p(:,2) .* real(ifft(DFT.disc.T.*fft(p(:,2)))) ), ...
                            sum( p(:,3) .* real(ifft(DFT.disc.T.*fft(p(:,3)))) ), ...
                            sum( p(:,4) .* real(ifft(DFT.disc.T.*fft(p(:,4)))) ), ...
                            sum( p(:,5) .* real(ifft(DFT.disc.T.*fft(p(:,5)))) )];

    E_test              =  [sum(DFT.occ{1}.*E_test(1:3)),sum(DFT.occ{2}.*E_test(4:5))] * dx;
    obj.showResult('DFT_energyKinetic (real orbitals)',isreal(E) && all(abs(E-E_test) < 1e-10));

    DFT.KSO.KSOup(:,1)  =   DFT.KSO.KSOup(:,1) * exp(1i*.3);
    DFT.KSO.KSOup(:,2)  =   DFT.KSO.KSOup(:,2) * exp(1i*.1);
    DFT.KSO.KSOup(:,3)  =   DFT.KSO.KSOup(:,3) * exp(1i*.5);
    DFT.KSO.KSOdw(:,1)  =   DFT.KSO.KSOdw(:,1) * exp(1i*.8);
    DFT.KSO.KSOdw(:,2)  =   DFT.KSO.KSOdw(:,2) * exp(1i*.7);
    E                   =   DFT.disc.DFT_energyKinetic(DFT.occ,DFT.KSO);
    obj.showResult('DFT_energyKinetic (complex orbitals)',isreal(E) && all(abs(E-E_test) < 1e-10));

    % Velocity gauge
    DFT.disc.setTv(A);


    DFT.KSO.KSOup       =   p(:,1:3);
    DFT.KSO.KSOdw       =   p(:,4:5);
    E                   =   DFT.disc.DFT_energyKinetic(DFT.occ,DFT.KSO);
    
    E_test              =   E_test + .5*A^2 * [sum(DFT.occ{1}), sum(DFT.occ{2})];
    obj.showResult('DFT_energyKinetic (velocity gauge, real orbitals)',isreal(E) && all(abs(E-E_test) < 1e-10));

    DFT.KSO.KSOup(:,1)  =   DFT.KSO.KSOup(:,1) * exp(1i*.3);
    DFT.KSO.KSOup(:,2)  =   DFT.KSO.KSOup(:,2) * exp(1i*.1);
    DFT.KSO.KSOup(:,3)  =   DFT.KSO.KSOup(:,3) * exp(1i*.5);
    DFT.KSO.KSOdw(:,1)  =   DFT.KSO.KSOdw(:,1) * exp(1i*.8);
    DFT.KSO.KSOdw(:,2)  =   DFT.KSO.KSOdw(:,2) * exp(1i*.7);
    E                   =   DFT.disc.DFT_energyKinetic(DFT.occ,DFT.KSO);
    obj.showResult('DFT_energyKinetic (velocity gauge, complex orbitals)',isreal(E) && all(abs(E-E_test) < 1e-10));

    DFT.disc.setTv([]);


    % DFT_energyOrbital ===================================================
    x                   =   DFT.x(:);   dx  =   x(2)-x(1);
    s                   =  [2 1.5];
    p                   =  [exp(-x.^2*.25/s(1)^2),x.*exp(-x.^2*.25/s(1)^2),x.^2.*exp(-x.^2*.25/s(1)^2),...
                            exp(-x.^2*.25/s(2)^2),x.*exp(-x.^2*.25/s(2)^2),];
    p(:,1)              =   p(:,1) / sqrt(sum(p(:,1).^2)*dx);
    p(:,2)              =   p(:,2) / sqrt(sum(p(:,2).^2)*dx);
    p(:,3)              =   p(:,3) / sqrt(sum(p(:,3).^2)*dx);
    p(:,4)              =   p(:,4) / sqrt(sum(p(:,4).^2)*dx);
    p(:,5)              =   p(:,5) / sqrt(sum(p(:,5).^2)*dx);

    V                   =   Vext.getPotential;

    % Length gauge
    DFT.KSO.KSOup       =   p(:,1:3);
    DFT.KSO.KSOdw       =   p(:,4:5);
    [E,err]             =   DFT.disc.DFT_energyOrbital(V,DFT.KSO);
    
    Hp                  =  [real(ifft(DFT.disc.T.*fft(p(:,1)))) + V.Vup.*p(:,1), ...
                            real(ifft(DFT.disc.T.*fft(p(:,2)))) + V.Vup.*p(:,2), ...
                            real(ifft(DFT.disc.T.*fft(p(:,3)))) + V.Vup.*p(:,3), ...
                            real(ifft(DFT.disc.T.*fft(p(:,4)))) + V.Vdw.*p(:,4), ...
                            real(ifft(DFT.disc.T.*fft(p(:,5)))) + V.Vdw.*p(:,5)];
    E_test              =  [sum(Hp(:,1).*p(:,1)); sum(Hp(:,2).*p(:,2)); sum(Hp(:,3).*p(:,3)); ...
                            sum(Hp(:,4).*p(:,4)); sum(Hp(:,5).*p(:,5))]*dx;
    err_test            =   sqrt([sum(abs(Hp(:,1)-E_test(1)*p(:,1)).^2); ...
                                  sum(abs(Hp(:,2)-E_test(2)*p(:,2)).^2); ...
                                  sum(abs(Hp(:,3)-E_test(3)*p(:,3)).^2); ...
                                  sum(abs(Hp(:,4)-E_test(4)*p(:,4)).^2); ...
                                  sum(abs(Hp(:,5)-E_test(5)*p(:,5)).^2)]*dx);
    obj.showResult('DFT_energyOrbital (real orbitals)',  ...
                            max(abs(E{1}   - E_test(1:3))  ) < 1e-10   && ...
                            max(abs(err{1} - err_test(1:3))) < 1e-10   && ...
                            max(abs(E{2}   - E_test(4:5))  ) < 1e-10   && ...
                            max(abs(err{2} - err_test(4:5))) < 1e-10);
    
    DFT.KSO.KSOup(:,1)  =   DFT.KSO.KSOup(:,1) * exp(1i*.3);
    DFT.KSO.KSOup(:,2)  =   DFT.KSO.KSOup(:,2) * exp(1i*.1);
    DFT.KSO.KSOup(:,3)  =   DFT.KSO.KSOup(:,3) * exp(1i*.5);
    DFT.KSO.KSOdw(:,1)  =   DFT.KSO.KSOdw(:,1) * exp(1i*.8);
    DFT.KSO.KSOdw(:,2)  =   DFT.KSO.KSOdw(:,2) * exp(1i*.7);
    [E,err]             =   DFT.disc.DFT_energyOrbital(V,DFT.KSO);

    obj.showResult('DFT_energyOrbital (complex orbitals)',  ...
                            max(abs(E{1}   - E_test(1:3))  ) < 1e-10   && ...
                            max(abs(err{1} - err_test(1:3))) < 1e-10   && ...
                            max(abs(E{2}   - E_test(4:5))  ) < 1e-10   && ...
                            max(abs(err{2} - err_test(4:5))) < 1e-10);

    % Velocity gauge
    DFT.disc.setTv(A);

    DFT.KSO.KSOup       =   p(:,1:3);
    DFT.KSO.KSOdw       =   p(:,4:5);
    [E,err]             =   DFT.disc.DFT_energyOrbital(V,DFT.KSO);
    
    Hpv                 =  [ifft(DFT.disc.Tv.*fft(p(:,1))) + V.Vup.*p(:,1), ...
                            ifft(DFT.disc.Tv.*fft(p(:,2))) + V.Vup.*p(:,2), ...
                            ifft(DFT.disc.Tv.*fft(p(:,3))) + V.Vup.*p(:,3), ...
                            ifft(DFT.disc.Tv.*fft(p(:,4))) + V.Vdw.*p(:,4), ...
                            ifft(DFT.disc.Tv.*fft(p(:,5))) + V.Vdw.*p(:,5)];
    E_test              =   E_test + .5*A^2;
    err_test            =   sqrt([sum(abs(Hpv(:,1)-E_test(1)*p(:,1)).^2); ...
                                  sum(abs(Hpv(:,2)-E_test(2)*p(:,2)).^2); ...
                                  sum(abs(Hpv(:,3)-E_test(3)*p(:,3)).^2); ...
                                  sum(abs(Hpv(:,4)-E_test(4)*p(:,4)).^2); ...
                                  sum(abs(Hpv(:,5)-E_test(5)*p(:,5)).^2)]*dx);
    obj.showResult('DFT_energyOrbital (velocity gauge, real orbitals)',  ...
                            max(abs(E{1}   - E_test(1:3))  ) < 1e-10   && ...
                            max(abs(err{1} - err_test(1:3))) < 1e-10   && ...
                            max(abs(E{2}   - E_test(4:5))  ) < 1e-10   && ...
                            max(abs(err{2} - err_test(4:5))) < 1e-10);
    
    DFT.KSO.KSOup(:,1)  =   DFT.KSO.KSOup(:,1) * exp(1i*.3);
    DFT.KSO.KSOup(:,2)  =   DFT.KSO.KSOup(:,2) * exp(1i*.1);
    DFT.KSO.KSOup(:,3)  =   DFT.KSO.KSOup(:,3) * exp(1i*.5);
    DFT.KSO.KSOdw(:,1)  =   DFT.KSO.KSOdw(:,1) * exp(1i*.8);
    DFT.KSO.KSOdw(:,2)  =   DFT.KSO.KSOdw(:,2) * exp(1i*.7);
    [E,err]             =   DFT.disc.DFT_energyOrbital(V,DFT.KSO);

    obj.showResult('DFT_energyOrbital (velocity gauge, complex orbitals)',  ...
                            max(abs(E{1}   - E_test(1:3))  ) < 1e-10   && ...
                            max(abs(err{1} - err_test(1:3))) < 1e-10   && ...
                            max(abs(E{2}   - E_test(4:5))  ) < 1e-10   && ...
                            max(abs(err{2} - err_test(4:5))) < 1e-10);

    DFT.disc.setTv([]);

    % DFT_operatorHamiltonian =============================================
    x0                  =   5;
    s                   =   2;
    p                   =   exp(-(DFT.x(:)-x0).^2 * .5/s^2);
    D2p                 =  -.5 * ((DFT.x(:)-x0).^2/s^4 - 1/s^2) .* p;

    A                   =   .1;
    D2pv                =   D2p + .5*A^2*p + A*1i/s^2 * (DFT.x(:)-x0) .*p;

    V                   =   Vext.getPotential;
    V.add(0,1);
    Vupp                =   V.Vup .* p;
    Vdwp                =   V.Vdw .* p;

    % Length gauge
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,true);
    R                   =   max(abs(Hp-D2p-Vupp)) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (spin up, no symmetry)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,false);
    R                   =   max(abs(Hp-D2p-Vdwp)) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (spin down, no symmetry)' ,R);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,true,1);
    R                   =   max(abs(Hp-.5*(D2p+Vupp+flip(D2p+Vupp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (spin up, symmetrized)' ,R);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,false,1);
    R                   =   max(abs(Hp-.5*(D2p+Vdwp+flip(D2p+Vdwp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (spin down, symmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,true,-1);
    R                   =   max(abs(Hp-.5*(D2p+Vupp-flip(D2p+Vupp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (spin up, antisymmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,false,-1);
    R                   =   max(abs(Hp-.5*(D2p+Vdwp-flip(D2p+Vdwp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (spin down, antisymmetrized)' ,R);

    % Velocity gauge
    DFT.disc.setTv(A);
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,true);
    R                   =   max(abs(Hp-D2pv-Vupp)) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (velocity, spin up, no symmetry)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,false);
    R                   =   max(abs(Hp-D2pv-Vdwp)) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (velocity, spin down, no symmetry)' ,R);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,true,1);
    R                   =   max(abs(Hp-.5*(D2pv+Vupp+flip(D2pv+Vupp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (velocity, spin up, symmetrized)' ,R);

    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,false,1);
    R                   =   max(abs(Hp-.5*(D2pv+Vdwp+flip(D2pv+Vdwp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (velocity, spin down, symmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,true,-1);
    R                   =   max(abs(Hp-.5*(D2pv+Vupp-flip(D2pv+Vupp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (velocity, spin up, antisymmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorHamiltonian(V,p,false,-1);
    R                   =   max(abs(Hp-.5*(D2pv+Vdwp-flip(D2pv+Vdwp)))) < 1e-10;
    obj.showResult('DFT_operatorHamiltonian (velocity, spin down, antisymmetrized)' ,R);

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
    V.add(0,1);
    Vupp                =   V.Vup .* p;
    Vdwp                =   V.Vdw .* p;

    % (same for length and velocity gauges)
    Hp                  =   DFT.disc.DFT_operatorPotential(V,p,true);
    R                   =   max(abs(Hp-Vupp)) < 1e-10;
    obj.showResult('DFT_operatorPotential (spin up, no symmetry)',R);

    Hp                  =   DFT.disc.DFT_operatorPotential(V,p,false);
    R                   =   max(abs(Hp-Vdwp)) < 1e-10;
    obj.showResult('DFT_operatorPotential (spin down, no symmetry)',R);
    
    Hp                  =   DFT.disc.DFT_operatorPotential(V,p,true,1);
    R                   =   max(abs(Hp-.5*(Vupp+flip(Vupp)))) < 1e-10;
    obj.showResult('DFT_operatorPotential (spin up, symmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorPotential(V,p,false,1);
    R                   =   max(abs(Hp-.5*(Vdwp+flip(Vdwp)))) < 1e-10;
    obj.showResult('DFT_operatorPotential (spin down, symmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorPotential(V,p,true,-1);
    R                   =   max(abs(Hp-.5*(Vupp-flip(Vupp)))) < 1e-10;
    obj.showResult('DFT_operatorPotential (spin up, antisymmetrized)' ,R);
    
    Hp                  =   DFT.disc.DFT_operatorPotential(V,p,false,-1);
    R                   =   max(abs(Hp-.5*(Vdwp-flip(Vdwp)))) < 1e-10;
    obj.showResult('DFT_operatorPotential (spin down, antisymmetrized)' ,R);

end