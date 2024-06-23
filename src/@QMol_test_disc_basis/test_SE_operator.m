function test_SE_operator(obj)
%test_SE_operator unit tests for the operator overload in the class

    % Initialization
    obj.showSection('Schrodinger-equation operators and functionals');

    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];
    d_R                 =   QMol_disc_basis('xspan',x,'basis',v);   d_R.orthonormalizeBasis;

    SE_R                =   QMol_SE(...
                                'discretization',       d_R,                            ...
                                'numberWaveFunction',   3,                              ...
                                'potential',            QMol_SE_V('atom',               ...
                                   {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5),  ...
                                    QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3),  ...
                                    QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)})  );
    SE_R.initialize;

    v                   =   v.*[exp( 2i*x),exp( 3i*x),exp(-4i*x),exp(-2i*x)];
    d_C                 =   QMol_disc_basis('xspan',x,'basis',v);   d_C.orthonormalizeBasis;
    SE_C                =   QMol_SE(...
                                'disc',                 d_C,                            ...
                                'numberWaveFunction',   3,                              ...
                                'potential',            QMol_SE_V('atom',               ...
                                   {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5),  ...
                                    QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3),  ...
                                    QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)})  );
    SE_C.initialize;

    A                   =   .1;

    % DFT_energyKinetic ===================================================
    wfcn_R              =   d_R.SE_allocateWaveFunction(3,[],randStr,true);
    wfcn_C              =   d_R.SE_allocateWaveFunction(5,[],randStr,false);

    % Length gauge
    rWfcn               =   d_R.SE_reconstructWaveFunction(wfcn_R);
    E                   =   d_R.SE_energyKinetic(wfcn_R);
    E_test              =   sum( rWfcn.wfcn(:,1) .* real(ifft(d_R.T.*fft(rWfcn.wfcn(:,1)))) ) + ...
                            sum( rWfcn.wfcn(:,2) .* real(ifft(d_R.T.*fft(rWfcn.wfcn(:,2)))) ) + ...
                            sum( rWfcn.wfcn(:,3) .* real(ifft(d_R.T.*fft(rWfcn.wfcn(:,3)))) );
    E_test              =   E_test * d_R.dx;
    obj.showResult('SE_energyKinetic (real wave functions, real basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_C.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    E                   =   d_C.SE_energyKinetic(wfcn_R);
    E_test              =   real(sum( conj(rWfcn.wfcn(:,1)) .* ifft(d_C.T.*fft(rWfcn.wfcn(:,1)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,2)) .* ifft(d_C.T.*fft(rWfcn.wfcn(:,2)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,3)) .* ifft(d_C.T.*fft(rWfcn.wfcn(:,3)))) );
    E_test              =   E_test * d_C.dx;
    obj.showResult('SE_energyKinetic (real wave functions, complex basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_R.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    E                   =   d_R.SE_energyKinetic(wfcn_C);
    E_test              =   real(sum( conj(rWfcn.wfcn(:,1)) .* ifft(d_R.T.*fft(rWfcn.wfcn(:,1)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,2)) .* ifft(d_R.T.*fft(rWfcn.wfcn(:,2)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,3)) .* ifft(d_R.T.*fft(rWfcn.wfcn(:,3)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,4)) .* ifft(d_R.T.*fft(rWfcn.wfcn(:,4)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,5)) .* ifft(d_R.T.*fft(rWfcn.wfcn(:,5)))) );
    E_test              =   E_test * d_R.dx;
    obj.showResult('SE_energyKinetic (complex wave functions, real basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_C.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    E                   =   d_C.SE_energyKinetic(wfcn_C);
    E_test              =   real(sum( conj(rWfcn.wfcn(:,1)) .* ifft(d_C.T.*fft(rWfcn.wfcn(:,1)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,2)) .* ifft(d_C.T.*fft(rWfcn.wfcn(:,2)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,3)) .* ifft(d_C.T.*fft(rWfcn.wfcn(:,3)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,4)) .* ifft(d_R.T.*fft(rWfcn.wfcn(:,4)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,5)) .* ifft(d_R.T.*fft(rWfcn.wfcn(:,5)))) );
    E_test              =   E_test * d_C.dx;
    obj.showResult('SE_energyKinetic (complex wave functions, complex basis)',isreal(E) && abs(E-E_test) < 1e-10);

    % Velocity gauge
    d_R.setTv(A);           d_C.setTv(A);

    d_R.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    E                   =   d_R.SE_energyKinetic(wfcn_R);
    E_test              =   sum( rWfcn.wfcn(:,1) .* real(ifft(d_R.Tv.*fft(rWfcn.wfcn(:,1)))) ) + ...
                            sum( rWfcn.wfcn(:,2) .* real(ifft(d_R.Tv.*fft(rWfcn.wfcn(:,2)))) ) + ...
                            sum( rWfcn.wfcn(:,3) .* real(ifft(d_R.Tv.*fft(rWfcn.wfcn(:,3)))) );
    E_test              =   E_test * d_R.dx;
    obj.showResult('SE_energyKinetic (velocity gauge, real wfcn, real basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_C.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    E                   =   d_C.SE_energyKinetic(wfcn_R);
    E_test              =   real(sum( conj(rWfcn.wfcn(:,1)) .* ifft(d_C.Tv.*fft(rWfcn.wfcn(:,1)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,2)) .* ifft(d_C.Tv.*fft(rWfcn.wfcn(:,2)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,3)) .* ifft(d_C.Tv.*fft(rWfcn.wfcn(:,3)))) );
    E_test              =   E_test * d_C.dx;
    obj.showResult('SE_energyKinetic (velocity gauge, real wfcn, complex basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_R.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    E                   =   d_R.SE_energyKinetic(wfcn_C);
    E_test              =   real(sum( conj(rWfcn.wfcn(:,1)) .* ifft(d_R.Tv.*fft(rWfcn.wfcn(:,1)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,2)) .* ifft(d_R.Tv.*fft(rWfcn.wfcn(:,2)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,3)) .* ifft(d_R.Tv.*fft(rWfcn.wfcn(:,3)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,4)) .* ifft(d_R.Tv.*fft(rWfcn.wfcn(:,4)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,5)) .* ifft(d_R.Tv.*fft(rWfcn.wfcn(:,5)))) );
    E_test              =   E_test * d_R.dx;
    obj.showResult('SE_energyKinetic (velocity gauge, complex wfcn, real basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_C.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    E                   =   d_C.SE_energyKinetic(wfcn_C);
    E_test              =   real(sum( conj(rWfcn.wfcn(:,1)) .* ifft(d_C.Tv.*fft(rWfcn.wfcn(:,1)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,2)) .* ifft(d_C.Tv.*fft(rWfcn.wfcn(:,2)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,3)) .* ifft(d_C.Tv.*fft(rWfcn.wfcn(:,3)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,4)) .* ifft(d_R.Tv.*fft(rWfcn.wfcn(:,4)))) ) + ...
                            real(sum( conj(rWfcn.wfcn(:,5)) .* ifft(d_R.Tv.*fft(rWfcn.wfcn(:,5)))) );
    E_test              =   E_test * d_C.dx;
    obj.showResult('SE_energyKinetic (velocity gauge, complex wfcn, complex basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_R.setTv([]);          d_C.setTv([]);

    % NOTE: if this test fails, check the following methods first
    %       - DFT_allocateOrbital
    %       - DFT_reconstructOrbital

    % SE_energyWaveFunction ===============================================
    d_R.SE_allocateWaveFunction(3,wfcn_R,randStr,true);
    d_R.SE_allocateWaveFunction(5,wfcn_C,randStr,false);

    % Length gauge
    d_R.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    [E,err]             =   d_R.SE_energyWaveFunction(SE_R.V,wfcn_R);
    
    Hp                  =  [real(ifft(d_R.T.*fft(rWfcn.wfcn(:,1)))) + SE_R.V.V.*rWfcn.wfcn(:,1), ...
                            real(ifft(d_R.T.*fft(rWfcn.wfcn(:,2)))) + SE_R.V.V.*rWfcn.wfcn(:,2), ...
                            real(ifft(d_R.T.*fft(rWfcn.wfcn(:,3)))) + SE_R.V.V.*rWfcn.wfcn(:,3)];
    E_test              =  [sum(Hp(:,1).*rWfcn.wfcn(:,1)); ...
                            sum(Hp(:,2).*rWfcn.wfcn(:,2)); ...
                            sum(Hp(:,3).*rWfcn.wfcn(:,3))]*d_R.dx;

    rWfcn.set('wfcn',Hp);   pWfcn   =   d_R.SE_projectWaveFunction(rWfcn);           % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pWfcn.wfcn(:,1)-E_test(1)*wfcn_R.wfcn(:,1)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,2)-E_test(2)*wfcn_R.wfcn(:,2)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,3)-E_test(3)*wfcn_R.wfcn(:,3)).^2)]);
    obj.showResult('SE_energyWaveFunction (real wave functions, real basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    
    d_C.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    [E,err]             =   d_C.SE_energyWaveFunction(SE_C.V,wfcn_R);
    
    Hp                  =  [ifft(d_C.T.*fft(rWfcn.wfcn(:,1))) + SE_C.V.V.*rWfcn.wfcn(:,1), ...
                            ifft(d_C.T.*fft(rWfcn.wfcn(:,2))) + SE_C.V.V.*rWfcn.wfcn(:,2), ...
                            ifft(d_C.T.*fft(rWfcn.wfcn(:,3))) + SE_C.V.V.*rWfcn.wfcn(:,3)];
    E_test              =  [sum(Hp(:,1).*conj(rWfcn.wfcn(:,1))); ...
                            sum(Hp(:,2).*conj(rWfcn.wfcn(:,2))); ...
                            sum(Hp(:,3).*conj(rWfcn.wfcn(:,3)))]*d_C.dx;

    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pWfcn.wfcn(:,1)-E_test(1)*wfcn_R.wfcn(:,1)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,2)-E_test(2)*wfcn_R.wfcn(:,2)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,3)-E_test(3)*wfcn_R.wfcn(:,3)).^2)]);
    obj.showResult('SE_energyWaveFunction (real wave functions, complex basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));


    d_R.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    [E,err]             =   d_R.SE_energyWaveFunction(SE_R.V,wfcn_C);
    
    Hp                  =  [ifft(d_R.T.*fft(rWfcn.wfcn(:,1))) + SE_R.V.V.*rWfcn.wfcn(:,1), ...
                            ifft(d_R.T.*fft(rWfcn.wfcn(:,2))) + SE_R.V.V.*rWfcn.wfcn(:,2), ...
                            ifft(d_R.T.*fft(rWfcn.wfcn(:,3))) + SE_R.V.V.*rWfcn.wfcn(:,3), ...
                            ifft(d_R.T.*fft(rWfcn.wfcn(:,4))) + SE_R.V.V.*rWfcn.wfcn(:,4), ...
                            ifft(d_R.T.*fft(rWfcn.wfcn(:,5))) + SE_R.V.V.*rWfcn.wfcn(:,5)];
    E_test              =  [sum(Hp(:,1).*conj(rWfcn.wfcn(:,1))); ...
                            sum(Hp(:,2).*conj(rWfcn.wfcn(:,2))); ...
                            sum(Hp(:,3).*conj(rWfcn.wfcn(:,3))); ...
                            sum(Hp(:,4).*conj(rWfcn.wfcn(:,4))); ...
                            sum(Hp(:,5).*conj(rWfcn.wfcn(:,5)))]*d_R.dx;

    rWfcn.set('wfcn',Hp);     d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pWfcn.wfcn(:,1)-E_test(1)*wfcn_C.wfcn(:,1)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,2)-E_test(2)*wfcn_C.wfcn(:,2)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,3)-E_test(3)*wfcn_C.wfcn(:,3)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,4)-E_test(4)*wfcn_C.wfcn(:,4)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,5)-E_test(5)*wfcn_C.wfcn(:,5)).^2)]);
    obj.showResult('SE_energyWaveFunction (complex wave functions, real basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    
    d_C.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    [E,err]             =   d_C.SE_energyWaveFunction(SE_C.V,wfcn_C);
    
    Hp                  =  [ifft(d_C.T.*fft(rWfcn.wfcn(:,1))) + SE_C.V.V.*rWfcn.wfcn(:,1), ...
                            ifft(d_C.T.*fft(rWfcn.wfcn(:,2))) + SE_C.V.V.*rWfcn.wfcn(:,2), ...
                            ifft(d_C.T.*fft(rWfcn.wfcn(:,3))) + SE_C.V.V.*rWfcn.wfcn(:,3), ...
                            ifft(d_R.T.*fft(rWfcn.wfcn(:,4))) + SE_R.V.V.*rWfcn.wfcn(:,4), ...
                            ifft(d_R.T.*fft(rWfcn.wfcn(:,5))) + SE_R.V.V.*rWfcn.wfcn(:,5)];
    E_test              =  [sum(Hp(:,1).*conj(rWfcn.wfcn(:,1))); ...
                            sum(Hp(:,2).*conj(rWfcn.wfcn(:,2))); ...
                            sum(Hp(:,3).*conj(rWfcn.wfcn(:,3))); ...
                            sum(Hp(:,4).*conj(rWfcn.wfcn(:,4))); ...
                            sum(Hp(:,5).*conj(rWfcn.wfcn(:,5)))]*d_C.dx;

    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pWfcn.wfcn(:,1)-E_test(1)*wfcn_C.wfcn(:,1)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,2)-E_test(2)*wfcn_C.wfcn(:,2)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,3)-E_test(3)*wfcn_C.wfcn(:,3)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,4)-E_test(4)*wfcn_C.wfcn(:,4)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,5)-E_test(5)*wfcn_C.wfcn(:,5)).^2)]);
    obj.showResult('SE_energyWaveFunction (complex wave functions, complex basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    % Velocity gauge
    d_R.setTv(A);           d_C.setTv(A);

    d_R.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    [E,err]             =   d_R.SE_energyWaveFunction(SE_R.V,wfcn_R);
    
    Hp                  =  [ifft(d_R.Tv.*fft(rWfcn.wfcn(:,1))) + SE_R.V.V.*rWfcn.wfcn(:,1), ...
                            ifft(d_R.Tv.*fft(rWfcn.wfcn(:,2))) + SE_R.V.V.*rWfcn.wfcn(:,2), ...
                            ifft(d_R.Tv.*fft(rWfcn.wfcn(:,3))) + SE_R.V.V.*rWfcn.wfcn(:,3)];
    E_test              =  [sum(Hp(:,1).*conj(rWfcn.wfcn(:,1))); ...
                            sum(Hp(:,2).*conj(rWfcn.wfcn(:,2))); ...
                            sum(Hp(:,3).*conj(rWfcn.wfcn(:,3)))]*d_R.dx;

    rWfcn.set('wfcn',Hp);   d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pWfcn.wfcn(:,1)-E_test(1)*wfcn_R.wfcn(:,1)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,2)-E_test(2)*wfcn_R.wfcn(:,2)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,3)-E_test(3)*wfcn_R.wfcn(:,3)).^2)]);
    obj.showResult('SE_energyWaveFunction (vel. g., real wfcn, real basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    
    d_C.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    [E,err]             =   d_C.SE_energyWaveFunction(SE_C.V,wfcn_R);
    
    Hp                  =  [ifft(d_C.Tv.*fft(rWfcn.wfcn(:,1))) + SE_C.V.V.*rWfcn.wfcn(:,1), ...
                            ifft(d_C.Tv.*fft(rWfcn.wfcn(:,2))) + SE_C.V.V.*rWfcn.wfcn(:,2), ...
                            ifft(d_C.Tv.*fft(rWfcn.wfcn(:,3))) + SE_C.V.V.*rWfcn.wfcn(:,3)];
    E_test              =  [sum(Hp(:,1).*conj(rWfcn.wfcn(:,1))); ...
                            sum(Hp(:,2).*conj(rWfcn.wfcn(:,2))); ...
                            sum(Hp(:,3).*conj(rWfcn.wfcn(:,3)))]*d_C.dx;

    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pWfcn.wfcn(:,1)-E_test(1)*wfcn_R.wfcn(:,1)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,2)-E_test(2)*wfcn_R.wfcn(:,2)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,3)-E_test(3)*wfcn_R.wfcn(:,3)).^2)]);
    obj.showResult('SE_energyWaveFunction (vel. g., real wfcn, complex basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));


    d_R.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    [E,err]             =   d_R.SE_energyWaveFunction(SE_R.V,wfcn_C);
    
    Hp                  =  [ifft(d_R.Tv.*fft(rWfcn.wfcn(:,1))) + SE_R.V.V.*rWfcn.wfcn(:,1), ...
                            ifft(d_R.Tv.*fft(rWfcn.wfcn(:,2))) + SE_R.V.V.*rWfcn.wfcn(:,2), ...
                            ifft(d_R.Tv.*fft(rWfcn.wfcn(:,3))) + SE_R.V.V.*rWfcn.wfcn(:,3), ...
                            ifft(d_R.Tv.*fft(rWfcn.wfcn(:,4))) + SE_R.V.V.*rWfcn.wfcn(:,4), ...
                            ifft(d_R.Tv.*fft(rWfcn.wfcn(:,5))) + SE_R.V.V.*rWfcn.wfcn(:,5)];
    E_test              =  [sum(Hp(:,1).*conj(rWfcn.wfcn(:,1))); ...
                            sum(Hp(:,2).*conj(rWfcn.wfcn(:,2))); ...
                            sum(Hp(:,3).*conj(rWfcn.wfcn(:,3))); ...
                            sum(Hp(:,4).*conj(rWfcn.wfcn(:,4))); ...
                            sum(Hp(:,5).*conj(rWfcn.wfcn(:,5)))]*d_R.dx;

    rWfcn.set('wfcn',Hp);     d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pWfcn.wfcn(:,1)-E_test(1)*wfcn_C.wfcn(:,1)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,2)-E_test(2)*wfcn_C.wfcn(:,2)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,3)-E_test(3)*wfcn_C.wfcn(:,3)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,4)-E_test(4)*wfcn_C.wfcn(:,4)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,5)-E_test(5)*wfcn_C.wfcn(:,5)).^2)]);
    obj.showResult('SE_energyWaveFunction (vel. g., complex wfcn, real basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    
    d_C.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    [E,err]             =   d_C.SE_energyWaveFunction(SE_C.V,wfcn_C);
    
    Hp                  =  [ifft(d_C.Tv.*fft(rWfcn.wfcn(:,1))) + SE_C.V.V.*rWfcn.wfcn(:,1), ...
                            ifft(d_C.Tv.*fft(rWfcn.wfcn(:,2))) + SE_C.V.V.*rWfcn.wfcn(:,2), ...
                            ifft(d_C.Tv.*fft(rWfcn.wfcn(:,3))) + SE_C.V.V.*rWfcn.wfcn(:,3), ...
                            ifft(d_C.Tv.*fft(rWfcn.wfcn(:,4))) + SE_C.V.V.*rWfcn.wfcn(:,4), ...
                            ifft(d_C.Tv.*fft(rWfcn.wfcn(:,5))) + SE_C.V.V.*rWfcn.wfcn(:,5)];
    E_test              =  [sum(Hp(:,1).*conj(rWfcn.wfcn(:,1))); ...
                            sum(Hp(:,2).*conj(rWfcn.wfcn(:,2))); ...
                            sum(Hp(:,3).*conj(rWfcn.wfcn(:,3))); ...
                            sum(Hp(:,4).*conj(rWfcn.wfcn(:,4))); ...
                            sum(Hp(:,5).*conj(rWfcn.wfcn(:,5)))]*d_C.dx;

    rWfcn.set('wfcn',Hp);     d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pWfcn.wfcn(:,1)-E_test(1)*wfcn_C.wfcn(:,1)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,2)-E_test(2)*wfcn_C.wfcn(:,2)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,3)-E_test(3)*wfcn_C.wfcn(:,3)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,4)-E_test(4)*wfcn_C.wfcn(:,4)).^2); ...
                                  sum(abs(pWfcn.wfcn(:,5)-E_test(5)*wfcn_C.wfcn(:,5)).^2)]);
    obj.showResult('SE_energyWaveFunction (vel. g. , complex wfcn, complex basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));


    d_R.setTv([]);          d_C.setTv([]);

    % NOTE: if this test fails, check the following methods first
    %       - DFT_allocateOrbital
    %       - DFT_reconstructOrbital
    %       - DFT_projectOrbital

    % DFT_operatorHamiltonian =============================================
    d_R.SE_allocateWaveFunction(1,wfcn_R,randStr,true);
    d_R.SE_allocateWaveFunction(1,wfcn_C,randStr,false);

    % Length gauge
    d_R.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    Hp                  =   ifft(d_R.T.*fft(rWfcn.wfcn)) + SE_R.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_R.SE_operatorHamiltonian(SE_R.V,wfcn_R.wfcn);
    obj.showResult('SE_operatorHamiltonian (real wave functions, real basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_C.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    Hp                  =   ifft(d_C.T.*fft(rWfcn.wfcn)) + SE_C.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_C.SE_operatorHamiltonian(SE_C.V,wfcn_R.wfcn);
    obj.showResult('SE_operatorHamiltonian (real wave functions, complex basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_R.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    Hp                  =   ifft(d_R.T.*fft(rWfcn.wfcn)) + SE_R.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_R.SE_operatorHamiltonian(SE_R.V,wfcn_C.wfcn);
    obj.showResult('SE_operatorHamiltonian (complex wave functions, real basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_C.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    Hp                  =   ifft(d_C.T.*fft(rWfcn.wfcn)) + SE_C.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_C.SE_operatorHamiltonian(SE_C.V,wfcn_C.wfcn);
    obj.showResult('SE_operatorHamiltonian (complex wave functions, complex basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    % Velocity gauge
    d_R.setTv(A);           d_C.setTv(A);

    d_R.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    Hp                  =   ifft(d_R.Tv.*fft(rWfcn.wfcn)) + SE_R.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_R.SE_operatorHamiltonian(SE_R.V,wfcn_R.wfcn);
    obj.showResult('SE_operatorHamiltonian (vel. g., real wfcn, real basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_C.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    Hp                  =   ifft(d_C.Tv.*fft(rWfcn.wfcn)) + SE_C.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_C.SE_operatorHamiltonian(SE_C.V,wfcn_R.wfcn);
    obj.showResult('SE_operatorHamiltonian (vel. g., real wfcn, complex basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_R.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    Hp                  =   ifft(d_R.Tv.*fft(rWfcn.wfcn)) + SE_R.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_R.SE_operatorHamiltonian(SE_R.V,wfcn_C.wfcn);
    obj.showResult('SE_operatorHamiltonian (vel. g., complex wfcn, real basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_C.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    Hp                  =   ifft(d_C.Tv.*fft(rWfcn.wfcn)) + SE_C.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_C.SE_operatorHamiltonian(SE_C.V,wfcn_C.wfcn);
    obj.showResult('SE_operatorHamiltonian (vel. g., complex wfcn, complex basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_R.setTv([]);          d_C.setTv([]);

    % NOTE: if this test fails, check the following methods first
    %       - DFT_allocateOrbital
    %       - DFT_reconstructOrbital
    %       - DFT_projectOrbital

    % DFT_operatorKinetic =================================================
    d_R.SE_allocateWaveFunction(1,wfcn_R,randStr,true);
    d_R.SE_allocateWaveFunction(1,wfcn_C,randStr,false);

    % Length gauge
    d_R.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    Hp                  =   ifft(d_R.T.*fft(rWfcn.wfcn));
    rWfcn.set('wfcn',Hp);   d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_R.SE_operatorKinetic(wfcn_R.wfcn);
    obj.showResult('SE_operatorKinetic (real wave functions, real basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_C.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    Hp                  =   ifft(d_C.T.*fft(rWfcn.wfcn));
    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_C.SE_operatorKinetic(wfcn_R.wfcn);
    obj.showResult('SE_operatorKinetic (real wave functions, complex basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_R.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    Hp                  =   ifft(d_R.T.*fft(rWfcn.wfcn));
    rWfcn.set('wfcn',Hp);   d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_R.SE_operatorKinetic(wfcn_C.wfcn);
    obj.showResult('SE_operatorKinetic (complex wave functions, real basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_C.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    Hp                  =   ifft(d_C.T.*fft(rWfcn.wfcn));
    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_C.SE_operatorKinetic(wfcn_C.wfcn);
    obj.showResult('SE_operatorKinetic (complex wave functions, complex basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    % Velocity gauge
    d_R.setTv(A);           d_C.setTv(A);

    d_R.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    Hp                  =   ifft(d_R.Tv.*fft(rWfcn.wfcn));
    rWfcn.set('wfcn',Hp);   d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_R.SE_operatorKinetic(wfcn_R.wfcn);
    obj.showResult('SE_operatorKinetic (vel. g., real wfcn, real basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_C.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    Hp                  =   ifft(d_C.Tv.*fft(rWfcn.wfcn));
    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_C.SE_operatorKinetic(wfcn_R.wfcn);
    obj.showResult('SE_operatorKinetic (vel. g., real wfcn, complex basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_R.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    Hp                  =   ifft(d_R.Tv.*fft(rWfcn.wfcn));
    rWfcn.set('wfcn',Hp);   d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_R.SE_operatorKinetic(wfcn_C.wfcn);
    obj.showResult('SE_operatorKinetic (vel. g., complex wfcn, real basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_C.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    Hp                  =   ifft(d_C.Tv.*fft(rWfcn.wfcn));
    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_C.SE_operatorKinetic(wfcn_C.wfcn);
    obj.showResult('SE_operatorKinetic (vel. g., complex wfcn, complex basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_R.setTv([]);          d_C.setTv([]);

    % NOTE: if this test fails, check the following methods first
    %       - DFT_allocateOrbital
    %       - DFT_reconstructOrbital
    %       - DFT_projectOrbital


    % DFT_operatorPotential ===============================================
    d_R.SE_allocateWaveFunction(1,wfcn_R,randStr,true);
    d_R.SE_allocateWaveFunction(1,wfcn_C,randStr,false);

    d_R.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    Hp                  =   SE_R.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_R.SE_operatorPotential(SE_R.V,wfcn_R.wfcn);
    obj.showResult('SE_operatorPotential (real wave functions, real basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_C.SE_reconstructWaveFunction(wfcn_R,rWfcn);
    Hp                  =   SE_C.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_C.SE_operatorPotential(SE_C.V,wfcn_R.wfcn);
    obj.showResult('SE_operatorPotential (real wave functions, complex basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_R.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    Hp                  =   SE_R.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_R.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_R.SE_operatorPotential(SE_R.V,wfcn_C.wfcn);
    obj.showResult('SE_operatorPotential (complex wave functions, real basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    d_C.SE_reconstructWaveFunction(wfcn_C,rWfcn);
    Hp                  =   SE_C.V.V.*rWfcn.wfcn;
    rWfcn.set('wfcn',Hp);   d_C.SE_projectWaveFunction(rWfcn,pWfcn);              % Project Hp onto basis
    Hp                  =   d_C.SE_operatorPotential(SE_C.V,wfcn_C.wfcn);
    obj.showResult('SE_operatorPotential (complex wave functions, complex basis)',...
                            all(abs(Hp - pWfcn.wfcn) < 1e-10,'all'));

    % NOTE: if this test fails, check the following methods first
    %       - DFT_allocateOrbital
    %       - DFT_reconstructOrbital
    %       - DFT_projectOrbital

end