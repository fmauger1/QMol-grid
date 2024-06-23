function test_DFT_operator_spinRes(obj)
%test_operator unit tests for the operator overload in the class

    % Initialization
    obj.showSection('DFT operators and functionals (spin restricted)');

    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];
    d_R                 =   QMol_disc_basis('xspan',x,'basis',v);   d_R.orthonormalizeBasis;

    DFT_R               =   QMol_DFT_spinRes(...
                                'discretization',       d_R,                            ...
                                'occupation',           [2 1 2],                        ...
                                'externalPotential',    QMol_DFT_Vext('atom',           ...
                                   {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5),  ...
                                    QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3),  ...
                                    QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)}), ...
                                'HartreePotential',     QMol_DFT_Vh_conv,               ...
                                'exchangeCorrelationPotential',     []);
    DFT_R.initialize;

    v                   =   v.*[exp( 2i*x),exp( 3i*x),exp(-4i*x),exp(-2i*x)];
    d_C                 =   QMol_disc_basis('xspan',x,'basis',v);   d_C.orthonormalizeBasis;
    DFT_C               =   QMol_DFT_spinRes(...
                                'disc',                 d_C,                            ...
                                'occupation',           [2 1 1],                        ...
                                'externalPotential',    QMol_DFT_Vext('atom',           ...
                                   {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5),  ...
                                    QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3),  ...
                                    QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)}), ...
                                'HartreePotential',     QMol_DFT_Vh_conv,               ...
                                'exchangeCorrelationPotential',     []);
    DFT_C.initialize;

    A                   =   .1;

    % DFT_energyKinetic ===================================================
    KSO_R               =   d_R.DFT_allocateOrbital(3,[],randStr,true);
    KSO_C               =   d_R.DFT_allocateOrbital(3,[],randStr,false);

    % Length gauge
    rKSO                =   d_R.DFT_reconstructOrbital(KSO_R);
    E                   =   d_R.DFT_energyKinetic(DFT_R.occ,KSO_R);
    E_test              =  [sum( rKSO.KSO(:,1) .* real(ifft(d_R.T.*fft(rKSO.KSO(:,1)))) ), ...
                            sum( rKSO.KSO(:,2) .* real(ifft(d_R.T.*fft(rKSO.KSO(:,2)))) ), ...
                            sum( rKSO.KSO(:,3) .* real(ifft(d_R.T.*fft(rKSO.KSO(:,3)))) )];
    E_test              =   sum(DFT_R.occ.*E_test) * d_R.dx;
    obj.showResult('DFT_energyKinetic (real orbitals, real basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_C.DFT_reconstructOrbital(KSO_R,rKSO);
    E                   =   d_C.DFT_energyKinetic(DFT_C.occ,KSO_R);
    E_test              =  [real(sum( conj(rKSO.KSO(:,1)) .* ifft(d_C.T.*fft(rKSO.KSO(:,1)))) ), ...
                            real(sum( conj(rKSO.KSO(:,2)) .* ifft(d_C.T.*fft(rKSO.KSO(:,2)))) ), ...
                            real(sum( conj(rKSO.KSO(:,3)) .* ifft(d_C.T.*fft(rKSO.KSO(:,3)))) )];
    E_test              =   sum(DFT_C.occ.*E_test) * d_C.dx;
    obj.showResult('DFT_energyKinetic (real orbitals, complex basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_R.DFT_reconstructOrbital(KSO_C,rKSO);
    E                   =   d_R.DFT_energyKinetic(DFT_R.occ,KSO_C);
    E_test              =  [real(sum( conj(rKSO.KSO(:,1)) .* ifft(d_R.T.*fft(rKSO.KSO(:,1)))) ), ...
                            real(sum( conj(rKSO.KSO(:,2)) .* ifft(d_R.T.*fft(rKSO.KSO(:,2)))) ), ...
                            real(sum( conj(rKSO.KSO(:,3)) .* ifft(d_R.T.*fft(rKSO.KSO(:,3)))) )];
    E_test              =   sum(DFT_R.occ.*E_test) * d_R.dx;
    obj.showResult('DFT_energyKinetic (complex orbitals, real basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_C.DFT_reconstructOrbital(KSO_C,rKSO);
    E                   =   d_C.DFT_energyKinetic(DFT_C.occ,KSO_C);
    E_test              =  [real(sum( conj(rKSO.KSO(:,1)) .* ifft(d_C.T.*fft(rKSO.KSO(:,1)))) ), ...
                            real(sum( conj(rKSO.KSO(:,2)) .* ifft(d_C.T.*fft(rKSO.KSO(:,2)))) ), ...
                            real(sum( conj(rKSO.KSO(:,3)) .* ifft(d_C.T.*fft(rKSO.KSO(:,3)))) )];
    E_test              =   sum(DFT_C.occ.*E_test) * d_C.dx;
    obj.showResult('DFT_energyKinetic (complex orbitals, complex basis)',isreal(E) && abs(E-E_test) < 1e-10);

    % Velocity gauge
    d_R.setTv(A);           d_C.setTv(A);

    d_R.DFT_reconstructOrbital(KSO_R,rKSO);
    E                   =   d_R.DFT_energyKinetic(DFT_R.occ,KSO_R);
    E_test              =  [sum( rKSO.KSO(:,1) .* real(ifft(d_R.Tv.*fft(rKSO.KSO(:,1)))) ), ...
                            sum( rKSO.KSO(:,2) .* real(ifft(d_R.Tv.*fft(rKSO.KSO(:,2)))) ), ...
                            sum( rKSO.KSO(:,3) .* real(ifft(d_R.Tv.*fft(rKSO.KSO(:,3)))) )];
    E_test              =   sum(DFT_R.occ.*E_test) * d_R.dx;
    obj.showResult('DFT_energyKinetic (velocity gauge, real KSO, real basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_C.DFT_reconstructOrbital(KSO_R,rKSO);
    E                   =   d_C.DFT_energyKinetic(DFT_C.occ,KSO_R);
    E_test              =  [real(sum( conj(rKSO.KSO(:,1)) .* ifft(d_C.Tv.*fft(rKSO.KSO(:,1)))) ), ...
                            real(sum( conj(rKSO.KSO(:,2)) .* ifft(d_C.Tv.*fft(rKSO.KSO(:,2)))) ), ...
                            real(sum( conj(rKSO.KSO(:,3)) .* ifft(d_C.Tv.*fft(rKSO.KSO(:,3)))) )];
    E_test              =   sum(DFT_C.occ.*E_test) * d_C.dx;
    obj.showResult('DFT_energyKinetic (velocity gauge, real KSO, complex basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_R.DFT_reconstructOrbital(KSO_C,rKSO);
    E                   =   d_R.DFT_energyKinetic(DFT_R.occ,KSO_C);
    E_test              =  [real(sum( conj(rKSO.KSO(:,1)) .* ifft(d_R.Tv.*fft(rKSO.KSO(:,1)))) ), ...
                            real(sum( conj(rKSO.KSO(:,2)) .* ifft(d_R.Tv.*fft(rKSO.KSO(:,2)))) ), ...
                            real(sum( conj(rKSO.KSO(:,3)) .* ifft(d_R.Tv.*fft(rKSO.KSO(:,3)))) )];
    E_test              =   sum(DFT_R.occ.*E_test) * d_R.dx;
    obj.showResult('DFT_energyKinetic (velocity gauge, complex KSO, real basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_C.DFT_reconstructOrbital(KSO_C,rKSO);
    E                   =   d_C.DFT_energyKinetic(DFT_C.occ,KSO_C);
    E_test              =  [real(sum( conj(rKSO.KSO(:,1)) .* ifft(d_C.Tv.*fft(rKSO.KSO(:,1)))) ), ...
                            real(sum( conj(rKSO.KSO(:,2)) .* ifft(d_C.Tv.*fft(rKSO.KSO(:,2)))) ), ...
                            real(sum( conj(rKSO.KSO(:,3)) .* ifft(d_C.Tv.*fft(rKSO.KSO(:,3)))) )];
    E_test              =   sum(DFT_C.occ.*E_test) * d_C.dx;
    obj.showResult('DFT_energyKinetic (velocity gauge, complex KSO, complex basis)',isreal(E) && abs(E-E_test) < 1e-10);

    d_R.setTv([]);          d_C.setTv([]);

    % NOTE: if this test fails, check the following methods first
    %       - DFT_allocateOrbital
    %       - DFT_reconstructOrbital

    % DFT_energyOrbital ===================================================
    d_R.DFT_allocateOrbital(3,KSO_R,randStr,true);
    d_R.DFT_allocateOrbital(3,KSO_C,randStr,false);

    % Length gauge
    d_R.DFT_reconstructOrbital(KSO_R,rKSO);
    V                   =   DFT_R.Vext.getPotential;        V.initialize(d_R);
    [E,err]             =   d_R.DFT_energyOrbital(V,KSO_R);
    
    Hp                  =  [real(ifft(d_R.T.*fft(rKSO.KSO(:,1)))) + V.V.*rKSO.KSO(:,1), ...
                            real(ifft(d_R.T.*fft(rKSO.KSO(:,2)))) + V.V.*rKSO.KSO(:,2), ...
                            real(ifft(d_R.T.*fft(rKSO.KSO(:,3)))) + V.V.*rKSO.KSO(:,3)];
    E_test              =  [sum(Hp(:,1).*rKSO.KSO(:,1)); ...
                            sum(Hp(:,2).*rKSO.KSO(:,2)); ...
                            sum(Hp(:,3).*rKSO.KSO(:,3))]*d_R.dx;

    rKSO.set('KSO',Hp);     pKSO    =   d_R.DFT_projectOrbital(rKSO);           % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pKSO.KSO(:,1)-E_test(1)*KSO_R.KSO(:,1)).^2); ...
                                  sum(abs(pKSO.KSO(:,2)-E_test(2)*KSO_R.KSO(:,2)).^2); ...
                                  sum(abs(pKSO.KSO(:,3)-E_test(3)*KSO_R.KSO(:,3)).^2)]);
    obj.showResult('DFT_energyOrbital (real orbitals, real basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    
    d_C.DFT_reconstructOrbital(KSO_R,rKSO);
    DFT_C.Vext.getPotential(V);     V.initialize(d_C);
    [E,err]             =   d_C.DFT_energyOrbital(V,KSO_R);
    
    Hp                  =  [ifft(d_C.T.*fft(rKSO.KSO(:,1))) + V.V.*rKSO.KSO(:,1), ...
                            ifft(d_C.T.*fft(rKSO.KSO(:,2))) + V.V.*rKSO.KSO(:,2), ...
                            ifft(d_C.T.*fft(rKSO.KSO(:,3))) + V.V.*rKSO.KSO(:,3)];
    E_test              =  [sum(Hp(:,1).*conj(rKSO.KSO(:,1))); ...
                            sum(Hp(:,2).*conj(rKSO.KSO(:,2))); ...
                            sum(Hp(:,3).*conj(rKSO.KSO(:,3)))]*d_C.dx;

    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pKSO.KSO(:,1)-E_test(1)*KSO_R.KSO(:,1)).^2); ...
                                  sum(abs(pKSO.KSO(:,2)-E_test(2)*KSO_R.KSO(:,2)).^2); ...
                                  sum(abs(pKSO.KSO(:,3)-E_test(3)*KSO_R.KSO(:,3)).^2)]);
    obj.showResult('DFT_energyOrbital (real orbitals, complex basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));


    d_R.DFT_reconstructOrbital(KSO_C,rKSO);
    DFT_R.Vext.getPotential(V);     V.initialize(d_R);
    [E,err]             =   d_R.DFT_energyOrbital(V,KSO_C);
    
    Hp                  =  [ifft(d_R.T.*fft(rKSO.KSO(:,1))) + V.V.*rKSO.KSO(:,1), ...
                            ifft(d_R.T.*fft(rKSO.KSO(:,2))) + V.V.*rKSO.KSO(:,2), ...
                            ifft(d_R.T.*fft(rKSO.KSO(:,3))) + V.V.*rKSO.KSO(:,3)];
    E_test              =  [sum(Hp(:,1).*conj(rKSO.KSO(:,1))); ...
                            sum(Hp(:,2).*conj(rKSO.KSO(:,2))); ...
                            sum(Hp(:,3).*conj(rKSO.KSO(:,3)))]*d_R.dx;

    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pKSO.KSO(:,1)-E_test(1)*KSO_C.KSO(:,1)).^2); ...
                                  sum(abs(pKSO.KSO(:,2)-E_test(2)*KSO_C.KSO(:,2)).^2); ...
                                  sum(abs(pKSO.KSO(:,3)-E_test(3)*KSO_C.KSO(:,3)).^2)]);
    obj.showResult('DFT_energyOrbital (complex orbitals, real basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    
    d_C.DFT_reconstructOrbital(KSO_C,rKSO);
    DFT_C.Vext.getPotential(V);     V.initialize(d_C);
    [E,err]             =   d_C.DFT_energyOrbital(V,KSO_C);
    
    Hp                  =  [ifft(d_C.T.*fft(rKSO.KSO(:,1))) + V.V.*rKSO.KSO(:,1), ...
                            ifft(d_C.T.*fft(rKSO.KSO(:,2))) + V.V.*rKSO.KSO(:,2), ...
                            ifft(d_C.T.*fft(rKSO.KSO(:,3))) + V.V.*rKSO.KSO(:,3)];
    E_test              =  [sum(Hp(:,1).*conj(rKSO.KSO(:,1))); ...
                            sum(Hp(:,2).*conj(rKSO.KSO(:,2))); ...
                            sum(Hp(:,3).*conj(rKSO.KSO(:,3)))]*d_C.dx;

    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pKSO.KSO(:,1)-E_test(1)*KSO_C.KSO(:,1)).^2); ...
                                  sum(abs(pKSO.KSO(:,2)-E_test(2)*KSO_C.KSO(:,2)).^2); ...
                                  sum(abs(pKSO.KSO(:,3)-E_test(3)*KSO_C.KSO(:,3)).^2)]);
    obj.showResult('DFT_energyOrbital (complex orbitals, complex basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    % Velocity gauge
    d_R.setTv(A);           d_C.setTv(A);

    d_R.DFT_reconstructOrbital(KSO_R,rKSO);
    DFT_R.Vext.getPotential(V);     V.initialize(d_R);
    [E,err]             =   d_R.DFT_energyOrbital(V,KSO_R);
    
    Hp                  =  [ifft(d_R.Tv.*fft(rKSO.KSO(:,1))) + V.V.*rKSO.KSO(:,1), ...
                            ifft(d_R.Tv.*fft(rKSO.KSO(:,2))) + V.V.*rKSO.KSO(:,2), ...
                            ifft(d_R.Tv.*fft(rKSO.KSO(:,3))) + V.V.*rKSO.KSO(:,3)];
    E_test              =  [sum(Hp(:,1).*conj(rKSO.KSO(:,1))); ...
                            sum(Hp(:,2).*conj(rKSO.KSO(:,2))); ...
                            sum(Hp(:,3).*conj(rKSO.KSO(:,3)))]*d_R.dx;

    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pKSO.KSO(:,1)-E_test(1)*KSO_R.KSO(:,1)).^2); ...
                                  sum(abs(pKSO.KSO(:,2)-E_test(2)*KSO_R.KSO(:,2)).^2); ...
                                  sum(abs(pKSO.KSO(:,3)-E_test(3)*KSO_R.KSO(:,3)).^2)]);
    obj.showResult('DFT_energyOrbital (velocity gauge, real KSO, real basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    
    d_C.DFT_reconstructOrbital(KSO_R,rKSO);
    DFT_C.Vext.getPotential(V);     V.initialize(d_C);
    [E,err]             =   d_C.DFT_energyOrbital(V,KSO_R);
    
    Hp                  =  [ifft(d_C.Tv.*fft(rKSO.KSO(:,1))) + V.V.*rKSO.KSO(:,1), ...
                            ifft(d_C.Tv.*fft(rKSO.KSO(:,2))) + V.V.*rKSO.KSO(:,2), ...
                            ifft(d_C.Tv.*fft(rKSO.KSO(:,3))) + V.V.*rKSO.KSO(:,3)];
    E_test              =  [sum(Hp(:,1).*conj(rKSO.KSO(:,1))); ...
                            sum(Hp(:,2).*conj(rKSO.KSO(:,2))); ...
                            sum(Hp(:,3).*conj(rKSO.KSO(:,3)))]*d_C.dx;

    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pKSO.KSO(:,1)-E_test(1)*KSO_R.KSO(:,1)).^2); ...
                                  sum(abs(pKSO.KSO(:,2)-E_test(2)*KSO_R.KSO(:,2)).^2); ...
                                  sum(abs(pKSO.KSO(:,3)-E_test(3)*KSO_R.KSO(:,3)).^2)]);
    obj.showResult('DFT_energyOrbital (velocity gauge, real KSO, complex basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));


    d_R.DFT_reconstructOrbital(KSO_C,rKSO);
    DFT_R.Vext.getPotential(V);     V.initialize(d_R);
    [E,err]             =   d_R.DFT_energyOrbital(V,KSO_C);
    
    Hp                  =  [ifft(d_R.Tv.*fft(rKSO.KSO(:,1))) + V.V.*rKSO.KSO(:,1), ...
                            ifft(d_R.Tv.*fft(rKSO.KSO(:,2))) + V.V.*rKSO.KSO(:,2), ...
                            ifft(d_R.Tv.*fft(rKSO.KSO(:,3))) + V.V.*rKSO.KSO(:,3)];
    E_test              =  [sum(Hp(:,1).*conj(rKSO.KSO(:,1))); ...
                            sum(Hp(:,2).*conj(rKSO.KSO(:,2))); ...
                            sum(Hp(:,3).*conj(rKSO.KSO(:,3)))]*d_R.dx;

    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pKSO.KSO(:,1)-E_test(1)*KSO_C.KSO(:,1)).^2); ...
                                  sum(abs(pKSO.KSO(:,2)-E_test(2)*KSO_C.KSO(:,2)).^2); ...
                                  sum(abs(pKSO.KSO(:,3)-E_test(3)*KSO_C.KSO(:,3)).^2)]);
    obj.showResult('DFT_energyOrbital (velocity gauge, complex KSO, real basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));

    
    d_C.DFT_reconstructOrbital(KSO_C,rKSO);
    DFT_C.Vext.getPotential(V);     V.initialize(d_C);
    [E,err]             =   d_C.DFT_energyOrbital(V,KSO_C);
    
    Hp                  =  [ifft(d_C.Tv.*fft(rKSO.KSO(:,1))) + V.V.*rKSO.KSO(:,1), ...
                            ifft(d_C.Tv.*fft(rKSO.KSO(:,2))) + V.V.*rKSO.KSO(:,2), ...
                            ifft(d_C.Tv.*fft(rKSO.KSO(:,3))) + V.V.*rKSO.KSO(:,3)];
    E_test              =  [sum(Hp(:,1).*conj(rKSO.KSO(:,1))); ...
                            sum(Hp(:,2).*conj(rKSO.KSO(:,2))); ...
                            sum(Hp(:,3).*conj(rKSO.KSO(:,3)))]*d_C.dx;

    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    err_test            =   sqrt([sum(abs(pKSO.KSO(:,1)-E_test(1)*KSO_C.KSO(:,1)).^2); ...
                                  sum(abs(pKSO.KSO(:,2)-E_test(2)*KSO_C.KSO(:,2)).^2); ...
                                  sum(abs(pKSO.KSO(:,3)-E_test(3)*KSO_C.KSO(:,3)).^2)]);
    obj.showResult('DFT_energyOrbital (velocity gauge, complex KSO, complex basis)',...
                            all([abs(E-E_test), abs(err-err_test)] < 1e-10,'all'));


    d_R.setTv([]);          d_C.setTv([]);

    % NOTE: if this test fails, check the following methods first
    %       - DFT_allocateOrbital
    %       - DFT_reconstructOrbital
    %       - DFT_projectOrbital

    % DFT_operatorHamiltonian =============================================
    d_R.DFT_allocateOrbital(1,KSO_R,randStr,true);
    d_R.DFT_allocateOrbital(1,KSO_C,randStr,false);

    % Length gauge
    DFT_R.Vext.getPotential(V);

    d_R.DFT_reconstructOrbital(KSO_R,rKSO);     V.initialize(d_R);
    Hp                  =   ifft(d_R.T.*fft(rKSO.KSO)) + V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_R.DFT_operatorHamiltonian(V,KSO_R.KSO);
    obj.showResult('DFT_operatorHamiltonian (real orbitals, real basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_C.DFT_reconstructOrbital(KSO_R,rKSO);     V.initialize(d_C);
    Hp                  =   ifft(d_C.T.*fft(rKSO.KSO)) + V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_C.DFT_operatorHamiltonian(V,KSO_R.KSO);
    obj.showResult('DFT_operatorHamiltonian (real orbitals, complex basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));


    d_R.DFT_reconstructOrbital(KSO_C,rKSO);     V.initialize(d_R);
    Hp                  =   ifft(d_R.T.*fft(rKSO.KSO)) + V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_R.DFT_operatorHamiltonian(V,KSO_C.KSO);
    obj.showResult('DFT_operatorHamiltonian (complex orbitals, real basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_C.DFT_reconstructOrbital(KSO_C,rKSO);     V.initialize(d_C);
    Hp                  =   ifft(d_C.T.*fft(rKSO.KSO)) + V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_C.DFT_operatorHamiltonian(V,KSO_C.KSO);
    obj.showResult('DFT_operatorHamiltonian (complex orbitals, complex basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    % Velocity gauge
    d_R.setTv(A);           d_C.setTv(A);

    d_R.DFT_reconstructOrbital(KSO_R,rKSO);     V.initialize(d_R);
    Hp                  =   ifft(d_R.Tv.*fft(rKSO.KSO)) + V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_R.DFT_operatorHamiltonian(V,KSO_R.KSO);
    obj.showResult('DFT_operatorHamiltonian (vel. g., real KSO, real basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_C.DFT_reconstructOrbital(KSO_R,rKSO);     V.initialize(d_C);
    Hp                  =   ifft(d_C.Tv.*fft(rKSO.KSO)) + V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_C.DFT_operatorHamiltonian(V,KSO_R.KSO);
    obj.showResult('DFT_operatorHamiltonian (vel. g., real KSO, complex basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));


    d_R.DFT_reconstructOrbital(KSO_C,rKSO);     V.initialize(d_R);
    Hp                  =   ifft(d_R.Tv.*fft(rKSO.KSO)) + V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_R.DFT_operatorHamiltonian(V,KSO_C.KSO);
    obj.showResult('DFT_operatorHamiltonian (vel. g., complex KSO, real basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_C.DFT_reconstructOrbital(KSO_C,rKSO);     V.initialize(d_C);
    Hp                  =   ifft(d_C.Tv.*fft(rKSO.KSO)) + V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_C.DFT_operatorHamiltonian(V,KSO_C.KSO);
    obj.showResult('DFT_operatorHamiltonian (vel. g., complex KSO, complex basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_R.setTv([]);          d_C.setTv([]);

    % NOTE: if this test fails, check the following methods first
    %       - DFT_allocateOrbital
    %       - DFT_reconstructOrbital
    %       - DFT_projectOrbital

    % DFT_operatorKinetic =================================================
    d_R.DFT_allocateOrbital(1,KSO_R,randStr,true);
    d_R.DFT_allocateOrbital(1,KSO_C,randStr,false);

    % Length gauge
    d_R.DFT_reconstructOrbital(KSO_R,rKSO);
    Hp                  =   ifft(d_R.T.*fft(rKSO.KSO));
    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_R.DFT_operatorKinetic(KSO_R.KSO);
    obj.showResult('DFT_operatorKinetic (real orbitals, real basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_C.DFT_reconstructOrbital(KSO_R,rKSO);
    Hp                  =   ifft(d_C.T.*fft(rKSO.KSO));
    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_C.DFT_operatorKinetic(KSO_R.KSO);
    obj.showResult('DFT_operatorKinetic (real orbitals, complex basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));


    d_R.DFT_reconstructOrbital(KSO_C,rKSO);
    Hp                  =   ifft(d_R.T.*fft(rKSO.KSO));
    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_R.DFT_operatorKinetic(KSO_C.KSO);
    obj.showResult('DFT_operatorKinetic (complex orbitals, real basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_C.DFT_reconstructOrbital(KSO_C,rKSO);
    Hp                  =   ifft(d_C.T.*fft(rKSO.KSO));
    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_C.DFT_operatorKinetic(KSO_C.KSO);
    obj.showResult('DFT_operatorKinetic (complex orbitals, complex basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    % Velocity gauge
    d_R.setTv(A);           d_C.setTv(A);

    d_R.DFT_reconstructOrbital(KSO_R,rKSO);
    Hp                  =   ifft(d_R.Tv.*fft(rKSO.KSO));
    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_R.DFT_operatorKinetic(KSO_R.KSO);
    obj.showResult('DFT_operatorKinetic (vel. g., real KSO, real basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_C.DFT_reconstructOrbital(KSO_R,rKSO);
    Hp                  =   ifft(d_C.Tv.*fft(rKSO.KSO));
    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_C.DFT_operatorKinetic(KSO_R.KSO);
    obj.showResult('DFT_operatorKinetic (vel. g., real KSO, complex basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));


    d_R.DFT_reconstructOrbital(KSO_C,rKSO);
    Hp                  =   ifft(d_R.Tv.*fft(rKSO.KSO));
    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_R.DFT_operatorKinetic(KSO_C.KSO);
    obj.showResult('DFT_operatorKinetic (vel. g., complex KSO, real basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_C.DFT_reconstructOrbital(KSO_C,rKSO);
    Hp                  =   ifft(d_C.Tv.*fft(rKSO.KSO));
    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_C.DFT_operatorKinetic(KSO_C.KSO);
    obj.showResult('DFT_operatorKinetic (vel. g., complex KSO, complex basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_R.setTv([]);          d_C.setTv([]);

    % NOTE: if this test fails, check the following methods first
    %       - DFT_allocateOrbital
    %       - DFT_reconstructOrbital
    %       - DFT_projectOrbital


    % DFT_operatorPotential ===============================================
    d_R.DFT_allocateOrbital(1,KSO_R,randStr,true);
    d_R.DFT_allocateOrbital(1,KSO_C,randStr,false);

    % Length gauge
    DFT_R.Vext.getPotential(V);

    d_R.DFT_reconstructOrbital(KSO_R,rKSO);     V.initialize(d_R);
    Hp                  =   V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_R.DFT_operatorPotential(V,KSO_R.KSO);
    obj.showResult('DFT_operatorPotential (real orbitals, real basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_C.DFT_reconstructOrbital(KSO_R,rKSO);     V.initialize(d_C);
    Hp                  =   V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_C.DFT_operatorPotential(V,KSO_R.KSO);
    obj.showResult('DFT_operatorPotential (real orbitals, complex basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));


    d_R.DFT_reconstructOrbital(KSO_C,rKSO);     V.initialize(d_R);
    Hp                  =   V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_R.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_R.DFT_operatorPotential(V,KSO_C.KSO);
    obj.showResult('DFT_operatorPotential (complex orbitals, real basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    d_C.DFT_reconstructOrbital(KSO_C,rKSO);     V.initialize(d_C);
    Hp                  =   V.V.*rKSO.KSO;
    rKSO.set('KSO',Hp);     d_C.DFT_projectOrbital(rKSO,pKSO);              % Project Hp onto basis
    Hp                  =   d_C.DFT_operatorPotential(V,KSO_C.KSO);
    obj.showResult('DFT_operatorPotential (complex orbitals, complex basis)',...
                            all(abs(Hp - pKSO.KSO) < 1e-10,'all'));

    % NOTE: if this test fails, check the following methods first
    %       - DFT_allocateOrbital
    %       - DFT_reconstructOrbital
    %       - DFT_projectOrbital

end