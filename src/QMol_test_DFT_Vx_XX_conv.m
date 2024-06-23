classdef QMol_test_DFT_Vx_XX_conv < QMol_test
%QMol_test_DFT_Vx_XX_conv unit test for the QMol_DFT_Vx_XX_conv class

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

% WARNING: dipole acceleration is commented out to avoid warnings (from the
%          parent object) during the test. Put it back once test is
%          available.

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_test)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test_DFT_Vx_XX_conv\n'); 
    QMol_test_DFT_Vx_XX_conv.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test units
    obj.runSpinRestricted;
    obj.runSpinPolarized;
    obj.runSpinRestricted_basis;
    obj.runSpinPolarized_basis;
    
end
end
methods (Access=private)
function runSpinRestricted(obj) %==========================================
%runSpinRestricted unit tests for spin-restricted DFT
    
    % Initialization
    obj.showSection('Spin-restricted DFT');

    s_ee                =   pi;
    Vee                 =   @(x)    exp(-x.^2 * .5/s_ee^2);
    
    Vx                  =   QMol_DFT_Vx_XX_conv('interactionPotential',Vee);
    DFT                 =   QMol_DFT_spinRes(...
                                'xspan',                -15:.3:20,          ...
                                'occupation',           [2 1.5 1.3],        ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vx);
    DFT.initialize;

    % Kohn-Sham orbitals
    x                   =   DFT.x(:);   dx  =   x(2)-x(1);
    s                   =   2;
    p                   =  [exp(-x.^2*.25/s^2),x.*exp(-x.^2*.25/s^2),x.^2.*exp(-x.^2*.25/s^2)];
    p(:,1)              =   p(:,1) / sqrt(sum(p(:,1).^2)*dx);
    p(:,2)              =   p(:,2) / sqrt(sum(p(:,2).^2)*dx);
    p(:,3)              =   p(:,3) / sqrt(sum(p(:,3).^2)*dx);

    DFT.KSO.KSO         =   p;
    
    DFT.setPotentialKernel;

    % Exact-exchange potential
    v                   =   exp(-(x-2).^2*.25/1.5^2);
    v                   =   v / sqrt(sum(v.^2)*dx);

    Hv                  =   -.5*DFT.occ(1) * p(:,1) .* conv(p(:,1).*v,Vx.V,'same') * dx ...
                            -.5*DFT.occ(2) * p(:,2) .* conv(p(:,2).*v,Vx.V,'same') * dx ...
                            -.5*DFT.occ(3) * p(:,3) .* conv(p(:,3).*v,Vx.V,'same') * dx;
    
    R                   =   max(abs(Hv - Vx.applyPotential(v))) < 1e-10;
    obj.showResult('applyPotential' ,R);

    % Exact-exchange potential derivative (dipole acceleration)
% COMMENTED OUT TO AVOID WARNING IN RELEASE VERSION
%     p(:,1)              =   p(:,1) .* exp(.2i*DFT.x(:).^2);
%     p(:,2)              =   p(:,2) .* exp(.5i*DFT.x(:));
%     p(:,3)              =   p(:,3) .* exp(.2i*DFT.x(:)+1i);
% 
%     DFT.KSO.KSO         =   p;
%     DFT.setPotentialKernel;
% 
%     a                   =   DFT.occ(1)*Vx.applyPotentialDerivative('dipacc',p(:,1)) + ...
%                             DFT.occ(2)*Vx.applyPotentialDerivative('dipacc',p(:,2)) + ...
%                             DFT.occ(3)*Vx.applyPotentialDerivative('dipacc',p(:,3));
%     R                   =   isreal(a) && abs(a) < 1e-10;
%     obj.showResult('applyPotentialDerivative (dipole acceleration)' ,R);
    fprintf('    > applyPotentialDerivative (dipole acceleration) is untested      ****\n');

    % Exact-exchange energy
    E                   =   0;
    for k = 1:3, for l = 1:3                                                %#ok<ALIGN> 
        E               =   E - .25*DFT.occ(k)*DFT.occ(l) * sum( conj(p(:,k)).*p(:,l) .* conv(p(:,k).*conj(p(:,l)),Vx.V,'same') ) * dx^2;
    end, end
    
    R                   =   abs(E - Vx.getEnergy(v)) < 1e-10;
    obj.showResult('getEnergy' ,R);

end
function runSpinRestricted_basis(obj) %====================================
%runSpinRestricted unit tests for spin-restricted DFT
    
    % Initialization
    obj.showSection('Basis-set spin-restricted DFT');

    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];
    disc                =   QMol_disc_basis('xspan',x,'basis',v);   disc.orthonormalizeBasis;

    s_ee                =   pi;
    Vee                 =   @(x)    exp(-x.^2 * .5/s_ee^2);
    
    Vx                  =   QMol_DFT_Vx_XX_conv('interactionPotential',Vee);
    DFT                 =   QMol_DFT_spinRes(...
                                'disc',                 disc,               ...
                                'occupation',           [2 1.5 1.3],        ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vx);
    DFT.initialize;

    % Kohn-Sham orbitals
    p                   =   disc.DFT_normalizeOrbital(rand(randStr,[disc.basisSize,3])-.5);
    DFT.KSO.KSO         =   p;
    P                   =   disc.DFT_reconstructOrbital(DFT.KSO);
    
    DFT.setPotentialKernel;

    % Exact-exchange potential
    v                   =   disc.DFT_normalizeOrbital(rand(randStr,[disc.basisSize,1])-.5);
    V                   =   v(1)*disc.v(:,1) + v(2)*disc.v(:,2) + v(3)*disc.v(:,3);

    Hv                  =   -.5*DFT.occ(1) * P.KSO(:,1) .* conv(P.KSO(:,1).*V,Vx.V,'same') * disc.dx ...
                            -.5*DFT.occ(2) * P.KSO(:,2) .* conv(P.KSO(:,2).*V,Vx.V,'same') * disc.dx ...
                            -.5*DFT.occ(3) * P.KSO(:,3) .* conv(P.KSO(:,3).*V,Vx.V,'same') * disc.dx;

    R                   =   max(abs(Hv - Vx.applyPotential(V))) < 1e-10;
    obj.showResult('applyPotential (real basis)' ,R);

    % Exact-exchange energy
    E                   =   0;
    for k = 1:3, for l = 1:3                                                %#ok<ALIGN> 
        E               =   E - .25*DFT.occ(k)*DFT.occ(l) * sum( P.KSO(:,k).*P.KSO(:,l) .* conv(P.KSO(:,k).*P.KSO(:,l),Vx.V,'same') ) * disc.dx^2;
    end, end

    R                   =   abs(E - Vx.getEnergy(V)) < 1e-10;
    obj.showResult('getEnergy (real basis)' ,R);

    % Complex basis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)] .* ...
                            [exp( 2i*x),    exp( 3i*x),       exp(-4i*x),exp(-2i*x)];
    disc.set('basis',v);    disc.orthonormalizeBasis;
    DFT.reset();            DFT.initialize;

    % Kohn-Sham orbitals
    p                   =   disc.DFT_normalizeOrbital( ...
                                rand(randStr,[disc.basisSize,3])-.5 + ...
                                1i*(rand(randStr,[disc.basisSize,3])-.5) );
    DFT.KSO.KSO         =   p;
    P                   =   disc.DFT_reconstructOrbital(DFT.KSO);
    
    DFT.setPotentialKernel;

    % Exact-exchange potential
    v                   =   disc.DFT_normalizeOrbital(rand(randStr,[disc.basisSize,1])-.5);
    V                   =   v(1)*disc.v(:,1) + v(2)*disc.v(:,2) + v(3)*disc.v(:,3);

    Hv                  =   -.5*DFT.occ(1) * P.KSO(:,1) .* conv(conj(P.KSO(:,1)).*V,Vx.V,'same') * disc.dx ...
                            -.5*DFT.occ(2) * P.KSO(:,2) .* conv(conj(P.KSO(:,2)).*V,Vx.V,'same') * disc.dx ...
                            -.5*DFT.occ(3) * P.KSO(:,3) .* conv(conj(P.KSO(:,3)).*V,Vx.V,'same') * disc.dx;

    R                   =   max(abs(Hv - Vx.applyPotential(V))) < 1e-10;
    obj.showResult('applyPotential (complex basis)' ,R);

    % Exact-exchange energy
    E                   =   0;
    for k = 1:3, for l = 1:3                                                %#ok<ALIGN> 
        E               =   E - .25*DFT.occ(k)*DFT.occ(l) * sum( ...
                                conj(P.KSO(:,k)).*P.KSO(:,l) .* ...
                                conv(P.KSO(:,k).*conj(P.KSO(:,l)),Vx.V,'same') ) * disc.dx^2;
    end, end

    R                   =   abs(E - Vx.getEnergy(V)) < 1e-10;
    obj.showResult('getEnergy (complex basis)' ,R);

end
function runSpinPolarized(obj) %===========================================
%runSpinRestricted unit tests for spin-restricted DFT
    
    % Initialization
    obj.showSection('Spin-polarized DFT');

    s_ee                =   pi;
    Vee                 =   @(x)    exp(-x.^2 * .5/s_ee^2);
    
    Vx                  =   QMol_DFT_Vx_XX_conv('interactionPotential',Vee);
    DFT                 =   QMol_DFT_spinPol(...
                                'xspan',                -15:.3:20,          ...
                                'occupation',           {[1 .9],[.8 .7 .1]}, ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vx);
    DFT.initialize;
    
    % Kohn-Sham orbitals
    x                   =   DFT.x(:);   dx  =   x(2)-x(1);
    s                   =  [2 1.5];
    p                   =  [exp(-x.^2*.25/s(1)^2),x.*exp(-x.^2*.25/s(1)^2),x.^2.*exp(-x.^2*.25/s(1)^2),...
                            exp(-x.^2*.25/s(2)^2),x.*exp(-x.^2*.25/s(2)^2),];
    p(:,1)              =   p(:,1) / sqrt(sum(p(:,1).^2)*dx);
    p(:,2)              =   p(:,2) / sqrt(sum(p(:,2).^2)*dx);
    p(:,3)              =   p(:,3) / sqrt(sum(p(:,3).^2)*dx);
    p(:,4)              =   p(:,4) / sqrt(sum(p(:,4).^2)*dx);
    p(:,5)              =   p(:,5) / sqrt(sum(p(:,5).^2)*dx);

    DFT.KSO.KSOup       =   p(:,1:2);
    DFT.KSO.KSOdw       =   p(:,3:5);
    
    DFT.setPotentialKernel;

    % Exact-exchange potential
    v                   =   exp(-(x-2).^2*.25/1.5^2);
    v                   =   v / sqrt(sum(v.^2)*dx);

    Hv                  =   -DFT.occ{1}(1) * p(:,1) .* conv(p(:,1).*v,Vx.V,'same') * dx ...
                            -DFT.occ{1}(2) * p(:,2) .* conv(p(:,2).*v,Vx.V,'same') * dx;
    
    R                   =   max(abs(Hv - Vx.applyPotential(v,true))) < 1e-10;
    obj.showResult('applyPotential (spin up)' ,R);

    Hv                  =   -DFT.occ{2}(1) * p(:,3) .* conv(p(:,3).*v,Vx.V,'same') * dx ...
                            -DFT.occ{2}(2) * p(:,4) .* conv(p(:,4).*v,Vx.V,'same') * dx ...
                            -DFT.occ{2}(3) * p(:,5) .* conv(p(:,5).*v,Vx.V,'same') * dx;
    
    R                   =   max(abs(Hv - Vx.applyPotential(v,false))) < 1e-10;
    obj.showResult('applyPotential (spin down)' ,R);

    % Exact-exchange potential derivative (dipole acceleration)
% COMMENTED OUT TO AVOID WARNING IN RELEASE VERSION
%     p(:,1)              =   p(:,1) .* exp(.2i*DFT.x(:).^2);
%     p(:,2)              =   p(:,2) .* exp(.5i*DFT.x(:));
%     p(:,3)              =   p(:,3) .* exp(.2i*DFT.x(:)+1i);
%     p(:,4)              =   p(:,4) .* exp(.2i*DFT.x(:).^2);
%     p(:,5)              =   p(:,5) .* exp(.5i*DFT.x(:));
% 
%     DFT.KSO.KSOup       =   p(:,1:2);
%     DFT.KSO.KSOdw       =   p(:,3:5);
%     DFT.setPotentialKernel;
% 
%     a_up                =   DFT.occ{1}(1)*Vx.applyPotentialDerivative('dipacc',p(:,1),true) + ...
%                             DFT.occ{1}(2)*Vx.applyPotentialDerivative('dipacc',p(:,2),true);
%     a_dw                =   DFT.occ{2}(1)*Vx.applyPotentialDerivative('dipacc',p(:,3),false) + ...
%                             DFT.occ{2}(2)*Vx.applyPotentialDerivative('dipacc',p(:,4),false) + ...
%                             DFT.occ{2}(3)*Vx.applyPotentialDerivative('dipacc',p(:,5),false);
%     R                   =   isreal(a_up) && isreal(a_dw) && abs(a_up) < 1e-10 && abs(a_dw) < 1e-10;
%     obj.showResult('applyPotentialDerivative (dipole acceleration)' ,R);
    fprintf('    > applyPotentialDerivative (dipole acceleration) is untested      ****\n');

    % Exact-exchange energy
    E                   =   0;
    for k = 1:2, for l = 1:2                                                %#ok<ALIGN> 
        E               =   E - .5*DFT.occ{1}(k)*DFT.occ{1}(l) * sum( conj(p(:,k)).*p(:,l) .* conv(p(:,k).*conj(p(:,l)),Vx.V,'same') ) * dx^2;
    end, end
    for k = 3:5, for l = 3:5                                                %#ok<ALIGN> 
        E               =   E - .5*DFT.occ{2}(k-2)*DFT.occ{2}(l-2) * sum( conj(p(:,k)).*p(:,l) .* conv(p(:,k).*conj(p(:,l)),Vx.V,'same') ) * dx^2;
    end, end
    
    R                   =   abs(E - Vx.getEnergy(v)) < 1e-10;
    obj.showResult('getEnergy' ,R);


end
function runSpinPolarized_basis(obj) %=====================================
%runSpinRestricted unit tests for spin-restricted DFT
    
    % Initialization
    obj.showSection('Basis-set spin-polarized DFT');

    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];
    disc                =   QMol_disc_basis('xspan',x,'basis',v);   disc.orthonormalizeBasis;

    s_ee                =   pi;
    Vee                 =   @(x)    exp(-x.^2 * .5/s_ee^2);
    
    Vx                  =   QMol_DFT_Vx_XX_conv('interactionPotential',Vee);
    DFT                 =   QMol_DFT_spinPol(...
                                'disc',                 disc,               ...
                                'occupation',           {[1 .9],[.8 .7 .1]},...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'HartreePotential',     QMol_DFT_Vh_conv,   ...
                                'exchangeCorrelationPotential', Vx);
    DFT.initialize;

    % Kohn-Sham orbitals
    p                   =   disc.DFT_normalizeOrbital(rand(randStr,[disc.basisSize,5])-.5);
    DFT.KSO.KSOup       =   p(:,1:2);
    DFT.KSO.KSOdw       =   p(:,3:5);
    P                   =   disc.DFT_reconstructOrbital(DFT.KSO);
    
    DFT.setPotentialKernel;

    % Exact-exchange potential
    v                   =   disc.DFT_normalizeOrbital(rand(randStr,[disc.basisSize,1])-.5);
    V                   =   v(1)*disc.v(:,1) + v(2)*disc.v(:,2) + v(3)*disc.v(:,3);

    Hv                  =   -DFT.occ{1}(1) * P.KSOup(:,1) .* conv(P.KSOup(:,1).*V,Vx.V,'same') * disc.dx ...
                            -DFT.occ{1}(2) * P.KSOup(:,2) .* conv(P.KSOup(:,2).*V,Vx.V,'same') * disc.dx;
    R                   =   max(abs(Hv - Vx.applyPotential(V,true))) < 1e-10;
    obj.showResult('applyPotential (real basis, spin up)' ,R);


    Hv                  =   -DFT.occ{2}(1) * P.KSOdw(:,1) .* conv(P.KSOdw(:,1).*V,Vx.V,'same') * disc.dx ...
                            -DFT.occ{2}(2) * P.KSOdw(:,2) .* conv(P.KSOdw(:,2).*V,Vx.V,'same') * disc.dx ...
                            -DFT.occ{2}(3) * P.KSOdw(:,3) .* conv(P.KSOdw(:,3).*V,Vx.V,'same') * disc.dx;
    R                   =   max(abs(Hv - Vx.applyPotential(V,false))) < 1e-10;
    obj.showResult('applyPotential (real basis, spin down)' ,R);

    % Exact-exchange energy
    E                   =   0;
    for k = 1:2, for l = 1:2                                                %#ok<ALIGN> 
        E               =   E - .5*DFT.occ{1}(k)*DFT.occ{1}(l) * sum( P.KSOup(:,k).*P.KSOup(:,l) .* conv(P.KSOup(:,k).*P.KSOup(:,l),Vx.V,'same') ) * disc.dx^2;
    end, end
    for k = 1:3, for l = 1:3                                                %#ok<ALIGN> 
        E               =   E - .5*DFT.occ{2}(k)*DFT.occ{2}(l) * sum( P.KSOdw(:,k).*P.KSOdw(:,l) .* conv(P.KSOdw(:,k).*P.KSOdw(:,l),Vx.V,'same') ) * disc.dx^2;
    end, end

    R                   =   abs(E - Vx.getEnergy(V)) < 1e-10;
    obj.showResult('getEnergy (real basis)' ,R);

    % Complex basis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)] .* ...
                            [exp( 2i*x),    exp( 3i*x),       exp(-4i*x),exp(-2i*x)];
    disc.set('basis',v);    disc.orthonormalizeBasis;
    DFT.reset();            DFT.initialize;

    % Kohn-Sham orbitals
    p                   =   disc.DFT_normalizeOrbital( ...
                                rand(randStr,[disc.basisSize,5])-.5 + ...
                                1i*(rand(randStr,[disc.basisSize,5])-.5) );
    DFT.KSO.KSOup       =   p(:,1:2);
    DFT.KSO.KSOdw       =   p(:,3:5);
    P                   =   disc.DFT_reconstructOrbital(DFT.KSO);
    
    DFT.setPotentialKernel;

    % Exact-exchange potential
    v                   =   disc.DFT_normalizeOrbital(rand(randStr,[disc.basisSize,1])-.5);
    V                   =   v(1)*disc.v(:,1) + v(2)*disc.v(:,2) + v(3)*disc.v(:,3);

    Hv                  =   -DFT.occ{1}(1) * P.KSOup(:,1) .* conv(conj(P.KSOup(:,1)).*V,Vx.V,'same') * disc.dx ...
                            -DFT.occ{1}(2) * P.KSOup(:,2) .* conv(conj(P.KSOup(:,2)).*V,Vx.V,'same') * disc.dx;
    R                   =   max(abs(Hv - Vx.applyPotential(V,true))) < 1e-10;
    obj.showResult('applyPotential (real basis, spin up)' ,R);


    Hv                  =   -DFT.occ{2}(1) * P.KSOdw(:,1) .* conv(conj(P.KSOdw(:,1)).*V,Vx.V,'same') * disc.dx ...
                            -DFT.occ{2}(2) * P.KSOdw(:,2) .* conv(conj(P.KSOdw(:,2)).*V,Vx.V,'same') * disc.dx ...
                            -DFT.occ{2}(3) * P.KSOdw(:,3) .* conv(conj(P.KSOdw(:,3)).*V,Vx.V,'same') * disc.dx;
    R                   =   max(abs(Hv - Vx.applyPotential(V,false))) < 1e-10;
    obj.showResult('applyPotential (real basis, spin down)' ,R);

    % Exact-exchange energy
    E                   =   0;
    for k = 1:2, for l = 1:2                                                %#ok<ALIGN> 
        E               =   E - .5*DFT.occ{1}(k)*DFT.occ{1}(l) * sum( conj(P.KSOup(:,k)).*P.KSOup(:,l) .* conv(P.KSOup(:,k).*conj(P.KSOup(:,l)),Vx.V,'same') ) * disc.dx^2;
    end, end
    for k = 1:3, for l = 1:3                                                %#ok<ALIGN> 
        E               =   E - .5*DFT.occ{2}(k)*DFT.occ{2}(l) * sum( conj(P.KSOdw(:,k)).*P.KSOdw(:,l) .* conv(P.KSOdw(:,k).*conj(P.KSOdw(:,l)),Vx.V,'same') ) * disc.dx^2;
    end, end

    R                   =   abs(E - Vx.getEnergy(V)) < 1e-10;
    obj.showResult('getEnergy (real basis)' ,R);


end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

