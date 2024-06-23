classdef QMol_test_DFT_density < QMol_test
%QMol_test_DFT_density suite of unit tests for QMol_DFT_density

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
    fprintf('  * QMol_test_DFT_density\n'); 
    QMol_test_DFT_density.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test components
    obj.test_operator;
    obj.test_gradient;
    obj.test_package_methods;
    
end
end
%% Test components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=private)
function test_operator(obj) %==============================================
%test_operator unit tests for operator overload in the class
    
    % Initialization
    obj.showSection('Object comparison');
    
    % Spin restricted
    R                   =   rand(15,1);
    rho1                =   QMol_DFT_density('isSpinPol',false,'rho',R);
    rho2                =   QMol_DFT_density('isSpinPol',false,'rho',R + .5e-10);
    obj.showResult('== (eq operator, spin restricted)' ,rho1 == rho2);

    rho2.clear;     rho2.set('isSpinPol',false,'rho',R + 1.01e-10);
    obj.showResult('~= (value mismatch, spin restricted)' ,rho1 ~= rho2);

    rho2.clear;     rho2.set('isSpinPol',false,'rho',[]);
    obj.showResult('~= (size mismatch, spin restricted)' ,rho1 ~= rho2);

    % Spin polarized
    rho1.clear;     rho1.set('isSpinPol',true,'rhoUp',2*R       ,'rhoDw',R-1);
    rho2.clear;     rho2.set('isSpinPol',true,'rhoUp',2*R-.5e-10,'rhoDw',R-1+.5e-10);
    obj.showResult('== (eq operator, spin polarized)' ,rho1 == rho2);

    rho2.clear;     rho2.set('isSpinPol',true,'rhoUp',2*R-1.01e-10,'rhoDw',R-1+.5e-10);
    obj.showResult('~= (up-spin value mismatch, spin polarized)' ,rho1 ~= rho2);

    rho2.clear;     rho2.set('isSpinPol',true,'rhoUp',2*R-.5e-10,'rhoDw',R-1+1.01e-10);
    obj.showResult('~= (down-spin value mismatch, spin polarized)' ,rho1 ~= rho2);

    rho2.clear;     rho2.set('isSpinPol',true,'rhoUp',[],'rhoDw',R-1+.5e-10);
    obj.showResult('~= (up-spin size mismatch, spin polarized)' ,rho1 ~= rho2);

    rho2.clear;     rho2.set('isSpinPol',true,'rhoUp',2*R-.5e-10,'rhoDw',[]);
    obj.showResult('~= (down-spin size mismatch, spin polarized)' ,rho1 ~= rho2);


    rho1.clear;     rho1.set('isSpinPol',true,'rho',R,'rhoUp',2*R       ,'rhoDw',R-1);
    rho2.clear;     rho2.set('isSpinPol',false,'rho',R,'rhoUp',2*R-.5e-10,'rhoDw',R-1+.5e-10);
    obj.showResult('~= (spin polarized/restricted mismatch)' ,rho1 ~= rho2);

end
function test_gradient(obj) %==============================================
%test_gradient unit tests for the computation of the gradient
    
    % Initialization
    obj.showSection('Gradient computation');

    d_SR                =   QMol_disc('xspan',-15:.1:20);       
    d_SR.initialize(QMol_DFT_spinRes('discretization',d_SR));

    d_SP                =   QMol_disc('xspan',-15:.1:22);       
    d_SP.initialize(QMol_DFT_spinPol('discretization',d_SP));

    % Spin restricted
    rho = QMol_DFT_density('isSpinPol',false,'density',exp(-(d_SR.x(:)-3.5).^2/2/3));
    D_rho               =   -(d_SR.x(:)-3.5).*exp(-(d_SR.x(:)-3.5).^2/2/3)/3;
    rho.initialize(d_SR);

    R                   =   all(abs(rho.getGradient(1) - D_rho) < 1e-10,'all');
    obj.showResult('getGradient (spin restricted)',R);

    rho.set('density',2.*exp(-(d_SR.x(:)-3.5).^2/2/3));
    rho.initialize(d_SR);

    R                   =   all(abs(rho.getGradient(1) - 2*D_rho) < 1e-10,'all');
    obj.showResult('getGradient (after density update, spin restricted)',R);

    % Spin polarized
    rho = QMol_DFT_density('isSpinPol',true, ...
                'densityUp',exp(-(d_SP.x(:)-4.5).^2/2/3), ...
                'densityDown',d_SP.x(:).*exp(-(d_SP.x(:)-4.2).^2/2/3));
    D_rhoUp             =   -(d_SP.x(:)-4.5).*exp(-(d_SP.x(:)-4.5).^2/2/3)/3;
    D_rhoDw             =   -d_SP.x(:).*(d_SP.x(:)-4.2).*exp(-(d_SP.x(:)-4.2).^2/2/3)/3 ...
                            +exp(-(d_SP.x(:)-4.2).^2/2/3);
    rho.initialize(d_SP);

    R                   =   all(abs(rho.getGradient(1,true ) - D_rhoUp) < 1e-10,'all') && ...
                            all(abs(rho.getGradient(1,false) - D_rhoDw) < 1e-10,'all');
    obj.showResult('getGradient (spin polarized)',R);

    rho.set('densityUp',3*exp(-(d_SP.x(:)-4.5).^2/2/3), ...
            'densityDown',4*d_SP.x(:).*exp(-(d_SP.x(:)-4.2).^2/2/3));
    rho.initialize(d_SP);

    R                   =   all(abs(rho.getGradient(1,true ) - 3*D_rhoUp) < 1e-10,'all') && ...
                            all(abs(rho.getGradient(1,false) - 4*D_rhoDw) < 1e-10,'all');
    obj.showResult('getGradient (after density update, spin polarized)',R);

end
function test_package_methods(obj) %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%test_package_methods unit tests for QMol-grid package methods
    
    % Initialization
    obj.showSection('QMol-grid package methods');
    
    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility

    d_SR                =   QMol_disc('xspan',-15:.1:20);
    d_SR.initialize(QMol_DFT_spinRes('discretization',d_SR));

    d_SP                =   QMol_disc('xspan',-15:.1:22);
    d_SP.initialize(QMol_DFT_spinPol('discretization',d_SP));

    x                   =   d_SR.x(:);
    v                   =  [exp(-(x-2).^2),     (x-2).*exp(-(x-2).^2/2),    ...
                            exp(-(x-1).^2/.7),  (x-1).*exp(-(x-1).^2/1.5),  ...
                            exp(-x.^2),         x.*exp(-x.^2),              ...
                            exp(-(x+1).^2/2),   (x+1).*exp(-(x+1).^2/4)];
    d_SR_b              =   QMol_disc_basis('xspan',x.','basis',v);         d_SR_b.orthonormalizeBasis;
    d_SR_b.initialize(QMol_DFT_spinRes('discretization',d_SR_b));

    x                   =   d_SP.x(:);
    v                   =  [exp(-(x-2).^2),     (x-2).*exp(-(x-2).^2/2),    ...
                            exp(-(x-1).^2/.7),  (x-1).*exp(-(x-1).^2/1.5),  ...
                            exp(-x.^2),         x.*exp(-x.^2),              ...
                            exp(-(x+1).^2/2),   (x+1).*exp(-(x+1).^2/4)];
    d_SP_b              =   QMol_disc_basis('xspan',x.','basis',v);         d_SP_b.orthonormalizeBasis;
    d_SP_b.initialize(QMol_DFT_spinPol('discretization',d_SP_b));

    % Density operations (spin restricted) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    occ                 =   [1 2 pi/2 1/3 .5];
    KSO                 =   d_SR.DFT_allocateOrbital(5,[],randStr,true);
    rho                 =   d_SR.DFT_allocateDensity;
    
    R                   =   occ(1)*KSO.KSO(:,1).^2 + occ(2)*KSO.KSO(:,2).^2 + ...
                            occ(3)*KSO.KSO(:,3).^2 + occ(4)*KSO.KSO(:,4).^2 + ...
                            occ(5)*KSO.KSO(:,5).^2;
    rho.getDensity(occ,KSO);
    R                   =   max(abs(rho.density-R),[],'all') < 1e-10;
    obj.showResult('getDensity (spin restricted)' ,R && rho.isInit);

    N                   =   rho.getCharge;
    R                   =   abs(N-sum(rho.rho)*(d_SR.x(2)-d_SR.x(1))) < 1e-10;
    obj.showResult('getCharge (spin restricted)' ,R);

    % Density operations (spin polarized) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    occ                 =   {[1 2 pi/2 1/3 .5],[exp(1) log(5) 1]};
    d_SP.DFT_allocateOrbital([5 3],KSO,randStr,true);
    
    Rup                 =   occ{1}(1)*KSO.KSOup(:,1).^2 + occ{1}(2)*KSO.KSOup(:,2).^2 + ...
                            occ{1}(3)*KSO.KSOup(:,3).^2 + occ{1}(4)*KSO.KSOup(:,4).^2 + ...
                            occ{1}(5)*KSO.KSOup(:,5).^2;
    Rdw                 =   occ{2}(1)*KSO.KSOdw(:,1).^2 + occ{2}(2)*KSO.KSOdw(:,2).^2 + ...
                            occ{2}(3)*KSO.KSOdw(:,3).^2;
    rho.getDensity(occ,KSO);
    R                   =   max(abs(rho.densityUp  -Rup),[],'all') < 1e-10  && ...
                            max(abs(rho.densityDown-Rdw),[],'all') < 1e-10;
    obj.showResult('getDensity (spin polarized)' ,R && rho.isInit);

    N                   =   rho.getCharge;
    R                   =   abs(N(1)-sum(rho.rhoUp)*(d_SP.x(2)-d_SR.x(1))) < 1e-10 && ...
                            abs(N(2)-sum(rho.rhoDw)*(d_SP.x(2)-d_SR.x(1))) < 1e-10;
    obj.showResult('getCharge (spin polarized)' ,R);

    % Density operations (basis-set model) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    occ                 =   [1 2 pi/2 1/3 .5];
    proj                =   d_SR_b.DFT_allocateOrbital(5,[],randStr,true);
    KSO                 =   d_SR_b.DFT_reconstructOrbital(proj);
    rho                 =   d_SR_b.DFT_allocateDensity;
    
    R                   =   occ(1)*KSO.KSO(:,1).^2 + occ(2)*KSO.KSO(:,2).^2 + ...
                            occ(3)*KSO.KSO(:,3).^2 + occ(4)*KSO.KSO(:,4).^2 + ...
                            occ(5)*KSO.KSO(:,5).^2;
    rho.getDensity(occ,proj);
    R                   =   max(abs(rho.density-R),[],'all') < 1e-10;
    obj.showResult('getDensity (basis-set orbital, spin restricted)' ,R && rho.isInit);


    occ                 =   {[1 2 pi/2 1/3 .5],[exp(1) log(5) 1]};
    d_SP_b.DFT_allocateOrbital([5 3],proj,randStr,true);
    d_SP_b.DFT_reconstructOrbital(proj,KSO);
    
    Rup                 =   occ{1}(1)*KSO.KSOup(:,1).^2 + occ{1}(2)*KSO.KSOup(:,2).^2 + ...
                            occ{1}(3)*KSO.KSOup(:,3).^2 + occ{1}(4)*KSO.KSOup(:,4).^2 + ...
                            occ{1}(5)*KSO.KSOup(:,5).^2;
    Rdw                 =   occ{2}(1)*KSO.KSOdw(:,1).^2 + occ{2}(2)*KSO.KSOdw(:,2).^2 + ...
                            occ{2}(3)*KSO.KSOdw(:,3).^2;
    rho.getDensity(occ,proj);
    R                   =   max(abs(rho.densityUp  -Rup),[],'all') < 1e-10  && ...
                            max(abs(rho.densityDown-Rdw),[],'all') < 1e-10;
    obj.showResult('getDensity (basis-set orbital, spin polarized)' ,R && rho.isInit);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

