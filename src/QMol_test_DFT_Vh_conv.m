classdef QMol_test_DFT_Vh_conv < QMol_test
%QMol_test_DFT_Vh_conv suite of unit tests for QMol_DFT_Vh_conv

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
    fprintf('  * QMol_test_DFT_Vh_conv\n'); 
    QMol_test_DFT_Vh_conv.version;
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

    s_ee                =   pi;
    Vee                 =   @(x)    exp(-x.^2 * .5/s_ee^2);
    
    Vext                =   QMol_DFT_Vext;
    Vh                  =   QMol_DFT_Vh_conv('interactionPotential',Vee);
    DFT                 =   QMol_DFT_spinRes(...
                                'xspan',                -15:.1:20,          ...
                                'occupation',           [5 5 1.3],          ...
                                'externalPotential',    Vext,               ...
                                'HartreePotential',     Vh);
    DFT.initialize;
    rho                 =   DFT.disc.DFT_allocateDensity;

    % No SIC
    DFT.set('SIC','none');  DFT.initialize;

    x0                  =   5;
    s                   =   2;
    rho.set('rho',exp(-(DFT.x(:)-x0).^2 * .5/s^2));
     
    a                   =   .5/s_ee^2 + .5/s^2;
    b                   =   (DFT.x(:)-x0) * .5/s^2;         D_b     =   .5/s^2;
    c                   =   (DFT.x(:)-x0).^2 * .5/s^2;      D_c     =   (DFT.x(:)-x0)/s^2;
    V_th                =   exp(b.^2/a-c) * sqrt(pi/a);

    V                   =   Vh.getPotential(rho);
    R                   =   isa(V,'QMol_DFT_Vks')                           && ... potential object
                            max(abs(V.potential - V_th)) < 1e-10;                % potential discretization
    obj.showResult('getPotential (no SIC and new)' ,R);

    DV_th               =   (2*D_b.*b/a - D_c) .* V_th;
    DV                  =   Vh.getPotentialDerivative(1,rho);
    R                   =   isa(DV,'QMol_DFT_Vks_grad')                     && ... potential object
                            max(abs(DV.potentialGradient - DV_th)) < 1e-10;      % potential discretization
    obj.showResult('getPotentialDerivative (no SIC and new)' ,R);
    
    % ADSIC
    DFT.set('SIC','averageDensity'); DFT.initialize;
    rho.rho             =   sqrt(rho.rho / sum(DFT.occ));
    for k = 1:numel(DFT.occ),   DFT.KSO.KSO(:,k)=   rho.rho;    end              % use DFT orbitals

    N                   =   sum(DFT.occupation);
    Vh.getPotential([],V,true);
    R                   =   isa(V,'QMol_DFT_Vks')                           && ... potential object
                            max(abs(V.potential - V_th*(2*N-1)/N)) < 1e-10;      % potential discretization
    obj.showResult('getPotential (ADSIC and add to existing potential)' ,R);

    Vh.getPotentialDerivative(1,[],DV,true);
    R                   =   isa(DV,'QMol_DFT_Vks_grad')                             && ... potential object
                            max(abs(DV.potentialGradient - DV_th*(2*N-1)/N)) < 1e-10;    % potential discretization
    obj.showResult('getPotentialDerivative (ADSIC and add to existing potential)' ,R);

    % Energy
    DFT.getDensity(rho);
    DFT.set('SIC','none'); DFT.initialize;
    E_th                =   0.5*sum(rho.rho.*V_th)*(DFT.x(2)-DFT.x(1));
    R                   =   abs(Vh.getEnergy(rho)-E_th) < 1e-10;
    obj.showResult('getEnergy (no SIC)' ,R);

    DFT.set('SIC','averageDensity'); DFT.initialize;
    N                   =   sum(DFT.occupation);
    R                   =   abs(Vh.getEnergy-E_th*(N-1)/N) < 1e-10;         % use DFT orbitals
    obj.showResult('getEnergy (ADSIC)' ,R);

end
function runSpinPolarized(obj) %===========================================
%runSpinRestricted unit tests for spin-restricted DFT
    
    % Initialization
    obj.showSection('Spin-polarized DFT');

    s_ee                =   pi;
    Vee                 =   @(x)    exp(-x.^2 * .5/s_ee^2);
    
    Vext                =   QMol_DFT_Vext;
    Vh                  =   QMol_DFT_Vh_conv('interactionPotential',Vee);
    DFT                 =   QMol_DFT_spinPol(...
                                'xspan',                -15:.1:20,          ...
                                'occupation',           {[5 2.4],[2 2 .1]}, ...
                                'externalPotential',    Vext,               ...
                                'HartreePotential',     Vh);
    DFT.initialize;
    rho                 =   DFT.disc.DFT_allocateDensity;

    % No SIC
    DFT.set('SIC','none');  DFT.initialize;
    
    x0_u                =   5;
    s_u                 =   2;
    rho.set('rhoUp',exp(-(DFT.x(:)-x0_u).^2 * .5/s_u^2));
    
    x0_d                =   -3;
    s_d                 =   .9;
    rho.set('rhoDw',exp(-(DFT.x(:)-x0_d).^2 * .5/s_d^2));
    
    a_u                 =   .5/s_ee^2 + .5/s_u^2;
    b_u                 =   (DFT.x(:)-x0_u) * .5/s_u^2;     D_b_u   =   .5/s_u^2;
    c_u                 =   (DFT.x(:)-x0_u).^2 * .5/s_u^2;  D_c_u   =   (DFT.x(:)-x0_u)/s_u^2;
    a_d                 =   .5/s_ee^2 + .5/s_d^2;
    b_d                 =   (DFT.x(:)-x0_d) * .5/s_d^2;     D_b_d   =   .5/s_d^2;
    c_d                 =   (DFT.x(:)-x0_d).^2 * .5/s_d^2;  D_c_d   =   (DFT.x(:)-x0_d)/s_d^2;
    V_th                =   exp(b_u.^2/a_u-c_u) * sqrt(pi/a_u) + exp(b_d.^2/a_d-c_d) * sqrt(pi/a_d);
    
    V                   =   Vh.getPotential(rho);
    R                   =   isa(V,'QMol_DFT_Vks')                           && ... potential object
                            max(abs(V.potentialUp   - V_th)) < 1e-10        && ... potential values
                            max(abs(V.potentialDown - V_th)) < 1e-10;
    obj.showResult('getPotential (no SIC and new)' ,R);

    DV_th               =   (2*D_b_u.*b_u/a_u - D_c_u) .* exp(b_u.^2/a_u-c_u) * sqrt(pi/a_u) + ...
                            (2*D_b_d.*b_d/a_d - D_c_d) .* exp(b_d.^2/a_d-c_d) * sqrt(pi/a_d);
                            
    DV                  =   Vh.getPotentialDerivative(1,rho);
    R                   =   isa(DV,'QMol_DFT_Vks_grad')                         && ... potential object
                            max(abs(DV.potentialGradientUp   - DV_th)) < 1e-10  && ... potential discretization
                            max(abs(DV.potentialGradientDown - DV_th)) < 1e-10;
    obj.showResult('getPotentialDerivative (no SIC and new)' ,R);

    % ADSIC
    DFT.set('SIC','averageDensity'); DFT.initialize;
    V.add(0,1);

    rho.rhoUp           =   sqrt(rho.rhoUp / sum(DFT.occ{1}));
    for k = 1:numel(DFT.occ{1}),    DFT.KSO.KSOup(:,k)  =   rho.rhoUp;      end      % use DFT orbitals
    rho.rhoDw           =   sqrt(rho.rhoDw / sum(DFT.occ{2}));
    for k = 1:numel(DFT.occ{2}),    DFT.KSO.KSOdw(:,k)  =   rho.rhoDw;      end      % use DFT orbitals

    N                   =   sum(DFT.occupation{1})+sum(DFT.occupation{2});
    Vh.getPotential([],V,true);
    R                   =   isa(V,'QMol_DFT_Vks')                               && ... potential object
                            max(abs(V.potentialUp  -V_th*(2*N-1)/N  ))<1e-10    && ... potential values
                            max(abs(V.potentialDown-V_th*(2*N-1)/N-1))<1e-10;
    obj.showResult('getPotential (ADSIC and add to existing potential)' ,R);
     
    DV.add(1,1,0)
    DV                  =   Vh.getPotentialDerivative(1,[],DV,true);
    R                   =   isa(DV,'QMol_DFT_Vks_grad')                                 && ... potential object
                            max(abs(DV.potentialGradientUp  -DV_th*(2*N-1)/N-1))<1e-10  && ... potential values
                            max(abs(DV.potentialGradientDown-DV_th*(2*N-1)/N  ))<1e-10;
    obj.showResult('getPotentialDerivative (ADSIC and add to existing potential)' ,R);

    % Energy
    DFT.getDensity(rho);
    DFT.set('SIC','none'); DFT.initialize;
    E_th                =   0.5*sum((rho.rhoUp+rho.rhoDw).*V_th)*(DFT.x(2)-DFT.x(1));
    R                   =   abs(Vh.getEnergy(rho)-E_th) < 1e-10;
    obj.showResult('getEnergy (no SIC)' ,R);

    DFT.set('SIC','averageDensity'); DFT.initialize;
    N                   =   sum(DFT.occupation{1})+sum(DFT.occupation{2});
    R                   =   abs(Vh.getEnergy-E_th*(N-1)/N) < 1e-10;         % use DFT orbitals
    obj.showResult('getEnergy (ADSIC)' ,R);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

