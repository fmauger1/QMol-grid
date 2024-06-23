function test_DFT_miscellaneous(obj)
%test_DFT_miscellaneous unit tests for miscellaneous DFT methods

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('DFT miscellaneous functions');
    randStr             =   RandStream('dsfmt19937','Seed',0);
    
    d_SR                =   QMol_disc('xspan',-15:.1:20);
    d_SR.initialize(QMol_DFT_spinRes('discretization',d_SR));

    d_SP                =   QMol_disc('xspan',-15:.1:22);       
    d_SP.initialize(QMol_DFT_spinPol('discretization',d_SP));

    % Distance between two objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    R1                  =   rand(randStr,numel(d_SR.x),1);
    R2                  =   rand(randStr,numel(d_SR.x),1);
    rho1                =   QMol_DFT_density('isSpinPol',false,'rho',R1);
    rho2                =   QMol_DFT_density('isSpinPol',false,'rho',R2);
    
    d                   =   sqrt( sum((R1-R2).^2) * (d_SR.x(2)-d_SR.x(1)) );
    R                   =   abs(d - d_SR.DFT_dist(rho1,rho2)) < 1e-10;
    obj.showResult('DFT_dist (density, spin restricted)' ,R);

    clear rho1 rho2

    R1                  =   rand(randStr,numel(d_SP.x),1);
    R2                  =   rand(randStr,numel(d_SP.x),1);
    R3                  =   rand(randStr,numel(d_SP.x),1);
    R4                  =   rand(randStr,numel(d_SP.x),1);
    rho1                =   QMol_DFT_density('isSpinPol',true,'rho',[],'rhoUp',R1,'rhoDw',R2);
    rho2                =   QMol_DFT_density('isSpinPol',true,'rho',[],'rhoUp',R3,'rhoDw',R4);

    d                   =   sqrt( [sum((R1-R3).^2) sum((R2-R4).^2)] * (d_SP.x(2)-d_SP.x(1)) );
    R                   =   all(abs(d - d_SP.DFT_dist(rho1,rho2)) < 1e-10);
    obj.showResult('DFT_dist (density, spin polarized)' ,R);

    clear rho1 rho2

    % Mix density objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    R1                  =   rand(randStr,numel(d_SR.x),1);
    R2                  =   rand(randStr,numel(d_SR.x),1);
    R3                  =   rand(randStr,numel(d_SR.x),1);
    rho                 =   QMol_DFT_density('isSpinPol',false);
    rho1                =   QMol_DFT_density('isSpinPol',false,'rho',R1);
    rho2                =   QMol_DFT_density('isSpinPol',false,'rho',R2);
    rho3                =   QMol_DFT_density('isSpinPol',false,'rho',R3);

    d_SR.DFT_mix(rho,.2,rho1);
    R                   =   ~rho.isSpinPol                                  && ... Spin restricted
                            all(size(rho.rho) == [numel(d_SR.x),1])         && ... Results at the right place
                            isempty(rho.rhoUp) && isempty(rho.rhoDw)        && ... 
                            all(abs(rho.rho - 0.2*R1) < 1e-10);
    obj.showResult('DFT_mix (density, spin restricted)' ,R);

    d_SR.DFT_mix(rho,1,rho,.2,rho1,.3,rho2,.5,rho3);
    R                   =   ~rho.isSpinPol                                  && ... Spin restricted
                            all(size(rho.rho) == [numel(d_SR.x),1])         && ... Results at the right place
                            isempty(rho.rhoUp) && isempty(rho.rhoDw)        && ... 
                            all(abs(rho.rho - 0.4*R1-.3*R2-.5*R3) < 1e-10);
    obj.showResult('DFT_mix (density, in place, spin restricted)' ,R);
    
    clear rho rho1 rho2 rho3

    R1                  =   rand(randStr,numel(d_SP.x),1);
    R2                  =   rand(randStr,numel(d_SP.x),1);
    R3                  =   rand(randStr,numel(d_SP.x),1);
    R4                  =   rand(randStr,numel(d_SP.x),1);
    R5                  =   rand(randStr,numel(d_SP.x),1);
    R6                  =   rand(randStr,numel(d_SP.x),1);
    rho                 =   QMol_DFT_density('isSpinPol',true,'rho',[],'rhoUp',[],'rhoDw',[]);
    rho1                =   QMol_DFT_density('isSpinPol',true,'rho',[],'rhoUp',R1,'rhoDw',R2);
    rho2                =   QMol_DFT_density('isSpinPol',true,'rho',[],'rhoUp',R3,'rhoDw',R4);
    rho3                =   QMol_DFT_density('isSpinPol',true,'rho',[],'rhoUp',R5,'rhoDw',R6);

    d_SP.DFT_mix(rho,[.2,.3],rho1);
    R                   =   rho.isSpinPol                                   && ... Spin restricted
                            isempty(rho.rho)                                && ... Results at the right place
                            all(size(rho.rhoUp) == [numel(d_SP.x),1])       && ... 
                            all(size(rho.rhoDw) == [numel(d_SP.x),1])       && ... 
                            all(abs(rho.rhoUp - 0.2*R1) < 1e-10)            && ... Results are correct
                            all(abs(rho.rhoDw - 0.3*R2) < 1e-10);
    obj.showResult('DFT_mix (density, spin polarized)' ,R);

    d_SP.DFT_mix(rho,[1 2],rho,[.2 .3],rho1,[.3 .1],rho2,[.5 .7],rho3);
    R                   =   rho.isSpinPol                                   && ... Spin restricted
                            isempty(rho.rho)                                && ... Results at the right place
                            all(size(rho.rhoUp) == [numel(d_SP.x),1])       && ... 
                            all(size(rho.rhoDw) == [numel(d_SP.x),1])       && ... 
                            all(abs(rho.rhoUp - 0.4*R1-.3*R3-.5*R5) < 1e-10)&& ... Results are correct
                            all(abs(rho.rhoDw - 0.9*R2-.1*R4-.7*R6) < 1e-10);
    obj.showResult('DFT_mix (density, in place, spin polarized)' ,R);

    clear rho rho1 rho2 rho3

    % Mix potential objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    V1                  =   rand(randStr,numel(d_SR.x),1);
    V2                  =   rand(randStr,numel(d_SR.x),1);
    V3                  =   rand(randStr,numel(d_SR.x),1);
    Vks                 =   QMol_DFT_Vks;
    Vks1                =   QMol_DFT_Vks('V',V1);
    Vks2                =   QMol_DFT_Vks('V',V2);
    Vks3                =   QMol_DFT_Vks('V',V3);

    d_SR.DFT_mix(Vks,.2,Vks1);
    R                   =   all(size(Vks.V) == [numel(d_SR.x),1])           && ... Results at the right place
                            isempty(Vks.Vup) && isempty(Vks.Vdw)            && ... 
                            all(abs(Vks.V - 0.2*V1) < 1e-10);
    obj.showResult('DFT_mix (potential, spin restricted)' ,R);

    d_SR.DFT_mix(Vks,1,Vks,.2,Vks1,.3,Vks2,.5,Vks3);
    R                   =   all(size(Vks.V) == [numel(d_SR.x),1])           && ... Results at the right place
                            isempty(Vks.Vup) && isempty(Vks.Vdw)            && ... 
                            all(abs(Vks.V - 0.4*V1-.3*V2-.5*V3) < 1e-10);
    obj.showResult('DFT_mix (potential, in place, spin restricted)' ,R);

    clear Vks Vks1 Vks2 Vks3

    V1                  =   rand(randStr,numel(d_SP.x),1);
    V2                  =   rand(randStr,numel(d_SP.x),1);
    V3                  =   rand(randStr,numel(d_SP.x),1);
    V4                  =   rand(randStr,numel(d_SP.x),1);
    V5                  =   rand(randStr,numel(d_SP.x),1);
    V6                  =   rand(randStr,numel(d_SP.x),1);
    Vks                 =   QMol_DFT_Vks('V',[],'Vup',[],'Vdw',[]);
    Vks1                =   QMol_DFT_Vks('V',[],'Vup',V1,'Vdw',V2);
    Vks2                =   QMol_DFT_Vks('V',[],'Vup',V3,'Vdw',V4);
    Vks3                =   QMol_DFT_Vks('V',[],'Vup',V5,'Vdw',V6);

    d_SP.DFT_mix(Vks,[.2,.3],Vks1);
    R                   =   isempty(Vks.V)                                  && ... Results at the right place
                            all(size(Vks.Vup) == [numel(d_SP.x),1])         && ... 
                            all(size(Vks.Vdw) == [numel(d_SP.x),1])         && ... 
                            all(abs(Vks.Vup - 0.2*V1) < 1e-10)              && ... Results are correct
                            all(abs(Vks.Vdw - 0.3*V2) < 1e-10);
    obj.showResult('DFT_mix (potential, spin polarized)' ,R);

    d_SP.DFT_mix(Vks,[1 2],Vks,[.2 .3],Vks1,[.3 .1],Vks2,[.5 .7],Vks3);
    R                   =   isempty(Vks.V)                                  && ... Results at the right place
                            all(size(Vks.Vup) == [numel(d_SP.x),1])         && ... 
                            all(size(Vks.Vdw) == [numel(d_SP.x),1])         && ... 
                            all(abs(Vks.Vup - 0.4*V1-.3*V3-.5*V5) < 1e-10)  && ... Results are correct
                            all(abs(Vks.Vdw - 0.9*V2-.1*V4-.7*V6) < 1e-10);
    obj.showResult('DFT_mix (potential, in place, spin polarized)' ,R);

    clear Vks Vks1 Vks2 Vks3

    % Normalize Kohn-Sham orbitals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    x                   =   d_SR.x(:);
    p                   =   [exp(-(x-3).^2), x.*exp(-(x+1).^2/.7), cos(2*x).*exp(-x.^2)];
    p                   =   d_SR.DFT_normalizeOrbital(p);
    R                   =   max(abs(sum(abs(p).^2,1)*(x(2)-x(1)) - 1)) < 1e-10;
    obj.showResult('DFT_normalizeOrbital',R);

    % Generate a random orbital (for eigs) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dx                  =   d_SR.x(2)-d_SR.x(1);

    p                   =   d_SR.DFT_randomOrbital(randStr);
    R                   =   abs(sum(abs(p).^2)*dx - 1) < 1e-10;             % Normalized (no symmetry)
    obj.showResult('DFT_randomOrbital (no symmetry)',R);

    p                   =   d_SR.DFT_randomOrbital(randStr,1);
    R                   =   abs(sum(abs(p).^2)*dx - 1) < 1e-10   &&       ... Normalized 
                            max(abs(p-flip(p))) < 1e-10;                    % Symmetry condition
    obj.showResult('DFT_randomOrbital (symmetric)',R);

    p                   =   d_SR.DFT_randomOrbital(randStr,-1);
    R                   =   abs(sum(abs(p).^2)*dx - 1) < 1e-10   &&       ... Normalized
                            max(abs(p+flip(p))) < 1e-10;                    % Symmetry condition
    obj.showResult('DFT_randomOrbital (antisymmetric)',R);



    % Anderson mixing coefficient (density) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    R1                  =   rand(randStr,numel(d_SR.x),1);
    R2                  =   rand(randStr,numel(d_SR.x),1);
    R3                  =   rand(randStr,numel(d_SR.x),1);
    R4                  =   rand(randStr,numel(d_SR.x),1);
    rho1                =   QMol_DFT_density('isSpinPol',false,'rho',R1);
    rho2                =   QMol_DFT_density('isSpinPol',false,'rho',R2);
    rho3                =   QMol_DFT_density('isSpinPol',false,'rho',R3);
    rho4                =   QMol_DFT_density('isSpinPol',false,'rho',R4);

    b                   =   d_SR.DFT_SCF_AndersonMixCoeff(rho1,rho2,rho3,rho4);
    R                   =   sum( (R1-R3).*(R1-R3-R2+R4) )/sum( (R1-R3-R2+R4).^2 );
    obj.showResult('DFT_SCF_AndersonMixCoeff (density, spin restricted)' ,abs(b-R) < 1e-10);

    clear rho1 rho2 rho3 rho4

    R1                  =   rand(randStr,numel(d_SP.x),1);
    R2                  =   rand(randStr,numel(d_SP.x),1);
    R3                  =   rand(randStr,numel(d_SP.x),1);
    R4                  =   rand(randStr,numel(d_SP.x),1);
    R5                  =   rand(randStr,numel(d_SP.x),1);
    R6                  =   rand(randStr,numel(d_SP.x),1);
    R7                  =   rand(randStr,numel(d_SP.x),1);
    R8                  =   rand(randStr,numel(d_SP.x),1);
    rho1                =   QMol_DFT_density('isSpinPol',true,'rho',[],'rhoUp',R1,'rhoDw',R2);
    rho2                =   QMol_DFT_density('isSpinPol',true,'rho',[],'rhoUp',R3,'rhoDw',R4);
    rho3                =   QMol_DFT_density('isSpinPol',true,'rho',[],'rhoUp',R5,'rhoDw',R6);
    rho4                =   QMol_DFT_density('isSpinPol',true,'rho',[],'rhoUp',R7,'rhoDw',R8);

    b                   =   d_SP.DFT_SCF_AndersonMixCoeff(rho1,rho2,rho3,rho4);
    R                   =  [sum( (R1-R5).*(R1-R5-R3+R7) )/sum( (R1-R5-R3+R7).^2 ), ...
                            sum( (R2-R6).*(R2-R6-R4+R8) )/sum( (R2-R6-R4+R8).^2 )];
    obj.showResult('DFT_SCF_AndersonMixCoeff (density, spin polarized)' ,all(abs(b-R) < 1e-10));

    clear rho1 rho2 rho3 rho4

    % Anderson mixing coefficient (potential) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    V1                  =   rand(randStr,numel(d_SR.x),1);
    V2                  =   rand(randStr,numel(d_SR.x),1);
    V3                  =   rand(randStr,numel(d_SR.x),1);
    V4                  =   rand(randStr,numel(d_SR.x),1);
    Vks1                =   QMol_DFT_Vks('V',V1);
    Vks2                =   QMol_DFT_Vks('V',V2);
    Vks3                =   QMol_DFT_Vks('V',V3);
    Vks4                =   QMol_DFT_Vks('V',V4);

    b                   =   d_SR.DFT_SCF_AndersonMixCoeff(Vks1,Vks2,Vks3,Vks4);
    V                   =   sum( (V1-V3).*(V1-V3-V2+V4) )/sum( (V1-V3-V2+V4).^2 );
    obj.showResult('DFT_SCF_AndersonMixCoeff (potential, spin restricted)' ,abs(b-V) < 1e-10);

    clear Vks1 Vks2 Vks3 Vks4

    V1                  =   rand(randStr,numel(d_SP.x),1);
    V2                  =   rand(randStr,numel(d_SP.x),1);
    V3                  =   rand(randStr,numel(d_SP.x),1);
    V4                  =   rand(randStr,numel(d_SP.x),1);
    V5                  =   rand(randStr,numel(d_SP.x),1);
    V6                  =   rand(randStr,numel(d_SP.x),1);
    V7                  =   rand(randStr,numel(d_SP.x),1);
    V8                  =   rand(randStr,numel(d_SP.x),1);
    Vks1                =   QMol_DFT_Vks('V',[],'Vup',V1,'Vdw',V2);
    Vks2                =   QMol_DFT_Vks('V',[],'Vup',V3,'Vdw',V4);
    Vks3                =   QMol_DFT_Vks('V',[],'Vup',V5,'Vdw',V6);
    Vks4                =   QMol_DFT_Vks('V',[],'Vup',V7,'Vdw',V8);

    b                   =   d_SP.DFT_SCF_AndersonMixCoeff(Vks1,Vks2,Vks3,Vks4);
    V                   =  [sum( (V1-V5).*(V1-V5-V3+V7) )/sum( (V1-V5-V3+V7).^2 ), ...
                            sum( (V2-V6).*(V2-V6-V4+V8) )/sum( (V2-V6-V4+V8).^2 )];
    obj.showResult('DFT_SCF_AndersonMixCoeff (potential, spin polarized)' ,all(abs(b-V) < 1e-10));

    clear Vks1 Vks2 Vks3 Vks4

    % Kohn-Sham orbital size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showResult('DFT_sizeOrbital',all(d_SR.DFT_sizeOrbital == [numel(d_SR.x) 1]));


end