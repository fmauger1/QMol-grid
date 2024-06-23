function test_DFT_allocation(obj)
%test_DFT_allocation tests data/object allocation for DFT simulations

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];

    d_SR                =   QMol_disc_basis('xspan',x,'basis',v);
    d_SR.initialize(QMol_DFT_spinRes('discretization',d_SR));
    d_SR.orthonormalizeBasis;

    d_SP                =   QMol_disc_basis('xspan',x,'basis',v);
    d_SP.initialize(QMol_DFT_spinPol('discretization',d_SP));
    d_SP.orthonormalizeBasis;

    % Density allocation (spin restricted) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('DFT one-body-density (rho) allocation (spin restricted)');
    
    rho                 =   d_SR.DFT_allocateDensity([],true);
    R                   =   ~rho.isSpinPol                                  && ...
                            isempty(rho.rhoUp) && isempty(rho.rhoDw)        && ...
                            all(size(rho.rho) == [numel(d_SR.x) 1]);
    obj.showResult('DFT_allocateDensity (new rho)' ,R);
    
    % Reset density
    rho.set('isSpinPol',true,'rho',rand([numel(d_SR.x) 1]));

    d_SR.DFT_allocateDensity(rho);
    R                   =   ~rho.isSpinPol                                  && ... spin restricted
                            isempty(rho.density);                                % no discretization data <=> 0
    obj.showResult('DFT_allocateDensity (reset rho)' ,R);

    % Density allocation (spin polarized) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('DFT one-body-density (rho) allocation (spin polarized)');
    
    rho                 =   d_SP.DFT_allocateDensity([],true);
    R                   =   rho.isSpinPol                                   && ...
                            isempty(rho.rho)                                && ...
                            all(size(rho.rhoUp) == [numel(d_SP.x) 1])       && ...
                            all(size(rho.rhoDw) == [numel(d_SP.x) 1]);
    obj.showResult('DFT_allocateDensity (new rho)' ,R);
    
    % Reset density
    rho.set('isSpinPol',false,'rhoUp',rand([numel(d_SR.x) 1]),'rhoDw',rand([numel(d_SR.x) 1]));

    d_SP.DFT_allocateDensity(rho);
    R                   =   rho.isSpinPol                                   && ... spin restricted
                            isempty(rho.rhoUp)   &&   isempty(rho.rhoDw);        % no discretization data <=> 0
    obj.showResult('DFT_allocateDensity (reset rho)' ,R);

    % Allocate new KSO (restricted) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('DFT Kohn-Sham-orbital (KSO) allocation (spin restricted)');
    KSO_SR              =   d_SR.DFT_allocateOrbital(3);
    R                   =   all(size(KSO_SR.KSO) == [size(d_SR.v,2) 3]);
    obj.showResult('DFT_allocateOrbital (new KSO)' ,R);
    
    % Overwrite KSO on old domain
    d_SR.DFT_allocateOrbital(5,KSO_SR);
    R                   =   all(size(KSO_SR.KSO) == [size(d_SR.v,2) 5]);
    obj.showResult('DFT_allocateOrbital (overwrite KSO on old domain)' ,R);

    % Update KSO, fewer orbitals
    K                   =   rand(size(d_SR.v,2),5);
    KSO_SR.set('KSO',K);

    d_SR.DFT_allocateOrbital(3,KSO_SR);
    R                   =   all(KSO_SR.KSO == K(:,1:3),'all');
    obj.showResult('DFT_allocateOrbital (update KSO, fewer orbitals)' ,R);

    % Update KSO, more orbitals
    d_SR.DFT_allocateOrbital(6,KSO_SR);
    R                   =   all(KSO_SR.KSO(:,1:3) == K(:,1:3),'all')        && ... same KSO
                            all(KSO_SR.KSO(:,4:end) == 0,'all');
    obj.showResult('DFT_allocateOrbital (update KSO, more orbitals)' ,R);

    % Random allocation
    KSO_SR              =   d_SR.DFT_allocateOrbital(3,[],randStr);
    R                   =   all(size(KSO_SR.KSO) == [size(d_SR.v,2) 3])      && ...
                            max(abs(KSO_SR.KSO - d_SR.DFT_normalizeOrbital(KSO_SR.KSO)),[],'all') < 1e-10 && ...
                            ~isreal(KSO_SR.KSO);
    obj.showResult('DFT_allocateOrbital (random, complex, new KSO)' ,R);

    d_SR.DFT_allocateOrbital(5,KSO_SR,randStr,true);
    R                   =   all(size(KSO_SR.KSO) == [size(d_SR.v,2) 5])      && ...
                            max(abs(KSO_SR.KSO - d_SR.DFT_normalizeOrbital(KSO_SR.KSO)),[],'all') < 1e-10 && ...
                            isreal(KSO_SR.KSO);
    obj.showResult('DFT_allocateOrbital (random, real, overwrite KSO)' ,R);
    
    % Allocate new KSO (polarized) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('DFT Kohn-Sham-orbital (KSO) allocation (spin polarized)');
    KSO_SP              =   d_SP.DFT_allocateOrbital([3 2]);
    R                   =   all(size(KSO_SP.KSOup) == [size(d_SP.v,2) 3])    && ...
                            all(size(KSO_SP.KSOdw) == [size(d_SP.v,2) 2]);
    obj.showResult('DFT_allocateOrbital (new KSO)' ,R);
    
    % Overwrite KSO on old domain
    d_SP.DFT_allocateOrbital([5 7],KSO_SP);
    R                   =   all(size(KSO_SP.KSOup) == [size(d_SP.v,2) 5])    && ...
                            all(size(KSO_SP.KSOdw) == [size(d_SP.v,2) 7]);
    obj.showResult('DFT_allocateOrbital (overwrite KSO on old domain)' ,R);

    % Update KSO, fewer orbitals
    K                   =   rand(size(d_SP.v,2),7);
    KSO_SP.set('KSOup',K(:,1:5),'KSOdw',K);

    d_SP.DFT_allocateOrbital([3 5],KSO_SP);
    R                   =   all(KSO_SP.KSOup == K(:,1:3),'all')             && ...
                            all(KSO_SP.KSOdw == K(:,1:5),'all');
    obj.showResult('DFT_allocateOrbital (update KSO, fewer orbitals)' ,R);

    % Update KSO, more orbitals
    d_SP.DFT_allocateOrbital([6 8],KSO_SP);
    R                   =   all(KSO_SP.KSOup(:,1:3) == K(:,1:3),'all')      && ...
                            all(KSO_SP.KSOdw(:,1:5) == K(:,1:5),'all')      && ...
                            all(KSO_SP.KSOup(:,4:end) == 0,'all')           && ...
                            all(KSO_SP.KSOdw(:,6:end) == 0,'all');
    obj.showResult('DFT_allocateOrbital (update KSO, more orbitals)' ,R);

    % Random allocation
    KSO_SP              =   d_SP.DFT_allocateOrbital([3 4],[],randStr);
    R                   =   all(size(KSO_SP.KSOup) == [size(d_SP.v,2) 3])    && ...
                            all(size(KSO_SP.KSOdw) == [size(d_SP.v,2) 4])    && ...
                            max(abs(KSO_SP.KSOup - d_SP.DFT_normalizeOrbital(KSO_SP.KSOup)),[],'all') < 1e-10 && ...
                            max(abs(KSO_SP.KSOdw - d_SP.DFT_normalizeOrbital(KSO_SP.KSOdw)),[],'all') < 1e-10 && ...
                            ~isreal(KSO_SP.KSOup) && ~isreal(KSO_SP.KSOdw);
    obj.showResult('DFT_allocateOrbital (random, complex, new KSO)' ,R);

    d_SP.DFT_allocateOrbital([5 3],KSO_SP,randStr,true);
    R                   =   all(size(KSO_SP.KSOup) == [size(d_SP.v,2) 5])    && ...
                            all(size(KSO_SP.KSOdw) == [size(d_SP.v,2) 3])    && ...
                            max(abs(KSO_SP.KSOup - d_SP.DFT_normalizeOrbital(KSO_SP.KSOup)),[],'all') < 1e-10 && ...
                            max(abs(KSO_SP.KSOdw - d_SP.DFT_normalizeOrbital(KSO_SP.KSOdw)),[],'all') < 1e-10 && ...
                            isreal(KSO_SP.KSOup) && isreal(KSO_SP.KSOdw);
    obj.showResult('DFT_allocateOrbital (random, real, overwrite KSO)' ,R);

    % Potential allocation (spin restricted) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('DFT potential (V) allocation (spin restricted)');

    V                   =   d_SR.DFT_allocatePotential;
    R                   =   isempty(V.V) && isempty(V.Vup)                  && ... potential starts off empty
                            isempty(V.Vdw) && isempty(V.Vimp)               && ... 
                            isa(V,'QMol_DFT_Vks_basis') && isempty(V.mV)    && ...
                            isempty(V.mVup) &&  isempty(V.mVdw)             && ... 
                            ~V.isSpinPol;                                        % polarization status
    obj.showResult('DFT_allocatePotential (new V)' ,R);

    % Reset potential
    V.set('isSpinPol',true,'potential',rand([numel(d_SR.x) 1]));

    d_SR.DFT_allocatePotential(V);
    R                   =   ~V.isSpinPol                                    && ... spin restricted
                            isempty(V.potential);                                % no discretization data <=> 0
    obj.showResult('DFT_allocatePotential (reset V)' ,R);

    % Potential allocation (spin polarized) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('DFT potential (V) allocation (spin polarized)');
    
    V                   =   d_SP.DFT_allocatePotential;
    R                   =   isempty(V.V) && isempty(V.Vup)                  && ... potential starts off empty
                            isempty(V.Vdw) && isempty(V.Vimp)               && ...  
                            isa(V,'QMol_DFT_Vks_basis') && isempty(V.mV)    && ...
                            isempty(V.mVup) &&  isempty(V.mVdw)             && ...
                            V.isSpinPol;
    obj.showResult('DFT_allocatePotential (new V)' ,R);
    
    % Reset potential
    V.set('isSpinPol',false,'potentialUp',rand([numel(d_SR.x) 1]),'potentialDown',rand([numel(d_SR.x) 1]));

    d_SP.DFT_allocatePotential(V);
    R                   =   V.isSpinPol                                     && ... spin restricted
                            isempty(V.Vup)   &&   isempty(V.Vdw);                % no discretization data <=> 0
    obj.showResult('DFT_allocatePotential (reset V)' ,R);

    % Potential gradient allocation (spin restricted) ~~~~~~~~~~~~~~~~~~~~~
    fprintf('    > DFT_allocatePotentialGradient is not (yet) implemented          ****\n');

end