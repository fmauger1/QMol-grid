function test_DFT_orbital(obj)
%test_DFT_orbital units tests for Kohn-Sham orbitals facilities

    % Initialization ======================================================
    obj.showSection('Kohn-Sham orbitals');

    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)] .* ...
                            [exp( 2i*x),    exp( 3i*x),       exp(-4i*x),exp(-2i*x)];
    
    d_SR                =   QMol_disc_basis('xspan',x,'basis',v);    
    d_SR.initialize(QMol_DFT_spinRes('discretization',d_SR));
    d_SR.orthonormalizeBasis;

    d_SP                =   QMol_disc_basis('xspan',x,'basis',v);    
    d_SP.initialize(QMol_DFT_spinPol('discretization',d_SP));
    d_SP.orthonormalizeBasis;

    % DFT_projectOrbital ==================================================
    dx                  =   x(2)-x(1);

    R                   =   rand(randStr,[numel(d_SR.x) 3],'like',1i)-.5-.5i;
    KSO                 =   QMol_DFT_orbital('isSpinPol',false,'orbital',R);
    proj                =   d_SR.DFT_projectOrbital(KSO);

    P                   =   proj.KSO;
    for k = 1:size(P,1),   for l = 1:size(P,2) %#ok<ALIGN> 
        P(k,l)          =   P(k,l) - sum(R(:,l).*conj(d_SR.v(:,k)))*dx;
    end, end
    R                   =  ~proj.isSpinPol                                  && ... spin restricted flag
                            isempty(proj.orbitalUp)                         && ... no spin-polarized components
                            isempty(proj.orbitalDown)                       && ...
                            isa(proj,'QMol_DFT_orbital_basis')              && ... correct object type
                            all(abs(P) < 1e-10,'all');                           % accurate projection
    obj.showResult('DFT_projectOrbital (spin restricted)',R);


    Rup                 =   rand(randStr,[numel(d_SP.x) 3],'like',1i)-.5-.5i;
    Rdw                 =   rand(randStr,[numel(d_SP.x) 5],'like',1i)-.5-.5i;
    KSO                 =   QMol_DFT_orbital('isSpinPol',false,'orbitalUp',Rup,'orbitalDown',Rdw);
    proj                =   d_SR.DFT_projectOrbital(KSO);

    Pup                 =   proj.KSOup;
    Pdw                 =   proj.KSOdw;
    for k = 1:size(Pup,1),   for l = 1:size(Pup,2) %#ok<ALIGN> 
        Pup(k,l)        =   Pup(k,l) - sum(Rup(:,l).*conj(d_SR.v(:,k)))*dx;
    end, end
    for k = 1:size(Pdw,1),   for l = 1:size(Pdw,2) %#ok<ALIGN> 
        Pdw(k,l)        =   Pdw(k,l) - sum(Rdw(:,l).*conj(d_SR.v(:,k)))*dx;
    end, end
    R                   =   proj.isSpinPol                                  && ... spin polarized flag
                            isempty(proj.orbital)                           && ... No spin-restricted components
                            isa(proj,'QMol_DFT_orbital_basis')              && ... correct object type
                            all(abs(Pup) < 1e-10,'all')                     && ... accurate projection
                            all(abs(Pdw) < 1e-10,'all');
    obj.showResult('DFT_projectOrbital (spin polarized)',R);

    % DFT_reconstructOrbital ==============================================
    proj                =   d_SR.DFT_allocateOrbital(3);
    proj.set('KSO',rand(randStr,[size(v,2),3],'like',1i)-.5-.5i);
    KSO                 =   d_SR.DFT_reconstructOrbital(proj);

    P                   =   KSO.KSO;
    for k = 1:d_SR.basisSize, for l = 1:3                                   %#ok<ALIGN> 
        P(:,l)          =   P(:,l) - proj.KSO(k,l) * d_SR.v(:,k);
    end, end
    R                   =  ~KSO.isSpinPol                                   && ... spin restricted flag
                            isempty(KSO.orbitalUp)                          && ... No spin-polarized components
                            isempty(KSO.orbitalDown)                        && ...
                           ~isa(KSO,'QMol_DFT_orbital_basis')               && ... correct object type
                            all(abs(P) < 1e-10,'all');                           % accurate projection
    obj.showResult('DFT_reconstructOrbital (spin restricted)',R);


    proj                =   d_SP.DFT_allocateOrbital([3 5]);
    proj.set('KSOup',rand(randStr,[size(v,2),3],'like',1i)-.5-.5i, ...
             'KSOdw',rand(randStr,[size(v,2),5],'like',1i)-.5-.5i);
    KSO                 =   d_SP.DFT_reconstructOrbital(proj);

    Pup                 =   KSO.KSOup;
    Pdw                 =   KSO.KSOdw;
    for k = 1:d_SR.basisSize, for l = 1:3                                   %#ok<ALIGN> 
        Pup(:,l)        =   Pup(:,l) - proj.KSOup(k,l) * d_SP.v(:,k);
    end, end
    for k = 1:d_SR.basisSize, for l = 1:5                                   %#ok<ALIGN> 
        Pdw(:,l)        =   Pdw(:,l) - proj.KSOdw(k,l) * d_SP.v(:,k);
    end, end
    R                   =   KSO.isSpinPol                                   && ... spin polarized flag
                            isempty(KSO.orbital)                            && ... No spin-restricted components
                           ~isa(KSO,'QMol_DFT_orbital_basis')               && ... correct object type
                            all(abs(Pup) < 1e-10,'all')                     && ... accurate projection
                            all(abs(Pdw) < 1e-10,'all');
    obj.showResult('DFT_reconstructOrbital (spin polarized)',R);

end