function test_SE_waveFunction(obj)
%test_SE_waveFunction units tests for Kohn-Sham orbitals facilities

    % Initialization ======================================================
    obj.showSection('Wave functions');

    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)] .* ...
                            [exp( 2i*x),    exp( 3i*x),       exp(-4i*x),exp(-2i*x)];
    
    disc                =   QMol_disc_basis('xspan',x,'basis',v);    
    disc.initialize(QMol_SE('discretization',disc));
    disc.orthonormalizeBasis;

    % SE_projectWaveFunc ==================================================
    dx                  =   x(2)-x(1);

    R                   =   rand(randStr,[numel(disc.x) 3],'like',1i)-.5-.5i;
    wfcn                =   QMol_SE_wfcn('waveFunction',R);
    proj                =   disc.SE_projectWaveFunction(wfcn);

    P                   =   proj.wfcn;
    for k = 1:size(P,1),   for l = 1:size(P,2) %#ok<ALIGN> 
        P(k,l)          =   P(k,l) - sum(R(:,l).*conj(disc.v(:,k)))*dx;
    end, end
    R                   =   isa(proj,'QMol_SE_wfcn_basis')                  && ... correct object type
                            all(abs(P) < 1e-10,'all');                           % accurate projection
    obj.showResult('SE_projectWaveFunction',R);

    % DFT_reconstructOrbital ==============================================
    proj                =   disc.SE_allocateWaveFunction(3);
    proj.set('wfcn',rand(randStr,[size(v,2),3],'like',1i)-.5-.5i);
    wfcn                 =   disc.SE_reconstructWaveFunction(proj);

    P                   =   wfcn.wfcn;
    for k = 1:disc.basisSize, for l = 1:3                                   %#ok<ALIGN> 
        P(:,l)          =   P(:,l) - proj.wfcn(k,l) * disc.v(:,k);
    end, end
    R                   =  ~isa(wfcn,'QMol_SE_wfcn_basis')                  && ... correct object type (despite overloading)
                            all(abs(P) < 1e-10,'all');                           % accurate projection
    obj.showResult('SE_reconstructWaveFunction',R);

end