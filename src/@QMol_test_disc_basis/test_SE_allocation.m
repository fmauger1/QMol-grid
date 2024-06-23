function test_SE_allocation(obj)
%test_SE_allocation tests data/object allocation for Schrodinger-equation
%   simulations

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];

    disc                =   QMol_disc_basis('xspan',x,'basis',v);
    disc.initialize(QMol_SE('discretization',disc));
    disc.orthonormalizeBasis;

    
    % Allocate new wave function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('Schrodinger-equation wave-function allocation');
    wfcn                =   disc.SE_allocateWaveFunction(3);
    R                   =   all(size(wfcn.wfcn) == [size(disc.v,2) 3]);
    obj.showResult('SE_allocateWaveFunction (new wfcn)' ,R);
    
    % Overwrite wave functions on old domain
    disc.SE_allocateWaveFunction(5,wfcn);
    R                   =   all(size(wfcn.wfcn) == [size(disc.v,2) 5]);
    obj.showResult('SE_allocateWaveFunction (overwrite wfcn on old domain)' ,R);

    % Update wave function, fewer orbitals
    K                   =   rand(size(disc.v,2),5);
    wfcn.set('wfcn',K);

    disc.SE_allocateWaveFunction(3,wfcn);
    R                   =   all(wfcn.wfcn == K(:,1:3),'all');
    obj.showResult('SE_allocateWaveFunction (update wfcn, fewer wfcn)' ,R);

    % Update wfcn, more orbitals
    disc.SE_allocateWaveFunction(6,wfcn);
    R                   =   all(wfcn.wfcn(:,1:3) == K(:,1:3),'all')         && ... same KSO
                            all(wfcn.wfcn(:,4:end) == 0,'all');
    obj.showResult('SE_allocateWaveFunction (update wfcn, more wfcn)' ,R);

    % Random allocation
    wfcn                =   disc.SE_allocateWaveFunction(3,[],randStr);
    R                   =   all(size(wfcn.wfcn) == [size(disc.v,2) 3])      && ...
                            max(abs(wfcn.wfcn - disc.SE_normalizeWaveFunction(wfcn.wfcn)),[],'all') < 1e-10 && ...
                            ~isreal(wfcn.wfcn);
    obj.showResult('SE_allocateWaveFunction (random, complex, new wfcn)' ,R);

    disc.SE_allocateWaveFunction(5,wfcn,randStr,true);
    R                   =   all(size(wfcn.wfcn) == [size(disc.v,2) 5])      && ...
                            max(abs(wfcn.wfcn - disc.SE_normalizeWaveFunction(wfcn.wfcn)),[],'all') < 1e-10 && ...
                            isreal(wfcn.wfcn);
    obj.showResult('SE_allocateWaveFunction (random, real, overwrite wfcn)' ,R);
    
end