function test_SE_allocation(obj)
%test_SE_allocation tests data/object allocation for Schrodinger-equation 
%   simulations

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    disc                =   QMol_disc('xspan',-15:.1:20);       
    disc.initialize(QMol_SE('discretization',disc));

    % Allocate new wave function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('Schrodinger-equation (SE) wave function (wfcn) allocation');
    wfcn                =   disc.SE_allocateWaveFunction(3);
    R                   =   all(size(wfcn.wfcn) == [numel(disc.x) 3]);
    obj.showResult('SE_allocateWaveFunction (new wfcn)' ,R);
    
    % Overwrite KSO on old domain
    disc.SE_allocateWaveFunction(5,wfcn);
    R                   =   all(size(wfcn.wfcn) == [numel(disc.x) 5]);
    obj.showResult('SE_allocateWaveFunction (overwrite wfcn on old domain)' ,R);

    % Update KSO, fewer orbitals
    K                   =   rand(numel(disc.x),5);
    wfcn.set('wfcn',K);

    disc.SE_allocateWaveFunction(3,wfcn);
    R                   =   all(wfcn.wfcn == K(:,1:3),'all');
    obj.showResult('SE_allocateWaveFunction (update wfcn, fewer wfcn)' ,R);

    % Update KSO, more orbitals
    disc.SE_allocateWaveFunction(6,wfcn);
    R                   =   all(wfcn.wfcn(:,1:3) == K(:,1:3),'all')        && ... same KSO
                            all(wfcn.wfcn(:,4:end) == 0,'all');
    obj.showResult('SE_allocateWaveFunction (update wfcn, more wfcn)' ,R);

    % Random allocation
    randStr             =   RandStream('dsfmt19937','Seed',0);

    wfcn                =   disc.SE_allocateWaveFunction(3,[],randStr);
    R                   =   all(size(wfcn.wfcn) == [numel(disc.x) 3])      && ...
                            max(abs(wfcn.wfcn - disc.SE_normalizeWaveFunction(wfcn.wfcn)),[],'all') < 1e-10 && ...
                            ~isreal(wfcn.wfcn);
    obj.showResult('SE_allocateWaveFunction (random, complex, new wfcn)' ,R);

    disc.SE_allocateWaveFunction(5,wfcn,randStr,true);
    R                   =   all(size(wfcn.wfcn) == [numel(disc.x) 5])      && ...
                            max(abs(wfcn.wfcn - disc.SE_normalizeWaveFunction(wfcn.wfcn)),[],'all') < 1e-10 && ...
                            isreal(wfcn.wfcn);
    obj.showResult('SE_allocateWaveFunction (random, real, overwrite wfcn)' ,R);
    
end