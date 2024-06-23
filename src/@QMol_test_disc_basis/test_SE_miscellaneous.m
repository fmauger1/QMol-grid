function test_SE_miscellaneous(obj)
%test_SE_miscellaneous unit tests for miscellaneous SE methods

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('Schrodinger-equation miscellaneous functions');

    randStr             =   RandStream('dsfmt19937','Seed',0);
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];
    
    disc                =   QMol_disc_basis('xspan',x,'basis',v);
    disc.initialize(QMol_SE('discretization',disc));
    disc.orthonormalizeBasis;

    % Normalize wave functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    p                   =   rand(randStr,[size(v,2) 6],'like',1i);
    p                   =   disc.SE_normalizeWaveFunction(p);
    R                   =   max(abs(sum(abs(p).^2,1) - 1)) < 1e-10;
    obj.showResult('SE_normalizeWaveFunction',R);

    % Generate a random wave functions (for eigs) ~~~~~~~~~~~~~~~~~~~~~~~~~
    p                   =   disc.SE_randomWaveFunction(randStr);
    R                   =   abs(sum(abs(p).^2) - 1) < 1e-10;
    obj.showResult('SE_randomWaveFunction',R);

    % Wave function size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showResult('SE_sizeWaveFunction',all(disc.SE_sizeWaveFunction == [size(v,2) 1]));


end