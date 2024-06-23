function test_SE_miscellaneous(obj)
%test_SE_miscellaneous unit tests for Schrodinger-equation miscellaneous
%   methods

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('Schrodinger-equation miscellaneous functions');
    randStr             =   RandStream('dsfmt19937','Seed',0);
    
    disc                =   QMol_disc('xspan',-15:.1:20);
    disc.initialize(QMol_DFT_spinRes('discretization',disc));

    % Normalize Kohn-Sham orbitals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    x                   =   disc.x(:);
    p                   =   [exp(-(x-3).^2), x.*exp(-(x+1).^2/.7), cos(2*x).*exp(-x.^2)];
    p                   =   disc.SE_normalizeWaveFunction(p);
    R                   =   max(abs(sum(abs(p).^2,1)*(x(2)-x(1)) - 1)) < 1e-10;
    obj.showResult('SE_normalizeWaveFunction',R);

    % Generate a random orbital (for eigs) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dx                  =   disc.x(2)-disc.x(1);

    p                   =   disc.SE_randomWaveFunction(randStr);
    R                   =   abs(sum(abs(p).^2)*dx - 1) < 1e-10;             % Normalized (no symmetry)
    obj.showResult('SE_randomWaveFunction (no symmetry)',R);

    p                   =   disc.SE_randomWaveFunction(randStr,1);
    R                   =   abs(sum(abs(p).^2)*dx - 1) < 1e-10   &&       ... Normalized 
                            max(abs(p-flip(p))) < 1e-10;                    % Symmetry condition
    obj.showResult('SE_randomWaveFunction (symmetric)',R);

    p                   =   disc.SE_randomWaveFunction(randStr,-1);
    R                   =   abs(sum(abs(p).^2)*dx - 1) < 1e-10   &&       ... Normalized
                            max(abs(p+flip(p))) < 1e-10;                    % Symmetry condition
    obj.showResult('SE_randomWaveFunction (antisymmetric)',R);

    % Kohn-Sham orbital size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showResult('SE_sizeWaveFunction',all(disc.SE_sizeWaveFunction == [numel(disc.x) 1]));


end