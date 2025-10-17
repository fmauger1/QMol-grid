function test_computeGroundState(obj)

    % Initialization
    obj.showSection('Ground and first few excited states (computeGroundState)');

    randStr             =   RandStream('dsfmt19937','Seed',0);

    % Small real CI matrix
    M                   =   rand(randStr,[20 20])-.5;
    M                   =   M + M';
    CI                  =   QMol_CI_conv('CImatrix',M);

    CI.computeGroundState(3)
    [E,DE]              =   CI.getEnergy;
    E_th                =   sum(CI.waveFunction .* (M*CI.waveFunction),1);

    R                   =   max( E_th(:)- E(:)) < 1e-10   &&   max(abs(DE)) < 1e-10;
    obj.showResult('small real CI matrix',R);

    % Large real CI matrix
    M                   =   rand(randStr,[2000 2000])-.5;
    M                   =   M + M';
    CI.set('CImatrix',M);

    CI.computeGroundState(3)
    [E,DE]              =   CI.getEnergy;
    E_th                =   sum(CI.waveFunction .* (M*CI.waveFunction),1);

    R                   =   max( E_th(:)- E(:)) < 1e-10   &&   max(abs(DE)) < 1e-10;
    obj.showResult('large real CI matrix',R);

    % Small complex CI matrix
    M                   =   rand(randStr,[15 15],like=1i)-.5-.5i;
    M                   =   M + M';
    CI                  =   QMol_CI_conv('CImatrix',M);

    CI.computeGroundState(3)
    [E,DE]              =   CI.getEnergy;
    E_th                =   sum(conj(CI.waveFunction) .* (M*CI.waveFunction),1);

    R                   =   max( E_th(:)- E(:)) < 1e-10   &&   max(abs(DE)) < 1e-10;
    obj.showResult('small complex CI matrix',R);

    % Large complex CI matrix
    M                   =   rand(randStr,[1500 1500],like=1i)-.5-.5i;
    M                   =   M + M';
    CI                  =   QMol_CI_conv('CImatrix',M);

    CI.computeGroundState(3)
    [E,DE]              =   CI.getEnergy;
    E_th                =   sum(conj(CI.waveFunction) .* (M*CI.waveFunction),1);

    R                   =   max( E_th(:)- E(:)) < 1e-10   &&   max(abs(DE)) < 1e-10;
    obj.showResult('large complex CI matrix',R);

end

