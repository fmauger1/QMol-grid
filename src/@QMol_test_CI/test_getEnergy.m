function test_getEnergy(obj)

    % Initialization
    obj.showSection('Wave function energy (getEnergy)');

    randStr             =   RandStream('dsfmt19937','Seed',0);

    % Real CI matrix, real wave function
    M                   =   rand(randStr,[10 10])-.5;
    M                   =   M + M';

    W                   =   rand(randStr,[10 5])-.5;
    W                   =   W./sqrt(sum(abs(W).^2,1));

    E_th                =   sum(W .* (M*W),1);
    DE_th               =   sqrt(sum(abs(M*W - E_th.*W).^2,1));

    CI                  =   QMol_CI_conv('CImatrix',M,'waveFunction',W);
    [E,DE]              =   CI.getEnergy;
    R_E                 =   max( E_th(:)- E(:)) < 1e-10;
    R_DE                =   max(DE_th(:)-DE(:)) < 1e-10;
    obj.showResult('real CI matrix, real wave function',R_E && R_DE);

    % Real CI matrix, complex wave function
    W                   =   rand(randStr,[10 3],like=1i)-.5-.5i;
    W                   =   W./sqrt(sum(abs(W).^2,1));

    E_th                =   sum(conj(W) .* (M*W),1);
    DE_th               =   sqrt(sum(abs(M*W - E_th.*W).^2,1));

    CI.set('waveFunction',W);
    [E,DE]              =   CI.getEnergy;
    R_E                 =   max( E_th(:)- E(:)) < 1e-10;
    R_DE                =   max(DE_th(:)-DE(:)) < 1e-10;
    obj.showResult('real CI matrix, complex wave function',R_E && R_DE);

    % Complex CI matrix, real wave function
    M                   =   rand(randStr,[15 15],like=1i)-.5-.5i;
    M                   =   M + M';

    W                   =   rand(randStr,[15 6])-.5;
    W                   =   W./sqrt(sum(abs(W).^2,1));

    E_th                =   sum(W .* (M*W),1);
    DE_th               =   sqrt(sum(abs(M*W - E_th.*W).^2,1));

    CI.set('CImatrix',M,'waveFunction',W);
    [E,DE]              =   CI.getEnergy;
    R_E                 =   max( E_th(:)- E(:)) < 1e-10;
    R_DE                =   max(DE_th(:)-DE(:)) < 1e-10;
    obj.showResult('complex CI matrix, real wave function',R_E && R_DE);

    % Complex CI matrix, complex wave function
    W                   =   rand(randStr,[15 4],like=1i)-.5-.5i;
    W                   =   W./sqrt(sum(abs(W).^2,1));

    E_th                =   sum(conj(W) .* (M*W),1);
    DE_th               =   sqrt(sum(abs(M*W - E_th.*W).^2,1));

    CI.set('waveFunction',W);
    [E,DE]              =   CI.getEnergy;
    R_E                 =   max( E_th(:)- E(:)) < 1e-10;
    R_DE                =   max(DE_th(:)-DE(:)) < 1e-10;
    obj.showResult('complex CI matrix, complex wave function',R_E && R_DE);

end

