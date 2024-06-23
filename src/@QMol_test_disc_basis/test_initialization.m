function test_initialization(obj)
%test_initialization unit tests for the initialization of the class

    % Initialization ======================================================
    obj.showSection('Object initialization');

    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility
    X                   =   (-15:.1:20)';
    
    disc                =   QMol_disc_basis('xspan',X);

    % Basis set orthonormalization ========================================
    V                   =   [exp(-(X-2).^2),exp(-(X-1).^2/.7),exp(-X.^2),exp(-(X+1).^2/2)];
    I                   =   [exp( 2i*X),exp( 3i*X),exp(-4i*X),exp(-2i*X)];

    disc.set('basis',V);    disc.initialize([]);
    obj.showResult('basisSize' ,disc.basisSize == 4);

    disc.set('basis',V);    disc.orthonormalizeBasis('overlap');
    R                   =   true;
    for k = 1:4, for l = 1:k                                                %#ok<ALIGN> 
        R               =   R   &&   abs(sum(conj(disc.v(:,k)).*disc.v(:,l))*disc.dx - (k == l)) < 1e-10;
    end, end
    obj.showResult('orthonormalizeBasis (overlap, real)' ,R);

    disc.set('basis',V.*I); disc.orthonormalizeBasis('overlap');
    R                   =   true;
    for k = 1:4, for l = 1:k                                                %#ok<ALIGN> 
        R               =   R   &&   abs(sum(conj(disc.v(:,k)).*disc.v(:,l))*disc.dx - (k == l)) < 1e-10;
    end, end
    obj.showResult('orthonormalizeBasis (overlap, complex)' ,R);

    disc.set('basis',V);    disc.orthonormalizeBasis('Gram–Schmidt');
    R                   =   true;
    for k = 1:4, for l = 1:k                                                %#ok<ALIGN> 
        R               =   R   &&   abs(sum(conj(disc.v(:,k)).*disc.v(:,l))*disc.dx - (k == l)) < 1e-10;
    end, end
    obj.showResult('orthonormalizeBasis (Gram-Schmidt, real)' ,R);

    disc.set('basis',V.*I); disc.orthonormalizeBasis('Gram–Schmidt');
    R                   =   true;
    for k = 1:4, for l = 1:k                                                %#ok<ALIGN> 
        R               =   R   &&   abs(sum(conj(disc.v(:,k)).*disc.v(:,l))*disc.dx - (k == l)) < 1e-10;
    end, end
    obj.showResult('orthonormalizeBasis (Gram-Schmidt, complex)' ,R);

    disc.set('basis',V);    disc.orthonormalizeBasis('Gram–Schmidt',[4 2 3 1]);
    R                   =   true;
    for k = 1:4, for l = 1:k                                                %#ok<ALIGN> 
        R               =   R   &&   abs(sum(conj(disc.v(:,k)).*disc.v(:,l))*disc.dx - (k == l)) < 1e-10;
    end, end
    obj.showResult('orthonormalizeBasis (ordered Gram-Schmidt, real)' ,R);

    disc.set('basis',V.*I); disc.orthonormalizeBasis('Gram–Schmidt',[4 2 3 1]);
    R                   =   true;
    for k = 1:4, for l = 1:k                                                %#ok<ALIGN> 
        R               =   R   &&   abs(sum(conj(disc.v(:,k)).*disc.v(:,l))*disc.dx - (k == l)) < 1e-10;
    end, end
    obj.showResult('orthonormalizeBasis (ordered Gram-Schmidt, complex)' ,R);
    
    % Test gradient and Laplacian =========================================
    disc.set('basis',[]);   disc.initialize([]);

    x_0                 =   2;
    s                   =   pi/2;
    
    F                   =   X                      .* exp(-(X-x_0).^2*.5/s^2);
    DF                  =   (1+x_0/s^2*X-X.^2/s^2) .* exp(-(X-x_0).^2*.5/s^2);
    TF                  =   -.5*( x_0/s^2-2*X/s^2 - (1+x_0/s^2*X-X.^2/s^2).*(X-x_0)/s^2 ) .* exp(-(X-x_0).^2*.5/s^2);
    
    obj.showResult('D (gradient)' ,max(abs( real(ifft(disc.D.*fft(F)))-DF )) < 1e-10);
    obj.showResult('T (Laplacian)',max(abs( real(ifft(disc.T.*fft(F)))-TF )) < 1e-10);

    A                   =   .1;
    disc.setTv(A);  TvF =   TF - A*1i*DF + .5*A^2*F;
    obj.showResult('Tv (velocity gauge Laplacian)',max(abs( ifft(disc.Tv.*fft(F))-TvF )) < 1e-10);

    % Test gradient and Laplacian matrices ================================
    p                   =   (rand(randStr,[4,1])-.5) + 1i*(rand(randStr,[4,1])-.5);
    p                   =   p / sqrt(sum(abs(p).^2));

    disc.set('basis',V);    disc.orthonormalizeBasis;   disc.initialize([]);

    Dp                  =   disc.mD * p;
    Tp                  =   disc.mT * p;
    p                   =   p(1)*disc.v(:,1) + p(2)*disc.v(:,2)+...
                            p(3)*disc.v(:,3) + p(4)*disc.v(:,4);
    for k = 1:4
        Dp(k)           =   Dp(k) - sum(conj(disc.v(:,k)).*ifft(disc.D.*fft(p)))*disc.dx;
        Tp(k)           =   Tp(k) - sum(conj(disc.v(:,k)).*ifft(disc.T.*fft(p)))*disc.dx;
    end
    obj.showResult('mD (gradient matrix, real basis)',all(abs(Dp) < 1e-10,'all'));
    obj.showResult('mT (Laplacian matrix, real basis)',all(abs(Tp) < 1e-10,'all'));

    A                   =   .1;
    disc.setTv(A);  mTv =   disc.mT - A*1i*disc.mD + .5*A^2*eye(disc.nV);
    obj.showResult('mTv (velocity gauge Laplacian matrix, real basis)',max(abs( mTv-disc.mTv ),[],'all') < 1e-10);
    disc.setTv([]);


    p                   =   (rand(randStr,[4,1])-.5) + 1i*(rand(randStr,[4,1])-.5);
    p                   =   p / sqrt(sum(abs(p).^2));

    disc.set('basis',V.*I); disc.orthonormalizeBasis;   disc.initialize([]);

    Dp                  =   disc.mD * p;
    Tp                  =   disc.mT * p;
    p                   =   p(1)*disc.v(:,1) + p(2)*disc.v(:,2)+...
                            p(3)*disc.v(:,3) + p(4)*disc.v(:,4);
    for k = 1:4
        Dp(k)           =   Dp(k) - sum(conj(disc.v(:,k)).*ifft(disc.D.*fft(p)))*disc.dx;
        Tp(k)           =   Tp(k) - sum(conj(disc.v(:,k)).*ifft(disc.T.*fft(p)))*disc.dx;
    end
    obj.showResult('mD (gradient matrix, complex basis)',all(abs(Dp) < 1e-10,'all'));
    obj.showResult('mT (Laplacian matrix, complex basis)',all(abs(Tp) < 1e-10,'all'));

    A                   =   .1;
    disc.setTv(A);  mTv =   disc.mT - A*1i*disc.mD + .5*A^2*eye(disc.nV);
    obj.showResult('mTv (velocity gauge Laplacian matrix, complex basis)',max(abs( mTv-disc.mTv ),[],'all') < 1e-10);
end