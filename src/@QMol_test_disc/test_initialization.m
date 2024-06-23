function test_initialization(obj)
%test_initialization unit tests for the initialization of the class

    % Initialization
    obj.showSection('Object initialization');
    X                   =   -15:.1:20;
    
    disc                =   QMol_disc('xspan',X);
    disc.initialize([]);
    
    % Test gradient and Laplacian
    X                   =   X.';
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
end