function E = SE_energyKinetic(obj,wfcn)
%SE_kineticEnergy computes the kinetic energy of the member SE, and
%   associated wave function, object

    % Parse orbitals
    E                   =   0;
    IR                  =   isreal(wfcn.wfcn);

    if isempty(obj.Tv), for k = 1:size(wfcn.wfcn,2)                         %#ok<ALIGN> 
        if IR,  E   =   E +      sum(     wfcn.wfcn(:,k) .*real(ifft(obj.T .*fft(wfcn.wfcn(:,k)))));
        else,   E   =   E + real(sum(conj(wfcn.wfcn(:,k)).*     ifft(obj.T .*fft(wfcn.wfcn(:,k))) ));     end
    end, else, for k = 1:size(wfcn.wfcn,2)                                  %#ok<ALIGN> 
        if IR,  E   =   E + real(sum(     wfcn.wfcn(:,k) .*     ifft(obj.Tv.*fft(wfcn.wfcn(:,k)))));
        else,   E   =   E + real(sum(conj(wfcn.wfcn(:,k)).*     ifft(obj.Tv.*fft(wfcn.wfcn(:,k))) ));     end
    end, end

    % Finalize
    E               =   E*obj.dx;

end