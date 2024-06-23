function E = DFT_energyKinetic(obj,occ,KSO)
%DFT_kineticEnergy computes the kinetic energy of the member DFT, and
%   associated orbitals, object

    % Spin polarized -----------------------
    if KSO.isSpinPol
        % Parse spin-up orbitals
        Eup             =   0;
        IR              =   isreal(KSO.KSOup);

        if isempty(obj.Tv), for k = 1:numel(occ{1})                         %#ok<ALIGN> 
            if IR,  Eup =   Eup + occ{1}(k)*     sum(     KSO.KSOup(:,k) .*real(ifft(obj.T .*fft(KSO.KSOup(:,k)))));
            else,   Eup =   Eup + occ{1}(k)*real(sum(conj(KSO.KSOup(:,k)).*     ifft(obj.T .*fft(KSO.KSOup(:,k))) ));end
        end, else, for k = 1:numel(occ{1})                                  %#ok<ALIGN> 
            if IR,  Eup =   Eup + occ{1}(k)*real(sum(     KSO.KSOup(:,k) .*     ifft(obj.Tv.*fft(KSO.KSOup(:,k)))));
            else,   Eup =   Eup + occ{1}(k)*real(sum(conj(KSO.KSOup(:,k)).*     ifft(obj.Tv.*fft(KSO.KSOup(:,k))) ));end
        end, end

        % Parse spin-down orbitals
        Edw             =   0;
        IR              =   isreal(KSO.KSOdw);

        if isempty(obj.Tv), for k = 1:numel(occ{2})                         %#ok<ALIGN> 
            if IR,  Edw =   Edw + occ{2}(k)*     sum(     KSO.KSOdw(:,k) .*real(ifft(obj.T .*fft(KSO.KSOdw(:,k)))));
            else,   Edw =   Edw + occ{2}(k)*real(sum(conj(KSO.KSOdw(:,k)).*     ifft(obj.T .*fft(KSO.KSOdw(:,k))) ));  end
        end, else, for k = 1:numel(occ{2})                                  %#ok<ALIGN> 
            if IR,  Edw =   Edw + occ{2}(k)*real(sum(     KSO.KSOdw(:,k) .*     ifft(obj.Tv.*fft(KSO.KSOdw(:,k)))));
            else,   Edw =   Edw + occ{2}(k)*real(sum(conj(KSO.KSOdw(:,k)).*     ifft(obj.Tv.*fft(KSO.KSOdw(:,k))) ));  end
        end, end

        % Put it all together
        E               =   [Eup Edw] * obj.dx;
    % Spin restricted ----------------------
    else
        % Parse orbitals
        E               =   0;
        IR              =   isreal(KSO.KSO);

        if isempty(obj.Tv), for k = 1:numel(occ)                            %#ok<ALIGN> 
            if IR,  E   =   E + occ(k)*     sum(     KSO.KSO(:,k) .*real(ifft(obj.T .*fft(KSO.KSO(:,k)))));
            else,   E   =   E + occ(k)*real(sum(conj(KSO.KSO(:,k)).*     ifft(obj.T .*fft(KSO.KSO(:,k))) ));     end
        end, else, for k = 1:numel(occ)                                     %#ok<ALIGN> 
            if IR,  E   =   E + occ(k)*real(sum(     KSO.KSO(:,k) .*     ifft(obj.Tv.*fft(KSO.KSO(:,k)))));
            else,   E   =   E + occ(k)*real(sum(conj(KSO.KSO(:,k)).*     ifft(obj.Tv.*fft(KSO.KSO(:,k))) ));     end
        end, end

        % Finalize
        E               =   E*obj.dx;
    end
end