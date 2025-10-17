function [E,DE] = getEnergy(obj)
%getEnergy compute the member wave function energy and error
%
%   [E,err] = obj.getEnergy computes the energies E and errors err for the
%     member-property waveFunction.

    % Compute energy components
    if nargout == 1
        % Only care about the energies
        E               =   NaN(size(obj.wfcn,2),1);
        for k = 1:size(obj.wfcn,2)
            E(k)        =   real(sum(conj(obj.wfcn(:,k)).*(obj.CI*obj.wfcn(:,k))));
        end
    else
        % Want both the energy and error
        E               =   NaN(size(obj.wfcn,2),1);
        DE              =   NaN(size(obj.wfcn,2),1);
        for k = 1:size(obj.wfcn,2)
            Hp          =   obj.CI*obj.wfcn(:,k);
            E(k)        =   real(sum(conj(obj.wfcn(:,k)).*Hp));
            DE(k)       =   sqrt(sum(abs(Hp-E(k)*obj.wfcn(:,k)).^2));
        end
    end
end

