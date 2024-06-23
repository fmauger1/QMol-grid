function [E,DE] = SE_energyWaveFunction(obj,V,wfcn)
%SE_energyWaveFunction computes the wave-function energies, defined as
%   E = < psi | H | psi > and error DE = || H | psi > - E | psi > ||
    
    
    % Initialization
    E                   =   NaN(size(wfcn.wfcn,2),1);
    DE                  =   NaN(size(wfcn.wfcn,2),1);
    
    % Parse spin-up orbitals
    for k = 1:size(wfcn.wfcn,2)
        % Orbital energy
        Hp          =   obj.SE_operatorHamiltonian(V,wfcn.wfcn(:,k),0);
        E(k)        =   real(sum(conj(wfcn.wfcn(:,k)) .* Hp));
        
        % Orbital error
        if nargout > 1
            DE(k)   =   sqrt(sum( abs( Hp - E(k)*wfcn.wfcn(:,k) ).^2 ));
        end
    end

end