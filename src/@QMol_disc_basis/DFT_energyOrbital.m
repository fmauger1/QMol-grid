function [E,DE] = DFT_energyOrbital(obj,V,KSO)
%DFT_orbitalEnergy computes the Kohn-Sham orbital energy, defined as
%   E = < psi | H | psi > and error DE = || H | psi > - E | psi > ||
    
    % Spin polarized -----------------------
    if KSO.isSpinPol
        % Initialization
        Eup     =   NaN(size(KSO.KSOup,2),1);       DEup    =   NaN(size(KSO.KSOup,2),1);
        Edw     =   NaN(size(KSO.KSOdw,2),1);       DEdw    =   NaN(size(KSO.KSOdw,2),1);
        
        % Parse spin-up orbitals
        for k = 1:size(KSO.KSOup,2)
            % Orbital energy
            Hp          =   obj.DFT_operatorHamiltonian(V,KSO.KSOup(:,k),true,0);
            Eup(k)      =   real(sum(conj(KSO.KSOup(:,k)) .* Hp));
            
            % Orbital error
            if nargout > 1
                DEup(k) =   sqrt(sum( abs( Hp - Eup(k)*KSO.KSOup(:,k) ).^2 ));
            end
        end
        
        % Parse spin-down orbitals
        for k = 1:size(KSO.KSOdw,2)
            % Orbital energy
            Hp          =   obj.DFT_operatorHamiltonian(V,KSO.KSOdw(:,k),false,0);
            Edw(k)      =   real(sum(conj(KSO.KSOdw(:,k)) .* Hp));
            
            % Orbital error
            if nargout > 1
                DEdw(k) =   sqrt(sum( abs( Hp - Edw(k)*KSO.KSOdw(:,k) ).^2 ));
            end
        end

        % Put it all together
        E   =   {Eup, Edw};     DE  =   {DEup,DEdw};

    % Spin restricted ----------------------
    else
        % Initialization
        E   =   NaN(size(KSO.KSO,2),1);             DE      =   NaN(size(KSO.KSO,2),1);
        
        % Parse spin-up orbitals
        for k = 1:size(KSO.KSO,2)
            % Orbital energy
            Hp          =   obj.DFT_operatorHamiltonian(V,KSO.KSO(:,k),0);
            E(k)        =   real(sum(conj(KSO.KSO(:,k)) .* Hp));
            
            % Orbital error
            if nargout > 1
                DE(k)   =   sqrt(sum( abs( Hp - E(k)*KSO.KSO(:,k) ).^2 ));
            end
        end
    end
end