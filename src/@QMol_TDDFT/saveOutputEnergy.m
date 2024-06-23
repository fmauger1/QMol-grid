function saveOutputEnergy(obj,k,t)
%saveOutputEnergy
    
    % Save the orbital energies
    if k == obj.oEKSO.ind(obj.oEKSO.n)
        % Orbital energies
        E               =   obj.DFT.getEnergy('KSO');                       if obj.DFT.isSpinPol
        obj.oEKSO.orbitalUp(  :,obj.oEKSO.n)    =   E{1};
        obj.oEKSO.orbitalDown(:,obj.oEKSO.n)    =   E{2};                   else
        obj.oEKSO.orbital(    :,obj.oEKSO.n)    =   E;                      end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oEKSO',k,t);        end

        % Update counter
        obj.oEKSO.n     =   obj.oEKSO.n + 1;
    end

    % Save the DFT energies
    if k == obj.oEDFT.ind(obj.oEDFT.n),                                     if obj.sEKSO
        % DFT-model energies
       [obj.oEDFT.total(              obj.oEDFT.n), ...
        obj.oEDFT.kinetic(          :,obj.oEDFT.n), ...
        obj.oEDFT.external(         :,obj.oEDFT.n), ...
        obj.oEDFT.Hartree(            obj.oEDFT.n), ...
        obj.oEDFT.exchangeCorrelation(obj.oEDFT.n)] =   obj.DFT.getEnergy('DFT',obj.DFT.rho);   % Density still good
                                                                            else
       [obj.oEDFT.total(              obj.oEDFT.n), ...
        obj.oEDFT.kinetic(          :,obj.oEDFT.n), ...
        obj.oEDFT.external(         :,obj.oEDFT.n), ...
        obj.oEDFT.Hartree(            obj.oEDFT.n), ...
        obj.oEDFT.exchangeCorrelation(obj.oEDFT.n)] =   obj.DFT.getEnergy('DFT');
                                                                            end
        % Esternal-field energies
       [obj.oEDFT.externalField(    :,obj.oEDFT.n), ...
        obj.oEDFT.autonomization(     obj.oEDFT.n)] =   obj.getExternalFieldEnergy(t);
        
        % Add external field
        if obj.sEDFT,       obj.addOutputExternalField('oEDFT',k,t);        end

        % Update total energy
        obj.oEDFT.total(obj.oEDFT.n)    =   obj.oEDFT.total(obj.oEDFT.n) + ...
                                            sum(obj.oEDFT.externalField(:,obj.oEDFT.n)) + ...
                                            obj.oEDFT.autonomization(     obj.oEDFT.n);

        % Update counter
        obj.oEDFT.n     =   obj.oEDFT.n + 1;
    end
    
end