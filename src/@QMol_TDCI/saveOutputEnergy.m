function saveOutputEnergy(obj,k,t)
%saveOutputEnergy
    
    % Common energy components
    if k == obj.oEWfcn.ind(obj.oEWfcn.n)   ||   k == obj.oECI.ind(obj.oECI.n)
        E               =   obj.CI.getEnergy;
    end

    % Save the wave function energies
    if k == obj.oEWfcn.ind(obj.oEWfcn.n)
        % Wave function energies
        obj.oEWfcn.waveFunction(:,obj.oEWfcn.n)    =   E;
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oEWfcn',k,t);       end

        % Update counter
        obj.oEWfcn.n    =   obj.oEWfcn.n + 1;
    end

    % Save the configuration interaction energies
    if k == obj.oECI.ind(obj.oECI.n)
        % Hamiltonian operator energy
        obj.oECI.Hamiltonian(obj.oECI.n)        =   sum(E);
                                                                            
        % External-field energies
       [obj.oECI.externalField(obj.oECI.n), ...
        obj.oECI.autonomization(obj.oECI.n)]    =   obj.getExternalFieldEnergy(t);
        
        % Add external field
        if obj.sECI,        obj.addOutputExternalField('oECI',k,t);         end

        % Update total energy
        obj.oECI.total(obj.oECI.n)              =   obj.oECI.Hamiltonian(obj.oECI.n) + ...
                                                    obj.oECI.externalField(obj.oECI.n) + ...
                                                    obj.oECI.autonomization(obj.oECI.n);

        % Update counter
        obj.oECI.n     =   obj.oECI.n + 1;
    end
    
end