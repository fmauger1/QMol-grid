function saveOutputEnergy(obj,k,t)
%saveOutputEnergy
    
    % Save the wave function energies
    if k == obj.oEWfcn.ind(obj.oEWfcn.n)
        % Wave function energies
        E               =   obj.SE.getEnergy('Wfcn');
        obj.oEWfcn.waveFunction(:,obj.oEWfcn.n)    =   E;
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oEWfcn',k,t);       end

        % Update counter
        obj.oEWfcn.n    =   obj.oEWfcn.n + 1;
    end

    % Save the Schrodinger-equation energies
    if k == obj.oESE.ind(obj.oESE.n)
        % Schrodinger-equation-model energies
       [obj.oESE.total(           obj.oESE.n), ...
        obj.oESE.kinetic(       :,obj.oESE.n), ...
        obj.oESE.potential(     :,obj.oESE.n)]   =   obj.SE.getEnergy('SE');   % Density still good
                                                                            
        % Esternal-field energies
       [obj.oESE.externalField( :,obj.oESE.n), ...
        obj.oESE.autonomization(  obj.oESE.n)]  =   obj.getExternalFieldEnergy(t);
        
        % Add external field
        if obj.sESE,        obj.addOutputExternalField('oESE',k,t);         end

        % Update total energy
        obj.oESE.total(obj.oESE.n)              =   obj.oESE.total(obj.oESE.n) + ...
                                                    sum(obj.oESE.externalField(:,obj.oESE.n)) + ...
                                                    obj.oESE.autonomization(     obj.oESE.n);

        % Update counter
        obj.oESE.n     =   obj.oESE.n + 1;
    end
    
end