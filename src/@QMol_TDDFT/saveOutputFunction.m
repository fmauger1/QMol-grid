function saveOutputFunction(obj,k,t)
%saveOutputFunction
    
    % Save output function of the density ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oFRho.ind(obj.oFRho.n)
        % Output function
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);

        ind             =   (obj.oFRho.n-1)*prod(obj.oFRho.shape)+1:obj.oFRho.n*prod(obj.oFRho.shape);
        obj.oFRho.result(ind)   =   obj.sFRho(obj.DFT.rho,t);
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oFRho',k,t);        end

        % Update counter
        obj.oFRho.n     =   obj.oFRho.n + 1;
    end

    % Save output function of the Kohn-Sham orbitals ~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oFKSO.ind(obj.oFKSO.n)
        % Output function
        ind             =   (obj.oFKSO.n-1)*prod(obj.oFKSO.shape)+1:obj.oFKSO.n*prod(obj.oFKSO.shape);
        obj.oFKSO.result(ind)   =   obj.sFKSO(obj.DFT.KSO,t);
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oFKSO',k,t);        end

        % Update counter
        obj.oFKSO.n     =   obj.oFKSO.n + 1;
    end
    
end