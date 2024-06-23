function saveOutputResults(obj,k,t)
%saveOutputResults
    
    % Save DFT object in MATLAB file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oDFT.ind(obj.oDFT.n)
        obj.saveOutputDFT(t);
    end

    % Save dipole/velocity/acceleration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oDip.ind(obj.oDip.n)   ||   k == obj.oVel.ind(obj.oVel.n)   ||   k == obj.oAcc.ind(obj.oAcc.n)
        obj.saveOutputDipole(k,t);       
    end

    % Save DFT and orbital energies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oEDFT.ind(obj.oEDFT.n)   ||   k == obj.oEKSO.ind(obj.oEKSO.n)
        obj.saveOutputEnergy(k,t);
    end

    % Save ionization statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oIon.ind(obj.oIon.n)
        obj.saveOutputIonization(k,t);
    end

    % Save orbital/density ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oRho.ind(obj.oRho.n)   ||   k == obj.oKSO.ind(obj.oKSO.n)   || k == obj.oKSOP.ind(obj.oKSOP.n)
        obj.saveOutputOrbitalDensity(k,t);
    end

    % Output functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oFRho.ind(obj.oFRho.n)   ||   k == obj.oFKSO.ind(obj.oFKSO.n)
        obj.saveOutputFunction(k,t);
    end

    % Additional implementation-specific output ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.saveOutputChildren(k,t);
end