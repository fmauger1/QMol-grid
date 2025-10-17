function saveOutputResults(obj,k,t)
%saveOutputResults
    
    % Save Schrodinger-equation object in MATLAB file ~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oCI.ind(obj.oCI.n)
        obj.saveOutputCI(t);
    end

    % Save dipole/velocity/acceleration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oDip.ind(obj.oDip.n)   ||   k == obj.oVel.ind(obj.oVel.n)   ||   k == obj.oAcc.ind(obj.oAcc.n)
        obj.saveOutputDipole(k,t);       
    end

    % Save DFT and orbital energies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oECI.ind(obj.oECI.n)   ||   k == obj.oEWfcn.ind(obj.oEWfcn.n)
        obj.saveOutputEnergy(k,t);
    end

    % Save ionization statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oIon.ind(obj.oIon.n)
        obj.saveOutputIonization(k,t);
    end

    % Save wave functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oWfcn.ind(obj.oWfcn.n)
        obj.saveOutputWaveFunction(k,t);
    end

    % Output functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oF.ind(obj.oF.n)
        obj.saveOutputFcn(k,t);
    end

    % Additional implementation-specific output ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.saveOutputChildren(k,t);
end