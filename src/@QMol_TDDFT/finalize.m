function finalize(obj)
%finalize
    
    % Clean output
    obj.setOutputDFT('clean');                                              % DFT object in separate files
    obj.setOutputDipole('clean');                                           % dipole/velocity/acceleration
    obj.setOutputEnergy('clean');                                           % DFT and orbital energies
    obj.setOutputIonization('clean');                                       % Ionization
    obj.setOutputOrbitalDensity('clean');                                   % Orbital/density
    obj.setOutputFunction('clean');                                         % Output function
    obj.setOutputRestart('clean');                                          % Restart file

    obj.setOutputChildren('clean');                                         % Implementation-specific output

end