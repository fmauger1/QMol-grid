function finalize(obj)
%finalize
    
    % Clean output
    obj.setOutputCI('clean');                                               % CI object in separate files
    obj.setOutputDipole('clean');                                           % dipole/velocity/acceleration
    obj.setOutputEnergy('clean');                                           % CI and wave function energies
    obj.setOutputIonization('clean');                                       % Ionization
    obj.setOutputWaveFunction('clean');                                     % Wave function(s)
    obj.setOutputFunction('clean');                                         % Output function
    obj.setOutputRestart('clean');                                          % Restart file

    obj.setOutputChildren('clean');                                         % Implementation-specific output

end