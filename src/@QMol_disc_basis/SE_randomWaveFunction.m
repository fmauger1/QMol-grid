function p = SE_randomWaveFunction(obj,randStr,~)
%SE_randomWaveFunction returns a random Schrodinger-equation wave function, 
%   using the input random-generator stream, within the specified symmetry 
%   group

    % Generate random wave functions
    p                   =   rand(randStr,obj.SE_sizeWaveFunction)-.5;       % Random wave function
    p                   =   obj.SE_normalizeWaveFunction(p);                % Normalize the wave function
end