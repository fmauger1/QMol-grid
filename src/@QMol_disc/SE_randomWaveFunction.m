function p = SE_randomWaveFunction(obj,randStr,S)
%SE_randomWaveFunction returns a random Schrodinger-equation wave function, 
%   using the input random-generator stream, within the specified symmetry 
%   group
    
    % Initialization
    if nargin < 3,  S   =   0;      end

    % Generate random orbitals
    p                   =   rand(randStr,obj.SE_sizeWaveFunction)-.5;       % Random wave function
    p                   =   obj.SE_applySymmetry(S,p);                      % Apply symmetry
    p                   =   obj.SE_normalizeWaveFunction(p);                % Normalize the wave function
end