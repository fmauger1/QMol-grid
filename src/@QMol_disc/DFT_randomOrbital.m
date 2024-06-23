function p = DFT_randomOrbital(obj,randStr,S)
%DFT_randomOrbital returns a random DFT Kohn-Sham orbitals, using the input
%   random-generator stream, within the specified symmetry group
    
    % Initialization
    if nargin < 3,  S   =   0;      end

    % Generate random orbitals
    p                   =   rand(randStr,obj.DFT_sizeOrbital)-.5;           % Random orbital
    p                   =   obj.DFT_applySymmetry(S,p);                     % Apply symmetry
    p                   =   obj.DFT_normalizeOrbital(p);                    % Normalize the orbital
end