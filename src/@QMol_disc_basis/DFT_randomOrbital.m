function p = DFT_randomOrbital(obj,randStr,~)
%DFT_randomOrbital returns a random DFT Kohn-Sham orbitals, using the input
%   random-generator stream, within the specified symmetry group

    % Generate random orbitals
    p                   =   rand(randStr,obj.DFT_sizeOrbital)-.5;           % Random orbital
    p                   =   obj.DFT_normalizeOrbital(p);                    % Normalize the orbital
end