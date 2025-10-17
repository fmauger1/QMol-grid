function cleanConfigurationBasis(obj)
%cleanConfigurationBasis cleans the configuration state basis
%   Use cleanConfigurationBasis to remove duplicate configuration states in
%   configurationBasis, including those related by a permutation of the
%   spin orbitals in the index vector. The cleaning process preserves the
%   order of (i) configuration states in the basis, keeping only the first
%   unique instance of each state, and (ii) the order of the spin-orbital
%   withing each configuration state. For instance, the configurationBasis
%       [-1 -2 1 2; -1 -3 1 2; -2 -1 1 2;-1 -2 1 3]
%   is cleaned to
%       [-1 -2 1 2; -1 -3 1 2; -1 -2 1 3]
%   since [-2 -1 1 2] is a duplicate of [-1 -2 1 2], and thus removed from
%   the set.

    % Sort configuration states to remove permutations
    CSB                 =   sort(obj.CSB,2);

    % Identify unique configuration states
    [~,ind]             =   unique(CSB,'rows');

    % Clean the configuration-state basis
    obj.CSB             =   obj.CSB(sort(ind),:);
end

