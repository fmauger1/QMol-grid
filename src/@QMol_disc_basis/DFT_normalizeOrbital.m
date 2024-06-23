function p = DFT_normalizeOrbital(~,p)
%DFT_normalizeOrbital normalizes the input (group) of Kohn-Sham orbitals
    
    for k = 1:size(p,2)
        p(:,k)          =   p(:,k) / sqrt(sum(abs(p(:,k)).^2));
    end
end