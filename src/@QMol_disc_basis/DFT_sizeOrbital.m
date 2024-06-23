function s = DFT_sizeOrbital(obj)
%DFT_sizeOrbital returns the size of a DFT Kohn-Sham orbital
    
    s                   =   [size(obj.v,2) 1];
end