function s = SE_sizeWaveFunction(obj)
%SE_sizeWaveFunction returns the size of a Schrodinger-equation wave
%   function
    
    s                   =   [size(obj.v,2) 1];
end