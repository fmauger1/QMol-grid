function p = SE_normalizeWaveFunction(obj,p)
%SE_normalizeWaveFunction normalizes the input (group) of wave function(s)
    
    for k = 1:size(p,2)
        p(:,k)          =   p(:,k) / sqrt(sum(abs(p(:,k)).^2)*obj.dx);
    end
end