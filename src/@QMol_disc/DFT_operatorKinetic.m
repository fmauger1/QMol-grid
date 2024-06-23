function Hp = DFT_operatorKinetic(obj,p,S)
%DFT_operatorKinetic kinetic operator T | psi > 
    
    % Initialization
    if nargin < 3,  S   =   0;  end

    % Kinetic operator
    if isempty(obj.Tv),     Hp  =   ifft(obj.T .* fft(p));
    else,                   Hp  =   ifft(obj.Tv.* fft(p));                  end

    % Symmetry (only is needed)
    if S ~= 0,              Hp  =   obj.DFT_applySymmetry(S,Hp);            end
end