function Hp = SE_operatorHamiltonian(obj,V,p,S)
%SE_operatorHamiltonian Schrodinger-equation-Hamiltonian operator H | psi >
    
    % Apply Hamiltonian operator
    if isempty(obj.Tv)                                                      % SE_operatorHamiltonian(obj,V,p,S)
        if nargin < 4,  Hp  =   ifft(obj.T .* fft(p)) + V.V.*p;    S = 0;
        else,           Hp  =   ifft(obj.T .* fft(p)) + V.V.*p;    end
    else                                                                    % SE_operatorHamiltonian(obj,V,p,S)
        if nargin < 4,  Hp  =   ifft(obj.Tv.* fft(p)) + V.V.*p;    S = 0;
        else,           Hp  =   ifft(obj.Tv.* fft(p)) + V.V.*p;    end
    end

    % Symmetry (only is needed)
    if S ~= 0,              Hp  =   obj.SE_applySymmetry(S,Hp);             end
end