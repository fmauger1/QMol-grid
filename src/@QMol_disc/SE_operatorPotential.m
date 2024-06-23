function Hp = SE_operatorPotential(obj,V,p,S)
%SE_operatorPotential potential operator V | psi >
    
    % Apply potential operator
    if nargin < 4,  Hp  =   V.V.*p;         S = 0;
    else,           Hp  =   V.V.*p;         end

    % Symmetry (only is needed)
    if S ~= 0,              Hp  =   obj.DFT_applySymmetry(S,Hp);            end
end