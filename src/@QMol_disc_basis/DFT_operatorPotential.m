function Hp = DFT_operatorPotential(~,V,p,isUp,~)
%DFT_operatorPotential potential operator V | psi >
    
    % Apply potential operator
    if V.isSpinPol,     Hp  =   V.applyPotential(p,isUp);
    else,               Hp  =   V.applyPotential(p);                        end

end