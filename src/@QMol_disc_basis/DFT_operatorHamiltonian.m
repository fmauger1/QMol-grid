function Hp = DFT_operatorHamiltonian(obj,V,p,isUp,~)
%DFT_operatorHamiltonian DFT-Hamiltonian operator H | psi >
    
    % Apply Hamiltonian operator
    if isempty(obj.mTv), if obj.QM.isSpinPol,   Hp  =   obj.mT *p + V.applyPotential(p,isUp);       %#ok<ALIGN> 
                         else,                  Hp  =   obj.mT *p + V.applyPotential(p);            end
    else,                if obj.QM.isSpinPol,   Hp  =   obj.mTv*p + V.applyPotential(p,isUp);
                        else,                   Hp  =   obj.mTv*p + V.applyPotential(p);            end, end
end