function Hp = SE_operatorHamiltonian(obj,V,p,~)
%SE_operatorHamiltonian Schrodinger-equation-Hamiltonian operator H | psi >
    
    % Kinetic operator
    if isempty(obj.Tv),     Hp  =   obj.mT *p;
    else,                   Hp  =   obj.mTv*p;                              end

    % Apply potential operator
    Hp                  =   Hp + V.mV*p;
    
end