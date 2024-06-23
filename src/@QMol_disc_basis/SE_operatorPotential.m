function Hp = SE_operatorPotential(~,V,p,~)
%SE_operatorPotential potential operator V | psi >
    
    % Apply potential operator
    Hp                  =   V.mV*p;
    
end