function Hp = SE_operatorKinetic(obj,p,~)
%DFT_operatorKinetic kinetic operator T | psi > 

    % Kinetic operator
    if isempty(obj.Tv),     Hp  =   obj.mT *p;
    else,                   Hp  =   obj.mTv*p;                              end
end