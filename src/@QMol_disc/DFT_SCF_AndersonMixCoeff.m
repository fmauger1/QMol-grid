function b = DFT_SCF_AndersonMixCoeff(obj,X1i,X2i,X1o,X2o)
%DFT_SCF_AndersonMixCoeff computes the mixing coefficient for the SCF
%   Anderson's scheme
    
    if isa(X1i,'QMol_DFT_density')
        if X1i.isSpinPol,   b   =  [sum((X1i.rhoUp-X1o.rhoUp).*(X1i.rhoUp-X1o.rhoUp-X2i.rhoUp+X2o.rhoUp))       ...
                                    sum((X1i.rhoDw-X1o.rhoDw).*(X1i.rhoDw-X1o.rhoDw-X2i.rhoDw+X2o.rhoDw))] ./   ...
                                   [sum( (X1i.rhoUp-X1o.rhoUp-X2i.rhoUp+X2o.rhoUp).^2 )                         ...
                                    sum( (X1i.rhoDw-X1o.rhoDw-X2i.rhoDw+X2o.rhoDw).^2 )];                       %#ok<ALIGN> 
        else,               b   =   sum((X1i.rho-X1o.rho).*(X1i.rho-X1o.rho-X2i.rho+X2o.rho)) /                 ...
                                    sum( (X1i.rho-X1o.rho-X2i.rho+X2o.rho).^2 );                                end
    elseif isa(X1i,'QMol_DFT_Vks')  % Only the explicit part of the potential is mixed
        if obj.QM.isSpinPol,b   =  [sum((X1i.Vup-X1o.Vup).*(X1i.Vup-X1o.Vup-X2i.Vup+X2o.Vup))       ...
                                    sum((X1i.Vdw-X1o.Vdw).*(X1i.Vdw-X1o.Vdw-X2i.Vdw+X2o.Vdw))] ./   ...
                                   [sum( (X1i.Vup-X1o.Vup-X2i.Vup+X2o.Vup).^2 )                     ...
                                    sum( (X1i.Vdw-X1o.Vdw-X2i.Vdw+X2o.Vdw).^2 )];                   %#ok<ALIGN> 
        else,               b   =   sum((X1i.V-X1o.V).*(X1i.V-X1o.V-X2i.V+X2o.V)) /                 ...
                                    sum( (X1i.V-X1o.V-X2i.V+X2o.V).^2 );                            end
    else
        if obj.QM.isSpinPol,b   =   NaN(1,2);
        else,               b   =   NaN;                                    end
        warning('QMol:disc:DFT_SCF_AndersonMixCoeff', ...
                ['Anderson''s mixing coefficient is not defined for ' class(X) ' objects.'])
    end
        
end