function d = DFT_dist(obj,X1,X2)
%DFT_dist computes the L2-norm distance between the input objects
    
    if isa(X1,'QMol_DFT_density')
        % Distance between two densities
        if X1.isSpinPol,    d   =  [sum((X1.rhoUp-X2.rhoUp).^2), ...
                                    sum((X1.rhoDw-X2.rhoDw).^2)];           %#ok<ALIGN> 
        else,               d   =   sum((X1.rho  -X2.rho  ).^2);            end

        d               =   sqrt(d*obj.dx);
    else
        if obj.QM.isSpinPol,d   =   NaN(1,2);
        else,               d   =   NaN;                                    end
        warning('QMol:disc:DFT_computeDistance', ...
                ['Distance is not defined for ' class(X) ' objects.'])
    end
end