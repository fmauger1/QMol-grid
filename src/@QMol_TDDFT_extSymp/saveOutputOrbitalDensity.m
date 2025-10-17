function saveOutputOrbitalDensity(obj,k,t)
%saveOutputOrbitalDensity

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.setOutputOrbital;

    % One-body density ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oRho.ind(obj.oRho.n)
        % Save density
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);

        ind             =   (obj.oRho.n-1)*prod(obj.oRho.shape)+1:obj.oRho.n*prod(obj.oRho.shape);
        if obj.DFT.isSpinPol,   obj.oRho.totalUp(ind)       =   obj.DFT.rho.rhoUp;  %#ok<ALIGN> 
                                obj.oRho.totalDown(ind)     =   obj.DFT.rho.rhoDw;
        else,                   obj.oRho.total(ind)         =   obj.DFT.rho.rho;    end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oRho',k,t);        end

        % Update counter
        obj.oRho.n      =   obj.oRho.n + 1;

    end
    % Kohn-Sham orbitals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oKSO.ind(obj.oKSO.n),                                       if obj.DFT.isSpinPol
        % Spin up
        ind             =   (obj.oKSO.n-1)*prod(obj.oKSO.shape{1})+1:obj.oKSO.n*prod(obj.oKSO.shape{1});
                            
        if obj.DFT.dim == 1,        obj.oKSO.orbitalUp(ind)     =   complex(obj.DFT.KSO.KSOup(:,    obj.oKSO.indexOrbitalUp));   %#ok<ALIGN> 
        elseif obj.DFT.dim == 2,    obj.oKSO.orbitalUp(ind)     =   complex(obj.DFT.KSO.KSOup(:,:,  obj.oKSO.indexOrbitalUp));
        elseif obj.DFT.dim == 3,    obj.oKSO.orbitalUp(ind)     =   complex(obj.DFT.KSO.KSOup(:,:,:,obj.oKSO.indexOrbitalUp));   end

        % Spin down
        ind             =   (obj.oKSO.n-1)*prod(obj.oKSO.shape{2})+1:obj.oKSO.n*prod(obj.oKSO.shape{2});
                            
        if obj.DFT.dim == 1,        obj.oKSO.orbitalDown(ind)   =   complex(obj.DFT.KSO.KSOdw(:,    obj.oKSO.indexOrbitalDown)); %#ok<ALIGN> 
        elseif obj.DFT.dim == 2,    obj.oKSO.orbitalDown(ind)   =   complex(obj.DFT.KSO.KSOdw(:,:,  obj.oKSO.indexOrbitalDown));
        elseif obj.DFT.dim == 3,    obj.oKSO.orbitalDown(ind)   =   complex(obj.DFT.KSO.KSOdw(:,:,:,obj.oKSO.indexOrbitalDown)); end
                                                                            else
        % Spin restricted
        ind             =   (obj.oKSO.n-1)*prod(obj.oKSO.shape)+1:obj.oKSO.n*prod(obj.oKSO.shape);
                            
        if obj.DFT.dim == 1,        obj.oKSO.orbital(ind)       =   complex(obj.DFT.KSO.KSO(:,    obj.oKSO.indexOrbital));       %#ok<ALIGN> 
        elseif obj.DFT.dim == 2,    obj.oKSO.orbital(ind)       =   complex(obj.DFT.KSO.KSO(:,:,  obj.oKSO.indexOrbital));
        elseif obj.DFT.dim == 3,    obj.oKSO.orbital(ind)       =   complex(obj.DFT.KSO.KSO(:,:,:,obj.oKSO.indexOrbital));       end
                                                                            end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oKSO',k,t);         end

        % Update counter
        obj.oKSO.n      =   obj.oKSO.n + 1;
    end

    % Kohn-Sham orbital projections ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oKSOP.ind(obj.oKSOP.n),                                     if obj.DFT.isSpinPol
        % Spin up
        ind             =   (obj.oKSOP.n-1)*prod(obj.oKSOP.shape{1})+1:obj.oKSOP.n*prod(obj.oKSOP.shape{1});

        if iscell(obj.oKSOP.basis),     obj.pKSO    =   obj.oKSOP.basis{1}.DFT_projectOrbital(obj.DFT.KSO,obj.pKSO);
        else,                           obj.pKSO    =   obj.oKSOP.basis.DFT_projectOrbital(   obj.DFT.KSO,obj.pKSO);    end

        obj.oKSOP.orbitalUp(ind)    =   complex(obj.pKSO.KSOup(:,obj.oKSOP.indexOrbitalUp));

        % Spin up
        ind             =   (obj.oKSOP.n-1)*prod(obj.oKSOP.shape{2})+1:obj.oKSOP.n*prod(obj.oKSOP.shape{2});

        if iscell(obj.oKSOP.basis),     obj.pKSO    =   obj.oKSOP.basis{2}.DFT_projectOrbital(obj.DFT.KSO,obj.pKSO);    end

        obj.oKSOP.orbitalDown(ind)  =   complex(obj.pKSO.KSOdw(:,obj.oKSOP.indexOrbitalDown));
                                                                            else
        % Spin restricted
        ind             =   (obj.oKSOP.n-1)*prod(obj.oKSOP.shape)+1:obj.oKSOP.n*prod(obj.oKSOP.shape);

        obj.pKSO        =   obj.oKSOP.basis.DFT_projectOrbital(obj.DFT.KSO,obj.pKSO);
        obj.oKSOP.orbital(ind)      =   complex(obj.pKSO.KSO(:,obj.oKSOP.indexOrbital));
                                                                            end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oKSOP',k,t);        end

        % Update counter
        obj.oKSOP.n     =   obj.oKSOP.n + 1;
    end
    
end