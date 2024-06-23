function saveOutputWaveFunction(obj,k,t)
%saveOutputWaveFunction

    % Wave functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oWfcn.ind(obj.oWfcn.n)
        % Save wave functions
        ind             =   (obj.oWfcn.n-1)*prod(obj.oWfcn.shape)+1:obj.oWfcn.n*prod(obj.oWfcn.shape);
                            
        if obj.SE.dim == 1,         obj.oWfcn.waveFunction(ind) =   complex(obj.SE.wfcn.wfcn(:,    obj.oWfcn.indexWaveFunction));       %#ok<ALIGN> 
        elseif obj.SE.dim == 2,     obj.oWfcn.waveFunction(ind) =   complex(obj.SE.wfcn.wfcn(:,:,  obj.oWfcn.indexWaveFunction));
        elseif obj.SE.dim == 3,     obj.oWfcn.waveFunction(ind) =   complex(obj.SE.wfcn.wfcn(:,:,:,obj.oWfcn.indexWaveFunction));       end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oWfcn',k,t);        end

        % Update counter
        obj.oWfcn.n     =   obj.oWfcn.n + 1;
    end

    % Wave function projections ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oWfcnP.ind(obj.oWfcnP.n)
        % Save wave function projection
        ind             =   (obj.oWfcnP.n-1)*prod(obj.oWfcnP.shape)+1:obj.oWfcnP.n*prod(obj.oWfcnP.shape);

        obj.pWfcn       =   obj.oWfcnP.basis.SE_projectWaveFunction(obj.SE.wfcn,obj.pWfcn);
        obj.oWfcnP.waveFunction(ind)    =   complex(obj.pWfcn.wfcn(:,obj.oWfcnP.indexWaveFunction));
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oWfcnP',k,t);       end

        % Update counter
        obj.oWfcnP.n    =   obj.oWfcnP.n + 1;
    end
    
end