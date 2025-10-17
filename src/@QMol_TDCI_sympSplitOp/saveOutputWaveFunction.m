function saveOutputWaveFunction(obj,k,t)
%saveOutputWaveFunction

    % Wave functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oWfcn.ind(obj.oWfcn.n)
        % Save wave functions
        ind             =   (obj.oWfcn.n-1)*prod(obj.oWfcn.shape)+1:obj.oWfcn.n*prod(obj.oWfcn.shape);
                            
        if isreal(obj.CI.wfcn), obj.oWfcn.waveFunction(ind) =   complex(obj.CI.wfcn(:,    obj.oWfcn.indexWaveFunction));
        else,                   obj.oWfcn.waveFunction(ind) =           obj.CI.wfcn(:,    obj.oWfcn.indexWaveFunction);     end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oWfcn',k,t);        end

        % Update counter
        obj.oWfcn.n     =   obj.oWfcn.n + 1;
    end

    
end