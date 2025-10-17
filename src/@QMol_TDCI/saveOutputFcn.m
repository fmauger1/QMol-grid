function saveOutputFcn(obj,k,t)
%saveOutputFcn
    
    % Save output function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if k == obj.oF.ind(obj.oF.n)
        % Output function
        ind             =   (obj.oF.n-1)*prod(obj.oF.shape)+1:obj.oF.n*prod(obj.oF.shape);
        obj.oF.result(ind)  =   obj.sF(obj.CI,t);
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oF',k,t);           end

        % Update counter
        obj.oF.n        =   obj.oF.n + 1;
    end
    
end