function setOutputFunction(obj,opt)
%setOutputFunction 

switch lower(opt)
case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Clean up any old data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.oFRho           =   [];
    obj.oFKSO           =   [];
    
    % Output function of the density ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if isa(obj.sFRho,'function_handle')
        % When to save
        obj.oFRho.ind   =   [obj.getOutputIndex(obj.sFRhoT),NaN];
        obj.oFRho.time  =   obj.tspan(obj.oFRho.ind(1:end-1));
        obj.oFRho.n     =   1;

        if obj.sEF,         obj.setOutputExternalField('oFRho','init');     end

        % Run output function on initial condition (to get the shape)
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);
        R               =   obj.sFRho(obj.DFT.rho,obj.tspan(1));

        S               =   size(R);
        obj.oFRho.shape =   S;
        if numel(S) > 1   &&   S(end) == 1,     S   =   S(1:end-1);         end

        obj.oFRho.result=   NaN([S numel(obj.oFRho.time)]);

        % Start time output to be kept?
        if obj.oFRho.ind(1) == 1
            % Copy results
            obj.oFRho.result(1:prod(obj.oFRho.shape))   =   R;

            % Add external field
            if obj.sEF,         obj.addOutputExternalField('oFRho',1,obj.tspan(1)); end

            % Update next to save
            if isscalar(obj.oFRho.ind),     obj.oFRho.n =  -1;              % done with output function of density
            else,                           obj.oFRho.n =   2;              end

        end
    else
        obj.oFRho.shape =   [];
        obj.oFRho.ind   =   NaN;
        obj.oFRho.n     =   1;
    end
    
    % Output function of the Kohn-Sham orbitals ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if isa(obj.sFKSO,'function_handle')
        % When to save
        obj.oFKSO.ind   =   [obj.getOutputIndex(obj.sFKSOT),NaN];
        obj.oFKSO.time  =   obj.tspan(obj.oFKSO.ind(1:end-1));
        obj.oFKSO.n     =   1;

        if obj.sEF,         obj.setOutputExternalField('oFKSO','init');     end

        % Run output function on initial condition
        R               =   obj.sFKSO(obj.DFT.KSO,obj.tspan(1));

        S               =   size(R);
        obj.oFKSO.shape =   S;
        if numel(S) > 1   &&   S(end) == 1,     S   =   S(1:end-1);         end

        obj.oFKSO.result=   NaN([S numel(obj.oFKSO.time)]);

        % Start time output to be kept?
        if obj.oFKSO.ind(1) == 1
            % Copy results
            obj.oFKSO.result(1:prod(obj.oFKSO.shape))   =   R;

            % Add external field
            if obj.sEF,         obj.addOutputExternalField('oFKSO',1,obj.tspan(1)); end

            % Update next to save
            if isscalar(obj.oFKSO.ind),     obj.oFKSO.n =  -1;              % done with output function of density
            else,                           obj.oFKSO.n =   2;              end

        end
    else
        obj.oFKSO.shape =   [];
        obj.oFKSO.ind   =   NaN;
        obj.oFKSO.n     =   1;
    end
    
case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main components
    if isa(obj.sFRho,'function_handle')                                     %#ok<ALIGN> 
            obj.oFRho   =   rmfield(obj.oFRho,{'ind','n','shape'});         if obj.sEF
            obj.setOutputExternalField('oFRho','clean');                    end
    else,   obj.oFRho   =   [];                                             end
    if isa(obj.sFKSO,'function_handle')                                     %#ok<ALIGN> 
            obj.oFKSO   =   rmfield(obj.oFKSO,{'ind','n','shape'});         if obj.sEF
            obj.setOutputExternalField('oFKSO','clean');                    end
    else,   obj.oFKSO   =   [];                                             end

otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unexpected option
    error('QMol:TDDFT:setOutputfunction',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end