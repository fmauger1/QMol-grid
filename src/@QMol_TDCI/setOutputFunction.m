function setOutputFunction(obj,opt)
%setOutputIonization 

switch lower(opt)
case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Clean up any old data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.oF              =   [];
    
    % Output function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if isa(obj.sF,'function_handle')
        % When to save
        obj.oF.ind      =   [obj.getOutputIndex(obj.sFT),NaN];
        obj.oF.time     =   obj.tspan(obj.oF.ind(1:end-1));
        obj.oF.n        =   1;

        if obj.sEF,         obj.setOutputExternalField('oF','init');        end

        % Run output function on initial condition
        R               =   obj.sF(obj.CI,obj.tspan(1));

        S               =   size(R);
        obj.oF.shape    =   S;
        if numel(S) > 1   &&   S(end) == 1,     S   =   S(1:end-1);         end

        obj.oF.result   =   NaN([S numel(obj.oF.time)]);

        % Start time output to be kept?
        if obj.oF.ind(1) == 1
            % Copy results
            obj.oF.result(1:prod(obj.oF.shape))   =   R;

            % Add external field
            if obj.sEF,         obj.addOutputExternalField('oF',1,obj.tspan(1)); end

            % Update next to save
            if isscalar(obj.oF.ind),        obj.oF.n    =  -1;              % done with output function
            else,                           obj.oF.n    =   2;              end

        end
    else
        obj.oF.shape    =   [];
        obj.oF.ind      =   NaN;
        obj.oF.n        =   1;
    end
    
case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main components
    if isa(obj.sF,'function_handle')                                        %#ok<ALIGN> 
            obj.oF      =   rmfield(obj.oF,{'ind','n','shape'});            if obj.sEF
            obj.setOutputExternalField('oF','clean');                       end
    else,   obj.oF      =   [];                                             end

otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unexpected option
    error('QMol:TDCI:setOutputfunction',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end