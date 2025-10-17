function setOutputEnergy(obj,opt)
%setOutputEnergy

switch lower(opt)
case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Clean up any old data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.oECI            =   [];
    obj.oEWfcn          =   [];
    
    % Schrodinger-equation energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sECI
        % When to save the DFT energies
        obj.oECI.ind    =   [obj.getOutputIndex(obj.sECIT),NaN];
        obj.oECI.time   =   obj.tspan(obj.oECI.ind(1:end-1));
        obj.oECI.n      =   1;

        if obj.sEF,         obj.setOutputExternalField('oECI','init');     end

        % Energy components
        obj.oECI.total              =   NaN(1,numel(obj.oECI.time));
        obj.oECI.Hamiltonian        =   NaN(1,numel(obj.oECI.time));
        obj.oECI.externalField      =   NaN(1,numel(obj.oECI.time));
        obj.oECI.autonomization     =   NaN(1,numel(obj.oECI.time));
    else
        obj.oECI.ind    =   NaN;
        obj.oECI.n      =   1;
    end
    
    % Wave function energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sEWfcn
        % When to save the DFT energies
        obj.oEWfcn.ind  =   [obj.getOutputIndex(obj.sEWfcnT) NaN];
        obj.oEWfcn.time =   obj.tspan(obj.oEWfcn.ind(1:end-1));
        obj.oEWfcn.n    =   1;

        if obj.sEF,         obj.setOutputExternalField('oEWfcn','init');     end
        
        % Energy components
        obj.oEWfcn.waveFunction         =   NaN(size(obj.CI.wfcn,2),numel(obj.oEWfcn.time));
    else
        obj.oEWfcn.ind  =   NaN;
        obj.oEWfcn.n    =   1;  
    end

case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main components
    if obj.sECI,    obj.oECI    =   rmfield(obj.oECI,{'ind','n'});          if obj.sEF %#ok<ALIGN> 
                    obj.setOutputExternalField('oECI','clean');             end
    else,           obj.oECI    =   [];                                     end
    if obj.sEWfcn,  obj.oEWfcn  =   rmfield(obj.oEWfcn,{'ind','n'});        if obj.sEF %#ok<ALIGN> 
                    obj.setOutputExternalField('oEWfcn','clean');           end
    else,           obj.oEWfcn  =   [];                                     end

otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unexpected option
    error('QMol:TDCI:setOutputEnergy',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end