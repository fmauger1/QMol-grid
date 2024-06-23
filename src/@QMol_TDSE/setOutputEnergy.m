function setOutputEnergy(obj,opt)
%setOutputEnergy

switch lower(opt)
case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Clean up any old data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.oESE            =   [];
    obj.oEWfcn          =   [];
    
    % Schrodinger-equation energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sESE
        % When to save the DFT energies
        obj.oESE.ind    =   [obj.getOutputIndex(obj.sESET),NaN];
        obj.oESE.time   =   obj.tspan(obj.oESE.ind(1:end-1));
        obj.oESE.n      =   1;

        if obj.sEF,         obj.setOutputExternalField('oESE','init');     end

        % Energy components
        obj.oESE.total                 =   NaN(1,numel(obj.oESE.time));
        obj.oESE.kinetic               =   NaN(1,numel(obj.oESE.time));
        obj.oESE.potential             =   NaN(1,numel(obj.oESE.time));
        obj.oESE.externalField         =   NaN(1,numel(obj.oESE.time));
        obj.oESE.autonomization        =   NaN(1,numel(obj.oESE.time));
    else
        obj.oESE.ind    =   NaN;
        obj.oESE.n      =   1;
    end
    
    % Wave function energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sEWfcn
        % When to save the DFT energies
        obj.oEWfcn.ind  =   [obj.getOutputIndex(obj.sEWfcnT) NaN];
        obj.oEWfcn.time =   obj.tspan(obj.oEWfcn.ind(1:end-1));
        obj.oEWfcn.n    =   1;

        if obj.sEF,         obj.setOutputExternalField('oEWfcn','init');     end
        
        % Energy components
        obj.oEWfcn.waveFunction         =   NaN(obj.SE.N,numel(obj.oEWfcn.time));
    else
        obj.oEWfcn.ind  =   NaN;
        obj.oEWfcn.n    =   1;  
    end

case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main components
    if obj.sESE,    obj.oESE    =   rmfield(obj.oESE,{'ind','n'});          if obj.sEF %#ok<ALIGN> 
                    obj.setOutputExternalField('oESE','clean');             end
    else,           obj.oESE    =   [];                                     end
    if obj.sEWfcn,  obj.oEWfcn  =   rmfield(obj.oEWfcn,{'ind','n'});        if obj.sEF %#ok<ALIGN> 
                    obj.setOutputExternalField('oEWfcn','clean');           end
    else,           obj.oEWfcn  =   [];                                     end

otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unexpected option
    error('QMol:TDSE:setOutputEnergy',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end