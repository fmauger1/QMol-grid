function setOutputEnergy(obj,opt)
%setOutputEnergy

switch lower(opt)
case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Clean up any old data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.oEDFT           =   [];
    obj.oEKSO           =   [];
    
    % DFT energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sEDFT
        % When to save the DFT energies
        obj.oEDFT.ind   =   [obj.getOutputIndex(obj.sEDFTT),NaN];
        obj.oEDFT.time  =   obj.tspan(obj.oEDFT.ind(1:end-1));
        obj.oEDFT.n     =   1;

        if obj.sEF,         obj.setOutputExternalField('oEDFT','init');     end

        % Energy components
        obj.oEDFT.total                 =   NaN(1,numel(obj.oEDFT.time));    if obj.DFT.isSpinPol
        obj.oEDFT.kinetic               =   NaN(2,numel(obj.oEDFT.time));
        obj.oEDFT.external              =   NaN(2,numel(obj.oEDFT.time));
        obj.oEDFT.externalField         =   NaN(2,numel(obj.oEDFT.time));    else
        obj.oEDFT.kinetic               =   NaN(1,numel(obj.oEDFT.time));
        obj.oEDFT.external              =   NaN(1,numel(obj.oEDFT.time));
        obj.oEDFT.externalField         =   NaN(1,numel(obj.oEDFT.time));    end
        obj.oEDFT.Hartree               =   NaN(1,numel(obj.oEDFT.time));
        obj.oEDFT.exchangeCorrelation   =   NaN(1,numel(obj.oEDFT.time));
        obj.oEDFT.autonomization        =   NaN(1,numel(obj.oEDFT.time));
    else
        obj.oEDFT.ind   =   NaN;
        obj.oEDFT.n     =   1;
    end
    
    % Orbital energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sEKSO
        % When to save the DFT energies
        obj.oEKSO.ind   =   [obj.getOutputIndex(obj.sEKSOT) NaN];
        obj.oEKSO.time  =   obj.tspan(obj.oEKSO.ind(1:end-1));
        obj.oEKSO.n     =   1;

        if obj.sEF,         obj.setOutputExternalField('oEKSO','init');     end,                 if obj.DFT.isSpinPol
        
        % Energy components
        obj.oEKSO.orbitalUp             =   NaN(numel(obj.DFT.occ{1}),numel(obj.oEKSO.time));
        obj.oEKSO.orbitalDown           =   NaN(numel(obj.DFT.occ{2}),numel(obj.oEKSO.time));    else
        obj.oEKSO.orbital               =   NaN(numel(obj.DFT.occ   ),numel(obj.oEKSO.time));    end
    else
        obj.oEKSO.ind   =   NaN;
        obj.oEKSO.n     =   1;  
    end

case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main components
    if obj.sEDFT,   obj.oEDFT   =   rmfield(obj.oEDFT,{'ind','n'});         if obj.sEF %#ok<ALIGN> 
                    obj.setOutputExternalField('oEDFT','clean');            end
    else,           obj.oEDFT   =   [];                                     end
    if obj.sEKSO,   obj.oEKSO   =   rmfield(obj.oEKSO,{'ind','n'});         if obj.sEF %#ok<ALIGN> 
                    obj.setOutputExternalField('oEKSO','clean');            end
    else,           obj.oEKSO   =   [];                                     end

otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unexpected option
    error('QMol:TDDFT:setOutputEnergy',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end