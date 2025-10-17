function setOutputIonization(obj,opt)
%setOutputIonization

switch lower(opt)
case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Clean up any old data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.oIon            =   [];
    
    % Ionization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sIon
        % When to save the ionization
        obj.oIon.ind    =   [obj.getOutputIndex(obj.sIonT),NaN];
        obj.oIon.time   =   obj.tspan(obj.oIon.ind(1:end-1));
        obj.oIon.n      =   1;

        if obj.sEF,         obj.setOutputExternalField('oIon','init');      end
        
        % Ionization
        if ischar(obj.sIWfcnI),                     if strcmpi(obj.sIWfcnI,'all')       % Save all wave functions
            % What to save
            obj.oIon.indexWaveFunction  =   1:size(obj.CI.wfcn,2);
            
            % Allocation
            obj.oIon.waveFunction   =   NaN(size(obj.CI.wfcn,2),numel(obj.oIon.time));
                                                    else
            warning('QMol:TDCI:outputIonization',...
                ['Unknown option ' obj.sIonI ' for ionization; feature disabled']);
            obj.oIon.indexWaveFunction  = [];       end
        elseif ~isempty(obj.sIWfcnI)                                                    % User-supplied indexes
            obj.oIon.waveFunction   =   NaN(numel(obj.sIWfcnI),numel(obj.oIon.time));
        else                                                                            % Default is all wave functions
            % What to save
            obj.oIon.indexWaveFunction  =   1:size(obj.CI.wfcn,2);
            
            % Allocation
            obj.oIon.waveFunction   =   NaN(size(obj.CI.wfcn,2),numel(obj.oIon.time));
        end
    else
        obj.oIon.ind    =   NaN;
        obj.oIon.n      =   1;
    end
    
case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main components
    if obj.sIon,    obj.oIon    =   rmfield(obj.oIon,{'ind','n'});          if obj.sEF %#ok<ALIGN> 
                    obj.setOutputExternalField('oIon','clean');             end
    else,           obj.oIon    =   [];                                     end

otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unexpected option
    error('QMol:TDCI:setOutputIonization',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end