function setOutputWaveFunction(obj,opt)
%setOutputOrbitalDensity
 
    switch lower(opt)
    case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%
     
        % Clean up any old data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        obj.oWfcn           =   [];
        
        % Wave functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if obj.sWfcn
            % When to save the wave functions
            obj.oWfcn.ind   =   [obj.getOutputIndex(obj.sWfcnT),NaN];
            obj.oWfcn.time  =   obj.tspan(obj.oWfcn.ind(1:end-1));
            obj.oWfcn.n     =   1;

            if obj.sEF,         obj.setOutputExternalField('oWfcn','init');  end
            
            S               =   size(obj.CI.CSB,1);
            if numel(S) > 1   &&   S(end) == 1,     S   =   S(1:end-1);     end
            
            % Indexes
            if isempty(obj.sWfcnI),         obj.oWfcn.indexWaveFunction =   [];
            elseif ischar(obj.sWfcnI),                                                                  if strcmpi(obj.sWfcnI,'all')
                                            obj.oWfcn.indexWaveFunction =   1:size(obj.CI.wfcn,2);      else
                                            obj.oWfcn.indexWaveFunction =   [];
                                            warning('QMol:TDCI:saveWfcn', ...
                                                ['Unknown option for the saving of the wave functions ' obj.sWfcnI ...
                                                '. No orbitals saved.']);                               end
            else,                           obj.oWfcn.indexWaveFunction =   obj.sWfcnI;
            end

            % Allocate output
            obj.oWfcn.shape         =   [S numel(obj.oWfcn.indexWaveFunction)];
            obj.oWfcn.waveFunction  =   NaN([obj.oWfcn.shape numel(obj.oWfcn.time)],'like',1i);

        else
            obj.oWfcn.shape =   [];
            obj.oWfcn.ind   =   NaN;
            obj.oWfcn.n     =   1;
        end
        
    case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if obj.sWfcn,   obj.oWfcn   =   rmfield(obj.oWfcn,{'ind','n','shape'}); if obj.sEF %#ok<ALIGN> 
                        obj.setOutputExternalField('oWfcn' ,'clean');           end
        else,           obj.oWfcn   =   [];                                     end
     
    otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Unexpected option
        error('QMol:TDCI:setOutputWaveFunction',['Unknown option ' opt]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end