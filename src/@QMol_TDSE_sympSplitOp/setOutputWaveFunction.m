function setOutputWaveFunction(obj,opt)
%setOutputOrbitalDensity
 
    switch lower(opt)
    case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%
     
        % Clean up any old data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        obj.oWfcn           =   [];
        obj.oWfcnP          =   [];
        
        % Wave functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if obj.sWfcn
            % When to save the wave functions
            obj.oWfcn.ind   =   [obj.getOutputIndex(obj.sWfcnT),NaN];
            obj.oWfcn.time  =   obj.tspan(obj.oWfcn.ind(1:end-1));
            obj.oWfcn.n     =   1;

            if obj.sEF,         obj.setOutputExternalField('oWfcn','init');  end
            
            S               =   obj.SE.disc.SE_sizeWaveFunction;
            if numel(S) > 1   &&   S(end) == 1,     S   =   S(1:end-1);     end
            
            % Indexes
            if isempty(obj.sWfcnI),         obj.oWfcn.indexWaveFunction =   [];
            elseif ischar(obj.sWfcnI),                                                                  if strcmpi(obj.sWfcnI,'all')
                                            obj.oWfcn.indexWaveFunction =   1:obj.SE.N;                 else
                                            obj.oWfcn.indexWaveFunction =   [];
                                            warning('QMol:TDSE:saveWfcn', ...
                                                ['Unknown option for the saving of the wave functions ' obj.sKSOI ...
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
        
        % Wave functions' projection ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if obj.sWfcnP
            % When to save the wave function projections
            obj.oWfcnP.ind  =   [obj.getOutputIndex(obj.sWfcnPT),NaN];
            obj.oWfcnP.time =   obj.tspan(obj.oWfcnP.ind(1:end-1));
            obj.oWfcnP.n    =   1;

            if obj.sEF,         obj.setOutputExternalField('oWfcnP','init');    end
            
            % Projection basis
            if isempty(obj.sWfcnB)                                      % Project on initial wave functions
                obj.oWfcnP.basis    =   QMol_disc_basis('x',obj.SE.disc.x,'v',obj.SE.wfcn.wfcn);
                obj.oWfcnP.basis.initialize;
            elseif isnumeric(obj.sWfcnB)
                obj.oWfcnP.basis    =   QMol_disc_basis('x',obj.SE.disc.x,'v',obj.sWfcnB);
                obj.oWfcnP.basis.initialize;
            else
                obj.oWfcnP.basis    =   obj.sWfcnB;
                obj.oWfcnP.basis.initialize;
            end

            % Indexes 
            S               =   obj.oWfcnP.basis.SE_sizeWaveFunction;
            if numel(S) > 1   &&   S(end) == 1,     S   =   S(1:end-1);   end

            if isempty(obj.sWfcnPI),        obj.oWfcnP.indexWaveFunction    =   [];
            elseif ischar(obj.sWfcnPI),                                                     if strcmpi(obj.sWfcnPI,'all')
                                            obj.oWfcnP.indexWaveFunction    =   1:obj.SE.N; else
                                            obj.oWfcnP.indexWaveFunction    =   [];
                                            warning('QMol:TDSE:saveWfcnProjection', ...
                                                ['Unknown option for the saving of the wave functions'' projection ' obj.sWfcnPI ...
                                                '. No wave function projections saved.']);  end
            else,                           obj.oWfcnP.indexWaveFunction    =   obj.sWfcnPI;
            end

                % Allocate output
                obj.oWfcnP.shape        =   [S numel(obj.oWfcnP.indexWaveFunction)];
                obj.oWfcnP.waveFunction =   NaN([obj.oWfcnP.shape numel(obj.oWfcnP.time)],'like',1i);

        else
            obj.oWfcnP.shape    =   [];
            obj.oWfcnP.ind      =   NaN;
            obj.oWfcnP.n        =   1;
        end
        
    case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if obj.sWfcn,   obj.oWfcn   =   rmfield(obj.oWfcn,{'ind','n','shape'}); if obj.sEF %#ok<ALIGN> 
                        obj.setOutputExternalField('oWfcn' ,'clean');           end
        else,           obj.oWfcn   =   [];                                     end
        if obj.sWfcnP,  obj.oWfcnP  =   rmfield(obj.oWfcnP,{'ind','n','shape'});if obj.sEF %#ok<ALIGN> 
                        obj.setOutputExternalField('oWfcnP','clean');           end
        else,           obj.oWfcnP   =   [];                                    end
     
    otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Unexpected option
        error('QMol:TDSE:setOutputWaveFunction',['Unknown option ' opt]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end