function setOutputDipole(obj,opt)
%setOutputDipole

switch lower(opt)
case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Clean up any old data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.oDip            =   [];
    obj.oVel            =   [];
    obj.oAcc            =   [];
    
    % Dipole signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sDip
        % When to save the dipole
        obj.oDip.ind    =   [obj.getOutputIndex(obj.sDipT),NaN];
        obj.oDip.time   =   obj.tspan(obj.oDip.ind(1:end-1));
        obj.oDip.n      =   1;

        if obj.sEF,         obj.setOutputExternalField('oDip','init');      end

        % Dipole
        if ischar(obj.sDipI),                       if strcmpi(obj.sDipI,'all')     % Save all wave functions
            % What to save
            obj.oDip.indexWaveFunction   =   1:obj.SE.N;
            
            % Allocation
            obj.oDip.waveFunction_x =   NaN(obj.SE.N,numel(obj.oDip.time)); if obj.SE.dim > 1 %#ok<ALIGN> 
            obj.oDip.waveFunction_y =   NaN(obj.SE.N,numel(obj.oDip.time)); if obj.SE.dim > 2
            obj.oDip.waveFunction_z =   NaN(obj.SE.N,numel(obj.oDip.time)); end, end
                                                    else
            warning('QMol:TDSE:outputDipole',...
                ['Unknown option ' obj.sDipI ' for dipole; feature disabled']);
            obj.oDip.indexWaveFunction  =   [];     end
        elseif ~isempty(obj.sDipI)                                                  % User-supplied indexes
            obj.oDip.indexWaveFunction  =   obj.sDipI(:).';
            obj.oDip.waveFunction_x      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.SE.dim > 1 %#ok<ALIGN> 
            obj.oDip.waveFunction_y      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.SE.dim > 2
            obj.oDip.waveFunction_z      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  end, end
        else                                                                        % Default is all wave functions
            % What to save
            obj.oDip.indexWaveFunction   =   1:obj.SE.N;
            
            % Allocation
            obj.oDip.waveFunction_x      =   NaN(obj.SE.N,numel(obj.oDip.time)); if obj.SE.dim > 1 %#ok<ALIGN> 
            obj.oDip.waveFunction_y      =   NaN(obj.SE.N,numel(obj.oDip.time)); if obj.SE.dim > 2
            obj.oDip.waveFunction_z      =   NaN(obj.SE.N,numel(obj.oDip.time)); end, end
        end
    else
        obj.oDip.ind    =   NaN;
        obj.oDip.n      =   1;
    end
    
    % Dipole velocity signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sVel
        % When to save the dipole velocity
        if ischar(obj.sVelT) && any(strcmpi(obj.sVelT,{'dipole','dip'}))    %#ok<ALIGN> 
                    obj.oVel.ind    =   [obj.getOutputIndex(obj.sDipT),NaN];% Same as for dipole
        else,       obj.oVel.ind    =   [obj.getOutputIndex(obj.sVelT),NaN];end % Own set of times
        obj.oVel.time   =   obj.tspan(obj.oVel.ind(1:end-1));
        obj.oVel.n      =   1;

        if obj.sEF,         obj.setOutputExternalField('oVel','init');      end

        % Dipole velocity
        if ischar(obj.sVelI),                       if strcmpi(obj.sVelI,'all')                   % Save all wave functions
            % What to save
            obj.oVel.indexWaveFunction   =   1:obj.SE.N;
            
            % Allocation
            obj.oVel.waveFunction_x      =   NaN(obj.SE.N,numel(obj.oVel.time)); if obj.SE.dim > 1 %#ok<ALIGN> 
            obj.oVel.waveFunction_y      =   NaN(obj.SE.N,numel(obj.oVel.time)); if obj.SE.dim > 2
            obj.oVel.waveFunction_z      =   NaN(obj.SE.N,numel(obj.oVel.time)); end, end
                                                    else
            warning('QMol:TDSE:outputDipoleVelocity',...
                ['Unknown option ' obj.sVelI ' for dipole velocity; feature disabled']);
            obj.oVel.indexWaveFunction   =   [];    end
        elseif ~isempty(obj.sVelI)                                          % User-supplied indexes
            obj.oVel.indexWaveFunction   =   obj.sVelI(:).';
            obj.oDip.waveFunction_x      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.SE.dim > 1 %#ok<ALIGN> 
            obj.oDip.waveFunction_y      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.SE.dim > 2
            obj.oDip.waveFunction_z      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  end, end
        else                                                                % Default is all wave functions
            % What to save
            obj.oVel.indexWaveFunction   =   1:obj.SE.N;
            
            % Allocation
            obj.oVel.waveFunction_x      =   NaN(obj.SE.N,numel(obj.oVel.time)); if obj.SE.dim > 1 %#ok<ALIGN> 
            obj.oVel.waveFunction_y      =   NaN(obj.SE.N,numel(obj.oVel.time)); if obj.SE.dim > 2
            obj.oVel.waveFunction_z      =   NaN(obj.SE.N,numel(obj.oVel.time)); end, end
        end
    else
        obj.oVel.ind    =   NaN;
        obj.oVel.n      =   1;
    end
    
    % Dipole acceleration signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sAcc
        % Initialize the potential derivative
        obj.SE.V.setDerivative;

        % When to save the dipole acceleration
        if ischar(obj.sAccT) && any(strcmpi(obj.sAccT,{'dipole','dip'}))    %#ok<ALIGN> 
                    obj.oAcc.ind    =   [obj.getOutputIndex(obj.sDipT),NaN];% Same as for dipole
        else,       obj.oAcc.ind    =   [obj.getOutputIndex(obj.sAccT),NaN];end % Own set of times
        obj.oAcc.time   =   obj.tspan(obj.oAcc.ind(1:end-1));
        obj.oAcc.n      =   1;

        if obj.sEF,         obj.setOutputExternalField('oAcc','init');      end

        % Dipole acceleration
        if ischar(obj.sAccI),                       if strcmpi(obj.sAccI,'all')                   % Save all orbitals
            % What to save
            obj.oAcc.indexWaveFunction   =   1:obj.SE.N;
            
            % Allocation
            obj.oAcc.waveFunction_x      =   NaN(obj.SE.N,numel(obj.oAcc.time)); if obj.SE.dim > 1 %#ok<ALIGN> 
            obj.oAcc.waveFunction_y      =   NaN(obj.SE.N,numel(obj.oAcc.time)); if obj.SE.dim > 2
            obj.oAcc.waveFunction_z      =   NaN(obj.SE.N,numel(obj.oAcc.time)); end, end
                                                    else
            warning('QMol:TDSE:outputDipoleAcceleration',...
                ['Unknown option ' obj.sAccI ' for dipole acceleration; feature disabled']);
            obj.oAcc.indexWaveFunction   =   [];    end
        elseif ~isempty(obj.sAccI)                                          % User-supplied indexes
            obj.oAcc.indexWaveFunction   =   obj.sAccI(:).';
            obj.oDip.waveFunction_x      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.SE.dim > 1 %#ok<ALIGN> 
            obj.oDip.waveFunction_y      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.SE.dim > 2
            obj.oDip.waveFunction_z      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  end, end
        else                                                                % Default is all wave functions
            % What to save
            obj.oAcc.indexWaveFunction   =   1:obj.SE.N;
            
            % Allocation
            obj.oAcc.waveFunction_x      =   NaN(obj.SE.N,numel(obj.oAcc.time)); if obj.SE.dim > 1 %#ok<ALIGN> 
            obj.oAcc.waveFunction_y      =   NaN(obj.SE.N,numel(obj.oAcc.time)); if obj.SE.dim > 2
            obj.oAcc.waveFunction_z      =   NaN(obj.SE.N,numel(obj.oAcc.time)); end, end
        end
    else
        obj.oAcc.ind    =   NaN;
        obj.oAcc.n      =   1;
    end

case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main components
    if obj.sDip,    obj.oDip    =   rmfield(obj.oDip,{'ind','n'});          if obj.sEF %#ok<ALIGN> 
                    obj.setOutputExternalField('oDip','clean');             end
    else,           obj.oDip    =   [];                                     end
    if obj.sVel,    obj.oVel    =   rmfield(obj.oVel,{'ind','n'});          if obj.sEF %#ok<ALIGN> 
                    obj.setOutputExternalField('oVel','clean');             end
    else,           obj.oVel    =   [];                                     end
    if obj.sAcc,    obj.oAcc    =   rmfield(obj.oAcc,{'ind','n'});          if obj.sEF %#ok<ALIGN> 
                    obj.setOutputExternalField('oAcc','clean');             end
    else,           obj.oAcc    =   [];                                     end

otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unexpected option
    error('QMol:TDSE:setOutputDipole',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end