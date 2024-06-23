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

        % Total dipole
        obj.oDip.total  =   NaN(obj.DFT.dim,numel(obj.oDip.time));           if obj.DFT.isSpinPol
        obj.oDip.totalUp=   NaN(obj.DFT.dim,numel(obj.oDip.time));
        obj.oDip.totalDown= NaN(obj.DFT.dim,numel(obj.oDip.time));           end

        % Orbital resolved dipole
        if ischar(obj.sDipI), if strcmpi(obj.sDipI,'all')                   % Save all orbitals
            if obj.DFT.isSpinPol
                % What to save
                obj.oDip.indexOrbitalUp=   1:numel(obj.DFT.occ{1});
                obj.oDip.indexOrbitalDown= 1:numel(obj.DFT.occ{2});

                % Spin-up allocation
                obj.oDip.orbitalUp_x    =   NaN(numel(obj.DFT.occ{1}),numel(obj.oDip.time)); if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbitalUp_y    =   NaN(numel(obj.DFT.occ{1}),numel(obj.oDip.time)); if obj.DFT.dim > 2
                obj.oDip.orbitalUp_z    =   NaN(numel(obj.DFT.occ{1}),numel(obj.oDip.time)); end, end

                % Spin down allocation
                obj.oDip.orbitalDown_x  =   NaN(numel(obj.DFT.occ{2}),numel(obj.oDip.time)); if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbitalDown_y  =   NaN(numel(obj.DFT.occ{2}),numel(obj.oDip.time)); if obj.DFT.dim > 2
                obj.oDip.orbitalDown_z  =   NaN(numel(obj.DFT.occ{2}),numel(obj.oDip.time)); end, end
            else
                % What to save
                obj.oDip.indexOrbital   =   1:numel(obj.DFT.occ);
                
                % Allocation
                obj.oDip.orbital_x      =   NaN(numel(obj.DFT.occ),numel(obj.oDip.time)); if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbital_y      =   NaN(numel(obj.DFT.occ),numel(obj.oDip.time)); if obj.DFT.dim > 2
                obj.oDip.orbital_z      =   NaN(numel(obj.DFT.occ),numel(obj.oDip.time)); end, end
            end,              else
                warning('QMol:TDDFT:outputDipole',...
                    ['Unknown option ' obj.sDipI ' for orbital-resolved dipole; feature disabled']);
                if obj.DFT.isSpinPol,   obj.oDip.indexOrbitalUp   =   [];     obj.oDip.indexOrbitalDown   =   [];
                else,                   obj.oDip.indexOrbital     =   [];
                end,          end
        elseif ~isempty(obj.sDipI)                                          % User-supplied indexes
            if obj.DFT.isSpinPol,       obj.oDip.indexOrbitalUp   =   obj.sDipI{1}(:).';
                                        obj.oDip.indexOrbitalDown =   obj.sDipI{2}(:).';
                obj.oDip.orbitalUp_x    =   NaN(numel(obj.sDipI{1}),numel(obj.oDip.time));   if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbitalUp_y    =   NaN(numel(obj.sDipI{1}),numel(obj.oDip.time));   if obj.DFT.dim > 2
                obj.oDip.orbitalUp_z    =   NaN(numel(obj.sDipI{1}),numel(obj.oDip.time));   end, end
                obj.oDip.orbitalDown_x  =   NaN(numel(obj.sDipI{2}),numel(obj.oDip.time));   if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbitalDown_y  =   NaN(numel(obj.sDipI{2}),numel(obj.oDip.time));   if obj.DFT.dim > 2
                obj.oDip.orbitalDown_z  =   NaN(numel(obj.sDipI{2}),numel(obj.oDip.time));   end, end
            else,                       obj.oDip.indexOrbital     =   obj.sDipI(:).';
                obj.oDip.orbital_x      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbital_y      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.DFT.dim > 2
                obj.oDip.orbital_z      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  end, end
            end
        else                                                                % No orbital resolving
            if obj.DFT.isSpinPol,       obj.oDip.indexOrbitalUp   =   [];     obj.oDip.indexOrbitalDown   =   [];
            else,                       obj.oDip.indexOrbital     =   [];
            end
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

        % Total dipole velocity
        obj.oVel.total  =   NaN(obj.DFT.dim,numel(obj.oVel.time));           if obj.DFT.isSpinPol
        obj.oVel.totalUp=   NaN(obj.DFT.dim,numel(obj.oVel.time));
        obj.oVel.totalDown= NaN(obj.DFT.dim,numel(obj.oVel.time));           end

        % Orbital resolved dipole velocity
        if ischar(obj.sVelI), if strcmpi(obj.sVelI,'all')                   % Save all orbitals
            if obj.DFT.isSpinPol
                % What to save
                obj.oVel.indexOrbitalUp =   1:numel(obj.DFT.occ{1});
                obj.oVel.indexOrbitalDown=  1:numel(obj.DFT.occ{2});

                % Spin-up allocation
                obj.oVel.orbitalUp_x    =   NaN(numel(obj.DFT.occ{1}),numel(obj.oVel.time)); if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oVel.orbitalUp_y    =   NaN(numel(obj.DFT.occ{1}),numel(obj.oVel.time)); if obj.DFT.dim > 2
                obj.oVel.orbitalUp_z    =   NaN(numel(obj.DFT.occ{1}),numel(obj.oVel.time)); end, end

                % Spin down allocation
                obj.oVel.orbitalDown_x  =   NaN(numel(obj.DFT.occ{2}),numel(obj.oVel.time)); if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oVel.orbitalDown_y  =   NaN(numel(obj.DFT.occ{2}),numel(obj.oVel.time)); if obj.DFT.dim > 2
                obj.oVel.orbitalDown_z  =   NaN(numel(obj.DFT.occ{2}),numel(obj.oVel.time)); end, end
            else
                % What to save
                obj.oVel.indexOrbital   =   1:numel(obj.DFT.occ);
                
                % Allocation
                obj.oVel.orbital_x      =   NaN(numel(obj.DFT.occ),numel(obj.oVel.time)); if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oVel.orbital_y      =   NaN(numel(obj.DFT.occ),numel(obj.oVel.time)); if obj.DFT.dim > 2
                obj.oVel.orbital_z      =   NaN(numel(obj.DFT.occ),numel(obj.oVel.time)); end, end
            end,              else
                warning('QMol:TDDFT:outputDipoleVelocity',...
                    ['Unknown option ' obj.sVelI ' for orbital-resolved dipole velocity; feature disabled']);
                if obj.DFT.isSpinPol,   obj.oVel.indexOrbitalUp   =   [];     obj.oVel.indexOrbitalDown   =   [];
                else,                   obj.oVel.indexOrbital     =   [];
                end,          end
        elseif ~isempty(obj.sVelI)                                          % User-supplied indexes
            if obj.DFT.isSpinPol,       obj.oVel.indexOrbitalUp   =   obj.sVelI{1}(:).';
                                        obj.oVel.indexOrbitalDown =   obj.sVelI{2}(:).';
                obj.oDip.orbitalUp_x    =   NaN(numel(obj.sDipI{1}),numel(obj.oDip.time));   if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbitalUp_y    =   NaN(numel(obj.sDipI{1}),numel(obj.oDip.time));   if obj.DFT.dim > 2
                obj.oDip.orbitalUp_z    =   NaN(numel(obj.sDipI{1}),numel(obj.oDip.time));   end, end
                obj.oDip.orbitalDown_x  =   NaN(numel(obj.sDipI{2}),numel(obj.oDip.time));   if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbitalDown_y  =   NaN(numel(obj.sDipI{2}),numel(obj.oDip.time));   if obj.DFT.dim > 2
                obj.oDip.orbitalDown_z  =   NaN(numel(obj.sDipI{2}),numel(obj.oDip.time));   end, end
            else,                       obj.oVel.indexOrbital     =   obj.sVelI(:).';
                obj.oDip.orbital_x      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbital_y      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.DFT.dim > 2
                obj.oDip.orbital_z      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  end, end
            end
        else                                                                % No orbital resolving
            if obj.DFT.isSpinPol,       obj.oVel.indexOrbitalUp   =   [];     obj.oVel.indexOrbitalDown   =   [];
            else,                       obj.oVel.indexOrbital     =   [];
            end
        end
    else
        obj.oVel.ind    =   NaN;
        obj.oVel.n      =   1;
    end
    
    % Dipole acceleration signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sAcc
        % When to save the dipole acceleration
        if ischar(obj.sAccT) && any(strcmpi(obj.sAccT,{'dipole','dip'}))    %#ok<ALIGN> 
                    obj.oAcc.ind    =   [obj.getOutputIndex(obj.sDipT),NaN];% Same as for dipole
        else,       obj.oAcc.ind    =   [obj.getOutputIndex(obj.sAccT),NaN];end % Own set of times
        obj.oAcc.time   =   obj.tspan(obj.oAcc.ind(1:end-1));
        obj.oAcc.n      =   1;

        if obj.sEF,         obj.setOutputExternalField('oAcc','init');      end

        % Total dipole acceleration
        obj.oAcc.total  =   NaN(obj.DFT.dim,numel(obj.oAcc.time));           if obj.DFT.isSpinPol
        obj.oAcc.totalUp=   NaN(obj.DFT.dim,numel(obj.oAcc.time));
        obj.oAcc.totalDown= NaN(obj.DFT.dim,numel(obj.oAcc.time));           end

        % Orbital resolved dipole acceleration
        if ischar(obj.sAccI), if strcmpi(obj.sAccI,'all')                   % Save all orbitals
            if obj.DFT.isSpinPol
                % What to save
                obj.oAcc.indexOrbitalUp =   1:numel(obj.DFT.occ{1});
                obj.oAcc.indexOrbitalDown=  1:numel(obj.DFT.occ{2});

                % Spin-up allocation
                obj.oAcc.orbitalUp_x    =   NaN(numel(obj.DFT.occ{1}),numel(obj.oAcc.time)); if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oAcc.orbitalUp_y    =   NaN(numel(obj.DFT.occ{1}),numel(obj.oAcc.time)); if obj.DFT.dim > 2
                obj.oAcc.orbitalUp_z    =   NaN(numel(obj.DFT.occ{1}),numel(obj.oAcc.time)); end, end

                % Spin down allocation
                obj.oAcc.orbitalDown_x  =   NaN(numel(obj.DFT.occ{2}),numel(obj.oAcc.time)); if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oAcc.orbitalDown_y  =   NaN(numel(obj.DFT.occ{2}),numel(obj.oAcc.time)); if obj.DFT.dim > 2
                obj.oAcc.orbitalDown_z  =   NaN(numel(obj.DFT.occ{2}),numel(obj.oAcc.time)); end, end
            else
                % What to save
                obj.oAcc.indexOrbital   =   1:numel(obj.DFT.occ);
                
                % Allocation
                obj.oAcc.orbital_x      =   NaN(numel(obj.DFT.occ),numel(obj.oAcc.time)); if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oAcc.orbital_y      =   NaN(numel(obj.DFT.occ),numel(obj.oAcc.time)); if obj.DFT.dim > 2
                obj.oAcc.orbital_z      =   NaN(numel(obj.DFT.occ),numel(obj.oAcc.time)); end, end
            end,              else
                warning('QMol:TDDFT:outputDipoleAcceleration',...
                    ['Unknown option ' obj.sAccI ' for orbital-resolved dipole acceleration; feature disabled']);
                if obj.DFT.isSpinPol,   obj.oAcc.indexOrbitalUp   =   [];     obj.oAcc.indexOrbitalDown   =   [];
                else,                   obj.oAcc.indexOrbital     =   [];
                end,          end
        elseif ~isempty(obj.sAccI)                                          % User-supplied indexes
            if obj.DFT.isSpinPol,       obj.oAcc.indexOrbitalUp   =   obj.sAccI{1}(:).';
                                        obj.oAcc.indexOrbitalDown =   obj.sAccI{2}(:).';
                obj.oDip.orbitalUp_x    =   NaN(numel(obj.sDipI{1}),numel(obj.oDip.time));   if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbitalUp_y    =   NaN(numel(obj.sDipI{1}),numel(obj.oDip.time));   if obj.DFT.dim > 2
                obj.oDip.orbitalUp_z    =   NaN(numel(obj.sDipI{1}),numel(obj.oDip.time));   end, end
                obj.oDip.orbitalDown_x  =   NaN(numel(obj.sDipI{2}),numel(obj.oDip.time));   if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbitalDown_y  =   NaN(numel(obj.sDipI{2}),numel(obj.oDip.time));   if obj.DFT.dim > 2
                obj.oDip.orbitalDown_z  =   NaN(numel(obj.sDipI{2}),numel(obj.oDip.time));   end, end
            else,                       obj.oAcc.indexOrbital     =   obj.sAccI(:).';
                obj.oDip.orbital_x      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.DFT.dim > 1 %#ok<ALIGN> 
                obj.oDip.orbital_y      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  if obj.DFT.dim > 2
                obj.oDip.orbital_z      =   NaN(numel(obj.sDipI),numel(obj.oDip.time));  end, end
            end
        else                                                                % No orbital resolving
            if obj.DFT.isSpinPol,       obj.oAcc.indexOrbitalUp   =   [];     obj.oAcc.indexOrbitalDown   =   [];
            else,                       obj.oAcc.indexOrbital     =   [];
            end
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
    error('QMol:TDDFT:setOutputDipole',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end