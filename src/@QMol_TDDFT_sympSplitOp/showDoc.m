function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the TDDFT solver

    % Test for LDA-type functionals
    if iscell(obj.DFT.Vxc),     for k = 1:numel(obj.DFT.Vxc),       if ~startsWith(obj.DFT.Vxc{k}.type,{'LDA','GGA'})   %#ok<ALIGN> 
        fprintf('\n  WARNING: ***************************************************************\n');
        fprintf('    %s-type DFT functional detected; propagation may not be accurate\n',obj.DFT.Vxc{k}.type);
        fprintf('    **********************************************************************\n\n'); end, end
    elseif ~isempty(obj.DFT.Vxc),                                   if ~startsWith(obj.DFT.Vxc.type,{'LDA','GGA'})      %#ok<ALIGN> 
        fprintf('\n  WARNING: ***************************************************************\n');
        fprintf('    %s-type DFT functional detected; propagation may not be accurate\n',obj.DFT.Vxc.type);
        fprintf('    **********************************************************************\n\n'); end
    end

    % Implementation of the TDDFT model
    fprintf('  * Symplectic split-operator propagation scheme [Mauger 2024],\n');
    if isempty(obj.splitMotif),             fprintf('    with a VTV split motif.\n');
    elseif strcmpi(obj.splitMotif,'vtv'),   fprintf('    with a VTV split motif.\n');
    elseif strcmpi(obj.splitMotif,'tvt'),   fprintf('    with a TVT split motif.\n');
    else,   warning('QMol:TDDFT:splitMotif',['Unknown motif for symplectic split ' obj.splitMotif])
    end
    obj.version;

    % External field
    if (isempty(obj.EF) || (isempty(obj.EF.electricField)   &&   isempty(obj.EF.potentialVector))) ...
                        && isempty(obj.FE)   &&   isempty(obj.FA)                           %#ok<ALIGN> 
                                                                    FG  =   0;              % No usable field
    elseif isempty(obj.EFG)                                             % Default mode (context based)
        if isempty(obj.EF)                                                      %#ok<ALIGN>
            if isa(obj.FE,'function_handle'),                       FG  =   1;  %#ok<ALIGN> % Length gauge
            elseif isa(obj.FA,'function_handle'),                   FG  =   2;              % Velocity gauge
            elseif ~isempty(obj.FE),                                FG  =   1;              % Length gauge
            elseif ~isempty(obj.FA),                                FG  =   2;              % Velocity gauge
            else,                                                   FG  =   0;  end         % No field to work with
        elseif isa(obj.EF.electricField,'function_handle'),         FG  =   1;              % Length gauge
        elseif isa(obj.EF.potentialVector,'function_handle'),       FG  =   2;              % Velocity gauge
        elseif ~isempty(obj.EF.electricField),                      FG  =   1;              % Length gauge
        elseif ~isempty(obj.EF.potentialVector),                    FG  =   2;              % Velocity gauge
        else,                                                       FG  =   0;  end         % No field to work with
    else, switch lower(obj.EFG)                         %#ok<ALIGN> 
        case {'length','length gauge','length_gauge','e'},          FG  =   1;
        case {'velocity','velocity gauge','velocity_gauge','a'},    FG  =   2;
        case {'none','off'},                                        FG  =   0;
        otherwise
            warning('QMol:TDDFT:gaugeField', ['Unknown gauge-field option ' obj.Gfield '; Default gauge mode used instead.'])

            if isempty(obj.EF),                                     FG  =   0;  %#ok<ALIGN> % No field to work with
            elseif isa(obj.EF.electricField,'function_handle'),     FG  =   1;              % Length gauge
            elseif isa(obj.EF.potentialVector,'function_handle'),   FG  =   2;              % Velocity gauge
            elseif ~isempty(obj.EF.electricField),                  FG  =   1;              % Length gauge
            elseif ~isempty(obj.EF.potentialVector),                FG  =   2;              % Velocity gauge
            else,                                                   FG  =   0;  end         % No field to work with
    end, end

    if FG == 0,     if ~isempty(obj.EF),    fprintf('  * WARNING: The external field is ignored in the propagation. ***********\n'); %#ok<ALIGN> 
                    else,                   fprintf('  * Field free propagation\n');    end
    elseif FG == 1,                         fprintf('  * Field driven propagation\n    dipole approximation, length gauge\n');
    elseif FG == 2,                         fprintf('  * Field driven propagation\n    dipole approximation, velocity gauge\n');
    else,   error('QMol:TDDFT:sympSplitOpGauge','Unexpected gauge mode');               end

    % References
    ref                 =   {'Mauger 2024'};

end