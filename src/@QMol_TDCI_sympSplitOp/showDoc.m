function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the TDCI solver

    % Implementation of the TDCI model
    fprintf('  * Symplectic split-operator propagation scheme [Mauger 2025],\n');
    if isempty(obj.splitMotif),             fprintf('    with a HDH split motif.\n');
    elseif strcmpi(obj.splitMotif,'hdh'),   fprintf('    with a HDH split motif.\n');
    elseif strcmpi(obj.splitMotif,'dhd'),   fprintf('    with a DHD split motif.\n');
    else,   warning('QMol:TDCI:splitMotif',['Unknown motif for symplectic split ' obj.splitMotif])
    end
    obj.version;

    % External field
    if isempty(obj.CI.DX),                                          FG  =   0; %#ok<ALIGN> % No dipole-coupling matrix
    elseif (isempty(obj.EF) || (isempty(obj.EF.electricField)   &&   isempty(obj.EF.potentialVector))) ...
                        && isempty(obj.FE)   &&   isempty(obj.FA)
                                                                    FG  =   0;              % No usable field
    elseif isempty(obj.EFG)                                             % Default mode (context based)
        if isempty(obj.EF)                                                      %#ok<ALIGN>
            if isa(obj.FE,'function_handle'),                       FG  =   1;  %#ok<ALIGN> % Length gauge
            elseif isa(obj.FA,'function_handle'),                   FG  =   1;              % Velocity gauge
            elseif ~isempty(obj.FE),                                FG  =   1;              % Length gauge
            elseif ~isempty(obj.FA),                                FG  =   1;              % Velocity gauge
            else,                                                   FG  =   0;  end         % No field to work with
        elseif isa(obj.EF.electricField,'function_handle'),         FG  =   1;              % Length gauge
        elseif isa(obj.EF.potentialVector,'function_handle'),       FG  =   1;              % Velocity gauge
        elseif ~isempty(obj.EF.electricField),                      FG  =   1;              % Length gauge
        elseif ~isempty(obj.EF.potentialVector),                    FG  =   1;              % Velocity gauge
        else,                                                       FG  =   0;  end         % No field to work with
    else, switch lower(obj.EFG)                         %#ok<ALIGN> 
        case {'length','length gauge','length_gauge','e'},          FG  =   1;
        case {'none','off'},                                        FG  =   0;
        otherwise
            warning('QMol:TDCI:gaugeField', ['Unknown gauge-field option ' obj.EFG '; Default gauge mode used instead.'])

            if isempty(obj.EF),                                     FG  =   0;  %#ok<ALIGN> % No field to work with
            elseif isa(obj.EF.electricField,'function_handle'),     FG  =   1;              % Length gauge
            elseif isa(obj.EF.potentialVector,'function_handle'),   FG  =   1;              % Velocity gauge
            elseif ~isempty(obj.EF.electricField),                  FG  =   1;              % Length gauge
            elseif ~isempty(obj.EF.potentialVector),                FG  =   1;              % Velocity gauge
            else,                                                   FG  =   0;  end         % No field to work with
    end, end

    if FG == 0,     if ~isempty(obj.EF),    fprintf('  * WARNING: The external field is ignored in the propagation. ***********\n'); %#ok<ALIGN> 
                    else,                   fprintf('  * Field free propagation\n');    end
    elseif FG == 1,                         fprintf('  * Field driven propagation\n    dipole approximation, length gauge\n');
    else,   error('QMol:TDCI:sympSplitOpGauge','Unexpected gauge mode');               end

    % References
    ref                 =   {'Mauger 2025'};

end