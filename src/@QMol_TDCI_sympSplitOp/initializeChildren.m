function initializeChildren(obj,isRst)
%initializeChildren

    % Common initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Miscellaneous
        obj.nWfcn       =   size(obj.CI.wfcn,2);

    if isRst
    % Restart ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Run-time variables
        obj.H0          =   obj.restart.H0;
        obj.D           =   obj.restart.D;

    else
    % From scratch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Integration motif
        if isempty(obj.splitMotif),             obj.isHDH   =   true;
        elseif strcmpi(obj.splitMotif,'hdh'),   obj.isHDH   =   true;
        elseif strcmpi(obj.splitMotif,'dhd'),   obj.isHDH   =   false;
        else,                                   obj.isHDH   =   true;
            warning('QMol:TDCI:splitMotif', ...
                ['Unknown motif for symplectic split ' obj.splitMotif '; default (HDH) used instead.'])
        end

        % CI+CAP and mask matrices
        obj.H0          =   -1i*obj.CI.CI;                                  if ~isempty(obj.ABC)
        if obj.ABC.isCAP,   obj.H0  =   obj.H0 - sign(obj.dt)*obj.ABC.getDampingMatrix;
        else,               obj.D   =   obj.ABC.getDampingMatrix;   end,    end

        % Remap external field
        if ~(isempty(obj.FA) && isempty(obj.FE) && isempty(obj.FDE)),                   if isempty(obj.EF)
            % Map to external field object
            obj.EF      =   QMol_extField_dipole('A',obj.FA,'E',obj.FE,'DE',obj.FDE);   if ~isempty(obj.FA)
            obj.EF.set('A', obj.FA);                                        end,        if ~isempty(obj.FE)
            obj.EF.set('E', obj.FE);                                        end,        if ~isempty(obj.FDE)
            obj.EF.set('DE',obj.FDE);                                       end,        end

            % Clear properties
            obj.FA      =   [];
            obj.FE      =   [];
            obj.FDE     =   [];
        end

        if ~isempty(obj.EF),        obj.EF.initialize(obj.CI.disc);        end

        % Choice of gauge
        if isempty(obj.CI.DX),                                          obj.FG  =   0; %#ok<ALIGN> % No dipole-coupling matrix
        elseif isempty(obj.EF) || (isempty(obj.EF.electricField)   &&   isempty(obj.EF.potentialVector))
                                                                        obj.FG  =   0;              % No usable field
        elseif isempty(obj.EFG)                                             % Default mode (context based)
            if isempty(obj.EF),                                         obj.FG  =   0;  %#ok<ALIGN> % No field to work with
            elseif isa(obj.EF.electricField,'function_handle'),         obj.FG  =   1;              % Length gauge
            elseif isa(obj.EF.potentialVector,'function_handle'),       obj.FG  =   1;              % Velocity gauge
            elseif ~isempty(obj.EF.electricField),                      obj.FG  =   1;              % Length gauge
            elseif ~isempty(obj.EF.potentialVector),                    obj.FG  =   1;              % Velocity gauge
            else,                                                       obj.FG  =   0;  end         % No field to work with
        else, switch lower(obj.EFG)                         %#ok<ALIGN> 
            case {'length','length gauge','length_gauge','e'},          obj.FG  =   1;
            case {'none','off'},                                        obj.FG  =   0;
            otherwise
                warning('QMol:TDDFT:gaugeField', ['Unknown gauge-field option ' obj.EFG '; Default gauge mode used instead.'])
    
                if isempty(obj.EF),                                     obj.FG  =   0;  %#ok<ALIGN> % No field to work with
                elseif isa(obj.EF.electricField,'function_handle'),     obj.FG  =   1;              % Length gauge
                elseif isa(obj.EF.potentialVector,'function_handle'),   obj.FG  =   1;              % Velocity gauge
                elseif ~isempty(obj.EF.electricField),                  obj.FG  =   1;              % Length gauge
                elseif ~isempty(obj.EF.potentialVector),                obj.FG  =   1;              % Velocity gauge
                else,                                                   obj.FG  =   0;  end         % No field to work with
        end, end
        
        % External fields
        if isempty(obj.EF)
            obj.uE      =   false;      obj.uA  =   false;      obj.uDE =   false;
        else
            obj.uE      =   isa(obj.EF.electricField,'function_handle');
            obj.uA      =   isa(obj.EF.potentialVector,'function_handle');
            obj.uDE     =   isa(obj.EF.electricFieldDerivative,'function_handle');
        end

        if obj.FG == 1,                                                     if ~obj.uE && ~obj.uA
            % What to use to compute fields (length gauge)
            if ~isempty(obj.EF.electricField),      obj.uE  =   true;   else,   obj.uA  =   true;   end
                                                                            end
            % Initialize relevant field
            obj.FE      =   obj.getFE( obj.tspan(1));                       if obj.sECI
            obj.FDE     =   obj.getFDE(obj.tspan(1));                       end

        end

        % Save CI energy
        obj.xi  =   0;
    end
    % Common initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end