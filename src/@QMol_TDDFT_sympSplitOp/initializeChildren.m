function initializeChildren(obj,isRst)
%initializeChildren

    % Common initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Miscellaneous
        if obj.DFT.isSpinPol,   obj.nKSO    =   [numel(obj.DFT.occ{1}),numel(obj.DFT.occ{2})];
        else,                   obj.nKSO    =   numel(obj.DFT.occ);             end

        obj.dv          =   obj.DFT.disc.dx;                if obj.DFT.dim > 1 %#ok<ALIGN> 
        obj.dv          =   obj.dv * obj.DFT.disc.dy;       if obj.DFT.dim > 2
        obj.dv          =   obj.dv * obj.DFT.disc.dz;       end, end

        % Run time variable
        obj.c_          =   [0 0];

    if isRst
    % Restart ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Projection basis
        if obj.sKSOP,                                                       if iscell(obj.oKSOP.basis)
            obj.oKSOP.basis{1}.initialize;
            obj.oKSOP.basis{2}.initialize;                                  else
            obj.oKSOP.basis.initialize;                                     end
        end

        % Run-time variables
        obj.c_          =   obj.restart.c_;
        obj.expT        =   obj.restart.expT;
        obj.expV        =   obj.restart.expV;
        obj.expVup      =   obj.restart.expVup;
        obj.expVdw      =   obj.restart.expVdw;


    else
    % From scratch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Check DFT functionals
        if iscell(obj.DFT.Vxc),     for k = 1:numel(obj.DFT.Vxc),       if ~startsWith(obj.DFT.Vxc{k}.type,{'LDA','GGA'})   %#ok<ALIGN> 
            warning('QMol:TDDFT:DFTfunctional', ...
                [obj.DFT.Vxc{k}.type '-type functional detected, propagation may not be accurate']);    end, end
        elseif ~isempty(obj.DFT.Vxc),                                   if ~startsWith(obj.DFT.Vxc.type,{'LDA','GGA'})      %#ok<ALIGN> 
            warning('QMol:TDDFT:DFTfunctional', ...
                [obj.DFT.Vxc.type '-type functional detected, propagation may not be accurate']);       end
        end

        % Integration motif
        if isempty(obj.splitMotif),             obj.isVTV   =   true;
        elseif strcmpi(obj.splitMotif,'vtv'),   obj.isVTV   =   true;
        elseif strcmpi(obj.splitMotif,'tvt'),   obj.isVTV   =   false;
        else,                                   obj.isVTV   =   true;
            warning('QMol:TDDFT:splitMotif', ...
                ['Unknown motif for symplectic split ' obj.splitMotif '; default (VTV) used instead.'])
        end

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

        if ~isempty(obj.EF),        obj.EF.initialize(obj.DFT.disc);        end

        % Choice of gauge
        if isempty(obj.EF) || (isempty(obj.EF.electricField)   &&   isempty(obj.EF.potentialVector)) %#ok<ALIGN> 
                                                                        obj.FG  =   0;              % No usable field
        elseif isempty(obj.EFG)                                             % Default mode (context based)
            if isempty(obj.EF),                                         obj.FG  =   0;  %#ok<ALIGN> % No field to work with
            elseif isa(obj.EF.electricField,'function_handle'),         obj.FG  =   1;              % Length gauge
            elseif isa(obj.EF.potentialVector,'function_handle'),       obj.FG  =   2;              % Velocity gauge
            elseif ~isempty(obj.EF.electricField),                      obj.FG  =   1;              % Length gauge
            elseif ~isempty(obj.EF.potentialVector),                    obj.FG  =   2;              % Velocity gauge
            else,                                                       obj.FG  =   0;  end         % No field to work with
        else, switch lower(obj.EFG)                         %#ok<ALIGN> 
            case {'length','length gauge','length_gauge','e'},          obj.FG  =   1;
            case {'velocity','velocity gauge','velocity_gauge','a'},    obj.FG  =   2;
            case {'none','off'},                                        obj.FG  =   0;
            otherwise
                warning('QMol:TDDFT:gaugeField', ['Unknown gauge-field option ' obj.Gfield '; Default gauge mode used instead.'])
    
                if isempty(obj.EF),                                     obj.FG  =   0;  %#ok<ALIGN> % No field to work with
                elseif isa(obj.EF.electricField,'function_handle'),     obj.FG  =   1;              % Length gauge
                elseif isa(obj.EF.potentialVector,'function_handle'),   obj.FG  =   2;              % Velocity gauge
                elseif ~isempty(obj.EF.electricField),                  obj.FG  =   1;              % Length gauge
                elseif ~isempty(obj.EF.potentialVector),                obj.FG  =   2;              % Velocity gauge
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
            obj.FE      =   obj.getFE( obj.tspan(1));                       if obj.sEDFT
            obj.FDE     =   obj.getFDE(obj.tspan(1));                       end

        elseif obj.FG == 2,                                                 if ~obj.uE && ~obj.uA
            % What to use to compute fields (velocity gauge)
            if ~isempty(obj.EF.potentialVector),    obj.uA  =   true;   else,   obj.uE  =   true;   end
                                                                            end, if obj.uA   ||   ~isempty(obj.EF.potentialVector) %#ok<ALIGN> 
            % Initialize relevant fields (trying to match any A input)
            obj.FA      =   obj.EF.potentialVector(obj.tspan(1));       else,   obj.FA  =   0;  end,    if obj.sAcc   ||   obj.sEDFT
            obj.FE      =   obj.getFE( obj.tspan(1));                                                       end
        end

        % Save DFT emergy
        obj.xi  =   0;
    end
    % Common initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Miscellaneous
        if (obj.FG == 1)   ||   obj.sDip,                                                   if     obj.DFT.dim == 1         %#ok<ALIGN> 
            obj.X               =   obj.DFT.disc.x(:);                                      elseif obj.DFT.dim == 2
            [obj.X,obj.Y]       =   meshgrid(obj.DFT.disc.x,obj.DFT.disc.y);                elseif obj.DFT.dim == 3
            [obj.X,obj.Y,obj.Z] =   meshgrid(obj.DFT.disc.x,obj.DFT.disc.y,obj.DFT.disc.z);
        end, end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end