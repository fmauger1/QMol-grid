function setOutputOrbitalDensity(obj,opt)
%setOutputOrbitalDensity
 
    switch lower(opt)
    case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%
     
        % Clean up any old data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        obj.oRho            =   [];
        obj.oKSO            =   [];
        obj.oKSOP           =   [];
        
        % One-body density ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if obj.sRho
            % When to save the one-body density
            obj.oRho.ind    =   [obj.getOutputIndex(obj.sRhoT),NaN];
            obj.oRho.time   =   obj.tspan(obj.oRho.ind(1:end-1));
            obj.oRho.n      =   1;

            if obj.sEF,         obj.setOutputExternalField('oRho','init');  end
     
            % Density output initialization
            S               =   obj.DFT.disc.DFT_sizeOrbital;
            obj.oRho.shape  =   S; 
            if numel(S) > 1   &&   S(end) == 1,     S   =   S(1:end-1);     end,if obj.DFT.isSpinPol
            
            obj.oRho.totalUp=   NaN([S numel(obj.oRho.time)]);
            obj.oRho.totalDown= NaN([S numel(obj.oRho.time)]);                  else
            obj.oRho.total  =   NaN([S numel(obj.oRho.time)]);                  end

        else
            obj.oRho.shape  =   [];
            obj.oRho.ind    =   NaN;
            obj.oRho.n      =   1;
        end
        
        % Kohn-Sham orbitals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if obj.sKSO
            % When to save the DFT energies
            obj.oKSO.ind    =   [obj.getOutputIndex(obj.sKSOT),NaN];
            obj.oKSO.time   =   obj.tspan(obj.oKSO.ind(1:end-1));
            obj.oKSO.n      =   1;

            if obj.sEF,         obj.setOutputExternalField('oKSO','init');  end
            
            S               =   obj.DFT.disc.DFT_sizeOrbital;
            if numel(S) > 1   &&   S(end) == 1,     S   =   S(1:end-1);     end

            % Kohn-Sham orbitals output initialization
            if obj.DFT.isSpinPol
                % Indexes
                if isempty(obj.sKSOI),          obj.oKSO.indexOrbitalUp     =   [];
                                                obj.oKSO.indexOrbitalDown   =   [];
                elseif ischar(obj.sKSOI),                                                                   if strcmpi(obj.sKSOI,'all')
                                                obj.oKSO.indexOrbitalUp     =   1:numel(obj.DFT.occ{1});
                                                obj.oKSO.indexOrbitalDown   =   1:numel(obj.DFT.occ{2});    else
                                                obj.oKSO.indexOrbitalUp     =   [];
                                                obj.oKSO.indexOrbitalDown   =   [];
                                                warning('QMol:TDDFT:saveKSO', ...
                                                    ['Unknown option for the saving of the Kohn-Sham orbitals ' obj.sKSOI ...
                                                    '. No orbitals saved.']);                               end
                elseif iscell(obj.sKSOI),       obj.oKSO.indexOrbitalUp     =   obj.sKSOI{1};
                                                obj.oKSO.indexOrbitalDown   =   obj.sKSOI{2};
                else,                           obj.oKSO.indexOrbitalUp     =   [];
                                                obj.oKSO.indexOrbitalDown   =   [];
                                                warning('QMol:TDDFT:saveKSO', ...
                                                    'Unable to determne which orbitals to save; no orbitals saved.');
                end

                % Allocate output
                obj.oKSO.shape          =   {[S numel(obj.oKSO.indexOrbitalUp)],[S numel(obj.oKSO.indexOrbitalDown)]};
                obj.oKSO.orbitalUp      =   NaN([obj.oKSO.shape{1} numel(obj.oKSO.time)],'like',1i);
                obj.oKSO.orbitalDown    =   NaN([obj.oKSO.shape{2} numel(obj.oKSO.time)],'like',1i);
            else
                % Indexes
                if isempty(obj.sKSOI),          obj.oKSO.indexOrbital       =   [];
                elseif ischar(obj.sKSOI),                                                                   if strcmpi(obj.sKSOI,'all')
                                                obj.oKSO.indexOrbital       =   1:numel(obj.DFT.occ);       else
                                                obj.oKSO.indexOrbital       =   [];
                                                warning('QMol:TDDFT:saveKSO', ...
                                                    ['Unknown option for the saving of the Kohn-Sham orbitals ' obj.sKSOI ...
                                                    '. No orbitals saved.']);                               end
                else,                           obj.oKSO.indexOrbital       =   obj.sKSOI;
                end

                % Allocate output
                obj.oKSO.shape          =   [S numel(obj.oKSO.indexOrbital)];
                obj.oKSO.orbital        =   NaN([obj.oKSO.shape numel(obj.oKSO.time)],'like',1i);

            end
        else
            obj.oKSO.shape  =   [];
            obj.oKSO.ind    =   NaN;
            obj.oKSO.n      =   1;
        end
        
        % Kohn-Sham orbitals' projection ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if obj.sKSOP
            % When to save the DFT energies
            obj.oKSOP.ind   =   [obj.getOutputIndex(obj.sKSOPT),NaN];
            obj.oKSOP.time  =   obj.tspan(obj.oKSOP.ind(1:end-1));
            obj.oKSOP.n     =   1;

            if obj.sEF,         obj.setOutputExternalField('oKSOP','init');     end
     
            % Kohn-Sham orbitals' projection output initialization
            if obj.DFT.isSpinPol
                % Projection basis
                if isempty(obj.sKSOPB)                                  % Project on initial KSO
                    obj.oKSOP.basis     =  {QMol_disc_basis('x',obj.DFT.disc.x,'v',obj.DFT.KSO.KSOup), ...
                                            QMol_disc_basis('x',obj.DFT.disc.x,'v',obj.DFT.KSO.KSOdw)};
                    obj.oKSOP.basis{1}.initialize;
                    obj.oKSOP.basis{2}.initialize;
                elseif isnumeric(obj.sKSOPB)                            % Projection basis as a matrix, same basis for up and down spins (assumed compatible)
                    obj.oKSOP.basis     =   QMol_disc_basis('x',obj.DFT.disc.x,'v',obj.sKSOPB);
                    obj.oKSOP.basis.initialize;
                elseif iscell(obj.sKSOPB) && isnumeric(obj.sKSOPB{1})   % Projection bases as a matrices (assumed compatible)
                    obj.oKSOP.basis     =  {QMol_disc_basis('x',obj.DFT.disc.x,'v',obj.sKSOPB{1}), ...
                                            QMol_disc_basis('x',obj.DFT.disc.x,'v',obj.sKSOPB{2})};
                    obj.oKSOP.basis{1}.initialize;
                    obj.oKSOP.basis{2}.initialize;
                else
                    obj.oKSOP.basis     =   obj.sKSOPB;                     if iscell(obj.oKSOP.basis)
                    obj.oKSOP.basis{1}.initialize;
                    obj.oKSOP.basis{2}.initialize;                          else
                    obj.oKSOP.basis.initialize;                             end
                end

                % Indexes 
                if iscell(obj.oKSOP.basis)
                    Sup         =   obj.oKSOP.basis{1}.DFT_sizeOrbital;
                    if numel(Sup) > 1   &&   Sup(end) == 1, Sup =   Sup(1:end-1);   end
                    Sdw         =   obj.oKSOP.basis{2}.DFT_sizeOrbital;
                    if numel(Sdw) > 1   &&   Sdw(end) == 1, Sdw =   Sdw(1:end-1);   end
                else
                    Sup         =   obj.oKSOP.basis.DFT_sizeOrbital;
                    if numel(Sup) > 1   &&   Sup(end) == 1, Sup =   Sup(1:end-1);   end
                    Sdw         =   Sup;
                end

                if isempty(obj.sKSOPI),         obj.oKSOP.indexOrbitalUp    =   [];
                                                obj.oKSOP.indexOrbitalDown  =   [];
                elseif ischar(obj.sKSOPI),                                                                  if strcmpi(obj.sKSOPI,'all')
                                                obj.oKSOP.indexOrbitalUp    =   1:numel(obj.DFT.occ{1});
                                                obj.oKSOP.indexOrbitalDown  =   1:numel(obj.DFT.occ{2});    else
                                                obj.oKSOP.indexOrbitalUp    =   [];
                                                obj.oKSOP.indexOrbitalDown  =   [];
                                                warning('QMol:TDDFT:saveKSOprojection', ...
                                                    ['Unknown option for the saving of the Kohn-Sham orbitals'' projection ' obj.sKSOPI ...
                                                    '. No orbital projections saved.']);                    end
                elseif iscell(obj.sKSOPI),      obj.oKSOP.indexOrbitalUp    =   obj.sKSOPI{1};
                                                obj.oKSOP.indexOrbitalDown  =   obj.sKSOPI{2};
                else,                           obj.oKSOP.indexOrbitalUp    =   [];
                                                obj.oKSOP.indexOrbitalDown  =   [];
                                                warning('QMol:TDDFT:saveKSOprojection', ...
                                                    'Unable to determne which orbitals'' projection to save; no orbital projections saved.');
                end

                % Allocate output
                obj.oKSOP.shape         =   {[Sup numel(obj.oKSOP.indexOrbitalUp)],[Sdw numel(obj.oKSOP.indexOrbitalDown)]};
                obj.oKSOP.orbitalUp     =   NaN([obj.oKSOP.shape{1} numel(obj.oKSOP.time)],'like',1i);
                obj.oKSOP.orbitalDown   =   NaN([obj.oKSOP.shape{2} numel(obj.oKSOP.time)],'like',1i);
            else
                % Projection basis
                if isempty(obj.sKSOPB)          % Project on initial KSO
                    obj.oKSOP.basis     =   QMol_disc_basis('x',obj.DFT.disc.x,'v',obj.DFT.KSO.KSO);
                    obj.oKSOP.basis.initialize;
                elseif isnumeric(obj.sKSOPB)    % Projection basis as a matrix (assumed compatible)
                    obj.oKSOP.basis     =   QMol_disc_basis('x',obj.DFT.disc.x,'v',obj.sKSOPB);
                    obj.oKSOP.basis.initialize;
                else                            % Projection basis as QMol_disc_basis
                    obj.oKSOP.basis     =   obj.sKSOPB;
                    obj.oKSOP.basis.initialize;
                end

                % Indexes 
                S               =   obj.oKSOP.basis.DFT_sizeOrbital;
                if numel(S) > 1   &&   S(end) == 1,     S   =   S(1:end-1);   end

                if isempty(obj.sKSOPI),         obj.oKSOP.indexOrbital      =   [];
                elseif ischar(obj.sKSOPI),                                                                  if strcmpi(obj.sKSOPI,'all')
                                                obj.oKSOP.indexOrbital      =   1:numel(obj.DFT.occ);       else
                                                obj.oKSOP.indexOrbital      =   [];
                                                warning('QMol:TDDFT:saveKSOprojection', ...
                                                    ['Unknown option for the saving of the Kohn-Sham orbitals'' projection ' obj.sKSOPI ...
                                                    '. No orbital projections saved.']);                    end
                else,                           obj.oKSOP.indexOrbital      =   obj.sKSOPI;
                end

                % Allocate output
                obj.oKSOP.shape         =   [S numel(obj.oKSOP.indexOrbital)];
                obj.oKSOP.orbital       =   NaN([obj.oKSOP.shape numel(obj.oKSOP.time)],'like',1i);

            end
        else
            obj.oKSOP.shape =   [];
            obj.oKSOP.ind   =   NaN;
            obj.oKSOP.n     =   1;
        end
        
    case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if obj.sRho,    obj.oRho    =   rmfield(obj.oRho,{'ind','n','shape'});  if obj.sEF %#ok<ALIGN> 
                        obj.setOutputExternalField('oRho' ,'clean');            end
        else,           obj.oRho    =   [];                                     end
        if obj.sKSO,    obj.oKSO    =   rmfield(obj.oKSO,{'ind','n','shape'});  if obj.sEF %#ok<ALIGN> 
                        obj.setOutputExternalField('oKSO' ,'clean');            end
        else,           obj.oKSO    =   [];                                     end
        if obj.sKSOP,   obj.oKSOP   =   rmfield(obj.oKSOP,{'ind','n','shape'}); if obj.sEF %#ok<ALIGN> 
                        obj.setOutputExternalField('oKSOP','clean');            end
        else,           obj.oKSOP   =   [];                                     end
     
    otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Unexpected option
        error('QMol:TDDFT:setOutputOrbitalDensity',['Unknown option ' opt]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end