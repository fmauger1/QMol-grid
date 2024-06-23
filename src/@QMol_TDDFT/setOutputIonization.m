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

        % Total ionization
        obj.oIon.total                  =   NaN(1,numel(obj.oIon.time));    if obj.DFT.isSpinPol
        obj.oIon.totalUp                =   NaN(1,numel(obj.oIon.time));
        obj.oIon.totalDown              =   NaN(1,numel(obj.oIon.time));    end
        
        % Orbital resolved ionization
        if ischar(obj.sIKSOI), if strcmpi(obj.sIKSOI,'all')                 % Save all orbitals
            if obj.DFT.isSpinPol
                % What to save
                obj.oIon.indexOrbitalUp=   1:numel(obj.DFT.occ{1});
                obj.oIon.indexOrbitalDown= 1:numel(obj.DFT.occ{2});

                % Allocation
                obj.oIon.orbitalUp      =   NaN(numel(obj.DFT.occ{1}),numel(obj.oIon.time));
                obj.oIon.orbitalDown    =   NaN(numel(obj.DFT.occ{2}),numel(obj.oIon.time));
            else
                % What to save
                obj.oIon.indexOrbital   =   1:numel(obj.DFT.occ);
                
                % Allocation
                obj.oIon.orbital        =   NaN(numel(obj.DFT.occ),numel(obj.oIon.time));
            end,              else
                warning('QMol:TDDFT:outputIonization',...
                    ['Unknown option ' obj.sIonI ' for orbital-resolved ionization; feature disabled']);
                if obj.DFT.isSpinPol,   obj.oIon.indexOrbitalUp   =   [];     obj.oIon.indexOrbitalDown   =   [];
                else,                   obj.oIon.indexOrbital     =   [];
                end,          end
        elseif ~isempty(obj.sIonI)                                          % User-supplied indexes
            if obj.DFT.isSpinPol,       obj.oIon.indexOrbitalUp   =   obj.sIonI{1}(:).';
                                        obj.oIon.indexOrbitalDown =   obj.sIonI{2}(:).';
                obj.oIon.orbitalUp      =   NaN(numel(obj.sIonI{1}),numel(obj.oIon.time));
                obj.oIon.orbitalDown    =   NaN(numel(obj.sIonI{2}),numel(obj.oIon.time));
            else,                       obj.oIon.indexOrbital     =   obj.sDipI(:).';
                obj.oIon.orbital        =   NaN(numel(obj.sIonI   ),numel(obj.oIon.time));
            end
        else                                                                % No orbital resolving
            if obj.DFT.isSpinPol,       obj.oIon.indexOrbitalUp   =   [];     obj.oIon.indexOrbitalDown   =   [];
            else,                       obj.oIon.indexOrbital     =   [];
            end
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
    error('QMol:TDDFT:setOutputIonization',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end