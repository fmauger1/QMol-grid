function mem = getMemoryProfileOrbitalDensity(obj,opt)
%getMemoryProfileOrbitalDensity
    
    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nargin < 2,  opt     =   false;  end

    if (obj.sRho || obj.sKSO || obj.sKSOP)   &&   opt
        QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham orbitals and one-body density', 0,1);
    end

    mem                 =   0;
    
    % Kohn-Sham orbitals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sKSO
        % How much memory
        ind             =   obj.getOutputIndex(obj.sKSOT);
        if ischar(obj.sKSOI),   if strcmpi(obj.sKSOI,'all'),                if obj.DFT.isSpinPol
            m           =   numel(obj.DFT.occ{1}) + numel(obj.DFT.occ{2});  else
            m           =   numel(obj.DFT.occ);                             end,    else
            m           =   0;  end
        elseif iscell(obj.sKSOI),                                           if obj.DFT.isSpinPol
            m           =   numel(obj.sKSOI{1}) + numel(obj.sKSOI{2});      end
        else,                                                               if ~obj.DFT.isSpinPol
            m           =   numel(obj.sKSOI);                               else
            m           =   0;                                              end
        end
        m               =   QMol_DFT_profiler.getMemoryFootprint(m*prod(obj.DFT.disc.DFT_sizeOrbital)*numel(ind),'imag');
        mem             =   mem + m;

        % Show results
        if opt,     QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham orbitals',m,2);  end
    end

    % Projection of the Kohn-Sham orbitals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sKSOP
        % Number of time saved
        ind             =   obj.getOutputIndex(obj.sKSOPT);     if obj.DFT.isSpinPol,   if isempty(obj.sKSOPB) %#ok<ALIGN> 

        % Projection basis
        dUp         =   QMol_disc_basis('x',obj.DFT.disc.x,'v',obj.DFT.KSO.KSOup);
        dDw         =   QMol_disc_basis('x',obj.DFT.disc.x,'v',obj.DFT.KSO.KSOdw);      elseif iscell(obj.sKSOPB)
        dUp         =   obj.sKSOPB{1};          dDw         =   obj.sKSOPB{2};          else
        dUp         =   obj.sKSOPB;             dDw         =   dUp;                    end
                                                                if ischar(obj.sKSOI),   if strcmpi(obj.sKSOI,'all') %#ok<ALIGN> 
        % Number of orbitals
        mUp         =   numel(obj.DFT.occ{1});  mDw         =   numel(obj.DFT.occ{2});  else
        mUp         =   0;                      mDw         =   0;                      end,    elseif iscell(obj.sKSOI)
        mUp         =   numel(obj.sKSOI{1});    mDw         =   numel(obj.sKSOI{1});            else
        mUp         =   0;                      mDw         =   0;                              end
        
        % Memory footprint
        m           =   prod(dUp.DFT_sizeOrbital) * mUp * numel(ind) + ...
                        prod(dDw.DFT_sizeOrbital) * mDw * numel(ind);
        m           =   QMol_DFT_profiler.getMemoryFootprint(m,'imag') + ...
                        dUp.getMemoryProfile(false) + dDw.getMemoryProfile(false);
        mem         =   mem + m;

        % Show result
        if opt,     QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham orbital projections',m,2);  end
                                                                else,                   if isempty(obj.sKSOPB)
        % Projection basis
        d           =   QMol_disc_basis('x',obj.DFT.disc.x,'v',obj.DFT.KSO.KSO);        else
        d           =   obj.sKSOPB;                                                     end
                                                                if ischar(obj.sKSOI),   if strcmpi(obj.sKSOI,'all') %#ok<ALIGN> 
        % Number of orbitals
        m           =   numel(obj.DFT.occ);     else,   m   =   0;                      end,    else
        m           =   numel(obj.sKSOI);                                                       end
        
        % Memory footprint
        m           =   prod(d.DFT_sizeOrbital) * m * numel(ind);
        m           =   QMol_DFT_profiler.getMemoryFootprint(m,'imag') + ...
                        d.getMemoryProfile(false);
        mem         =   mem + m;

        % Show result
        if opt,     QMol_DFT_profiler.showMemoryFootprint('Kohn-Sham orbital projections',m,2);  end
                                                end
    end

    % One-body density ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sRho
        % How much memory
        ind             =   obj.getOutputIndex(obj.sRhoT);
        m               =   QMol_DFT_profiler.getMemoryFootprint(prod(obj.DFT.disc.DFT_sizeOrbital)*numel(ind),'real');
        if obj.DFT.isSpinPol,       m   =   2*m;            end
        mem             =   mem + m;
        
        % Show results
        if opt,     QMol_DFT_profiler.showMemoryFootprint('One-body density',m,2);  end
    end

end