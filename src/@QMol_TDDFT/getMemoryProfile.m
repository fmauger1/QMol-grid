function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint of the TDDFT object with all its components initialized and
%   used.
    
    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nargin < 2,  opt     =   false;  end

    % DFT model memory footprint ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ~obj.DFT.isInit,     obj.DFT.initialize;                             end

    fprintf('\n    ------------------------ DFT model ------------------------\n');
    mem                 =   obj.DFT.getMemoryProfile(opt);

    % TDDFT propagator footprint ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('\n    --------------------- TDDFT propagator --------------------\n');
    mem                 =   mem + obj.getMemoryProfilePropagator(opt);

    if ~isempty(obj.EF)     % External field (returned to original state)
        disc            =   obj.EF.disc;        obj.EF.disc     =   obj.DFT.disc;
        mem             =   mem + obj.EF.getMemoryProfilePropagator(opt);
        obj.EF.disc     =   disc;
    end

    if ~isempty(obj.ABC)    % Absorber (returned to original state)
        DFT             =   obj.ABC.DFT;        obj.ABC.DFT     =   obj.DFT;
        mem             =   mem + obj.ABC.getMemoryProfile(opt);
        obj.ABC.DFT     =   DFT;
    end
    
    fprintf('\n    ----------------------- TDDFT output ----------------------\n');

    TSPAN               =   obj.tspan;      IREF        =   obj.iref;
    obj.tspan           =   uniquetol([obj.T(1):obj.dt:obj.T(end), obj.T(end)],1e-10,'DataScale',1);
    obj.iref            =   obj.getOutputIndex(obj.T);
    
    if obj.sEDFT   ||   obj.sEKSO %----------------------------------------
        if opt, QMol_DFT_profiler.showMemoryFootprint('DFT and orbital energies', 0,1); end

        % DFT energy
        ind             =   obj.getOutputIndex(obj.sEDFTT);
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind),'real');
        if obj.DFT.isSpinPol,   m   =   12*m;   else,   m   =   9*m;        end

        mem             =   mem + m;
        if opt, QMol_DFT_profiler.showMemoryFootprint('DFT energy', m,2); end

        % Orbital energy
        ind             =   obj.getOutputIndex(obj.sEKSOT);
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind),'real');
        if obj.DFT.isSpinPol,   m   =   (sum(obj.DFT.occ{1})+sum(obj.DFT.occ{2})+2)*m;   
        else,                   m   =   (sum(obj.DFT.occ                       )+2)*m;  end

        mem             =   mem + m;
        if opt, QMol_DFT_profiler.showMemoryFootprint('Orbital energy', m,2); end
    end

    if obj.sDip   ||   obj.sVel   ||   obj.sAcc %--------------------------
        if opt, QMol_DFT_profiler.showMemoryFootprint('Dipole, dipole velocity, and dipole acceleration', 0,1); end
                                                                            if obj.sDip
        % Total dipole
        ind             =   obj.getOutputIndex(obj.sDipT);
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind),'real');
        if obj.DFT.isSpinPol,   m   =   (3*obj.DFT.dim+2)*m;   
        else,                   m   =   (  obj.DFT.dim+2)*m;                end

        mem             =   mem + m;
        if opt, QMol_DFT_profiler.showMemoryFootprint('Total dipole signal', m,2); end,                 if ~isempty(obj.sDipI)

        % Orbital-resolved dipole
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sDipI), if strcmpi(obj.sDipI,'all')                   %#ok<ALIGN> 
            if obj.DFT.isSpinPol,   m   =   obj.DFT.dim*(sum(obj.DFT.occ{1})+sum(obj.DFT.occ{2}))*m;
            else,                   m   =   obj.DFT.dim*(sum(obj.DFT.occ                       ))*m;    end
            else,                   m   =   0;                                              end, else
            if obj.DFT.isSpinPol,   m   =   obj.DFT.dim*(sum(obj.sDipI{1})+sum(obj.sDipI{2})    )*m;
            else,                   m   =   obj.DFT.dim*(sum(obj.sDipI                        ) )*m;    end
        end

        mem             =   mem + m;                
        if opt && m>0,      QMol_DFT_profiler.showMemoryFootprint('Orbital-resolved dipole', m,2);  end, end
                                                                            end
                                                                            if obj.sVel
        % Total dipole velocity
        if ischar(obj.sVelT) && any(strcmpi(obj.sVelT,{'dipole','dip'}))    %#ok<ALIGN> 
                                ind =   obj.getOutputIndex(obj.sDipT);
        else,                   ind =   obj.getOutputIndex(obj.sVelT);      end
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind),'real');
        if obj.DFT.isSpinPol,   m   =   (3*obj.DFT.dim+2)*m;   
        else,                   m   =   (  obj.DFT.dim+2)*m;                end

        mem             =   mem + m;
        if opt, QMol_DFT_profiler.showMemoryFootprint('Total dipole-velocity signal', m,2); end,                if ~isempty(obj.sDipI)

        % Orbital-resolved dipole velocity
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sVelI), if strcmpi(obj.sVelI,'all')                   %#ok<ALIGN> 
            if obj.DFT.isSpinPol,   m   =   obj.DFT.dim*(sum(obj.DFT.occ{1})+sum(obj.DFT.occ{2}))*m;
            else,                   m   =   obj.DFT.dim*(sum(obj.DFT.occ                       ))*m;    end
            else,                   m   =   0;                                              end, else
            if obj.DFT.isSpinPol,   m   =   obj.DFT.dim*(sum(obj.sVelI{1})+sum(obj.sVelI{2})    )*m;
            else,                   m   =   obj.DFT.dim*(sum(obj.sVelI                        ) )*m;    end
        end

        mem             =   mem + m;                
        if opt && m>0,      QMol_DFT_profiler.showMemoryFootprint('Orbital-resolved dipole velocity', m,2); end, end
                                                                            end
                                                                            if obj.sAcc
        % Total dipole acceleration
        if ischar(obj.sAccT) && any(strcmpi(obj.sAccT,{'dipole','dip'}))    %#ok<ALIGN> 
                                ind =   obj.getOutputIndex(obj.sDipT);
        else,                   ind =   obj.getOutputIndex(obj.sAccT);      end
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind),'real');
        if obj.DFT.isSpinPol,   m   =   (3*obj.DFT.dim+2)*m;   
        else,                   m   =   (  obj.DFT.dim+2)*m;                end

        mem             =   mem + m;
        if opt, QMol_DFT_profiler.showMemoryFootprint('Total dipole-acceleration signal', m,2); end,                if ~isempty(obj.sDipI)

        % Orbital-resolved dipole acceleration
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sAccI), if strcmpi(obj.sAccI,'all')                   %#ok<ALIGN> 
            if obj.DFT.isSpinPol,   m   =   obj.DFT.dim*(sum(obj.DFT.occ{1})+sum(obj.DFT.occ{2}))*m;
            else,                   m   =   obj.DFT.dim*(sum(obj.DFT.occ                       ))*m;    end
            else,                   m   =   0;                                              end, else
            if obj.DFT.isSpinPol,   m   =   obj.DFT.dim*(sum(obj.sAccI{1})+sum(obj.sAccI{2})    )*m;
            else,                   m   =   obj.DFT.dim*(sum(obj.sAccI                        ) )*m;    end
        end

        mem             =   mem + m;                
        if opt && m>0,      QMol_DFT_profiler.showMemoryFootprint('Orbital-resolved dipole acceleration', m,2); end, end
                                                                            end
    end

    if obj.sIon %----------------------------------------------------------
        if opt, QMol_DFT_profiler.showMemoryFootprint('Ionization statistics', 0,1); end

        % Total ionization
        ind             =   obj.getOutputIndex(obj.sIonT);
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind),'real');
        if obj.DFT.isSpinPol,   m   =   5*m;    else,   m   =   4*m;        end

        mem             =   mem + m;
        if opt, QMol_DFT_profiler.showMemoryFootprint('Total ionization signal', m,2); end,                         if ~isempty(obj.sIKSOI)

        % Orbital-resolved ionization
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sIKSOI), if strcmpi(obj.sIKSOI,'all')                 %#ok<ALIGN> 
            if obj.DFT.isSpinPol,   m   =   (sum(obj.DFT.occ{1})+sum(obj.DFT.occ{2}))*m;
            else,                   m   =   (sum(obj.DFT.occ                       ))*m;    end
            else,                   m   =   0;                                              end, else
            if obj.DFT.isSpinPol,   m   =   (sum(obj.sIKSOI{1})+sum(obj.sIKSOI{2})  )*m;
            else,                   m   =   (sum(obj.sIKSOI                       ) )*m;    end
        end

        mem             =   mem + m;                
        if opt && m>0,      QMol_DFT_profiler.showMemoryFootprint('Orbital-resolved ionization signal', m,2);  end, end
    end

    mem                 =   mem + obj.getMemoryProfileOrbitalDensity(opt);

    if isa(obj.sFRho,'function_handle')   ||   isa(obj.sFKSO,'function_handle')
        if opt, QMol_DFT_profiler.showMemoryFootprint('Output functions', 0,1); end

        % Output function of the density
        ind             =   obj.getOutputIndex(obj.sFRhoT);

        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);
        R               =   obj.sFRho(obj.DFT.rho,obj.tspan(1));                                    if isreal(R)
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind)*(numel(R)+1),'real');   else
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind)* numel(R)   ,'imag') + ...
                            QMol_DFT_profiler.getMemoryFootprint(numel(ind)             ,'real');   end

        mem             =   mem + m;                
        if opt,             QMol_DFT_profiler.showMemoryFootprint('Output function of the density', m,2);  end


        % Output function of the orbitals
        ind             =   obj.getOutputIndex(obj.sFKSOT);

        R               =   obj.sFKSO(obj.DFT.KSO,obj.tspan(1));                                    if isreal(R)
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind)*(numel(R)+1),'real');   else
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(ind)* numel(R)   ,'imag') + ...
                            QMol_DFT_profiler.getMemoryFootprint(numel(ind)             ,'real');   end

        mem             =   mem + m;                
        if opt,             QMol_DFT_profiler.showMemoryFootprint('Output function of the orbitals', m,2);  end
    end


    % Children output
    mem                 =   mem + obj.getMemoryProfileChildren(opt);

    % Finalize ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('\n');
    obj.tspan           =   TSPAN;          obj.iref    =   IREF;

end