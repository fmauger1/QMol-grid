function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint of the TDCI object with all its components initialized and
%   used.
    
    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nargin < 2,  opt     =   false;  end

    % CI model memory footprint ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ~obj.CI.isInit,      obj.CI.initialize;                              end

    fprintf('\n    ------------------------- CI model ------------------------\n');
    mem                 =   obj.CI.getMemoryProfile(opt);

    % TDCI propagator footprint ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('\n    --------------------- TDCI propagator ---------------------\n');
    mem                 =   mem + obj.getMemoryProfilePropagator(opt);

    if ~isempty(obj.EF)     % External field (returned to original state)
        disc            =   obj.EF.disc;        obj.EF.disc     =   obj.CI.disc;
        mem             =   mem + obj.EF.getMemoryProfilePropagator(opt);
        obj.EF.disc     =   disc;
    end

    if ~isempty(obj.ABC)    % Absorber (returned to original state)
        CI              =   obj.ABC.CI;         obj.ABC.CI      =   obj.CI;
        mem             =   mem + obj.ABC.getMemoryProfile(opt);
        obj.ABC.CI      =   CI;
    end
    
    fprintf('\n    ----------------------- TDCI output -----------------------\n');

    TSPAN               =   obj.tspan;      IREF        =   obj.iref;
    obj.tspan           =   uniquetol([obj.T(1):obj.dt:obj.T(end), obj.T(end)],1e-10,'DataScale',1);
    obj.iref            =   obj.getOutputIndex(obj.T);
    
    if obj.sECI   ||   obj.sEWfcn %----------------------------------------
        if opt, QMol_DFT_profiler.showMemoryFootprint('Configuration-interaction and wave-function energies', 0,1); end

        % CI energy
        ind             =   obj.getOutputIndex(obj.sECIT);
        m               =   5*QMol_CI_profiler.getMemoryFootprint(numel(ind),'real');

        mem             =   mem + m;
        if opt, QMol_DFT_profiler.showMemoryFootprint('Configuration interaction energy', m,2); end

        % Wave function energy
        ind             =   obj.getOutputIndex(obj.sEWfcnT);
        m               =   QMol_CI_profiler.getMemoryFootprint(numel(ind),'real');
        m               =   sum(obj.CI.N)*m;

        mem             =   mem + m;
        if opt, QMol_CI_profiler.showMemoryFootprint('Wave function energy', m,2); end
    end

    if obj.sDip   ||   obj.sVel   ||   obj.sAcc %--------------------------
        if opt, QMol_DFT_profiler.showMemoryFootprint('Dipole, dipole velocity, and dipole acceleration', 0,1); end
                                                                            if obj.sDip
        % Dipole
        ind             =   obj.getOutputIndex(obj.sDipT);
        m               =   QMol_CI_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sDipI),                           if strcmpi(obj.sDipI,'all')
            m           =   obj.CI.dim*obj.CI.N*m;      else
            m           =   0;                          end 
        elseif isempty(obj.sDipI)
            m           =   obj.CI.dim*obj.CI.N*m;
        else
            m           =   obj.CI.dim*(sum(obj.sDipI) )*m;
        end
        m_dip           =   m;

        mem             =   mem + m;                
        if opt && m>0,      QMol_CI_profiler.showMemoryFootprint('Dipole', m,2);  end
                                                                            else
        m_dip           =   0;
                                                                            end
                                                                            if obj.sVel
        % Dipole velocity
        if ischar(obj.sVelT) && any(strcmpi(obj.sVelT,{'dipole','dip'}))    %#ok<ALIGN> 
                                ind =   obj.getOutputIndex(obj.sDipT);
        else,                   ind =   obj.getOutputIndex(obj.sVelT);      end
        m               =   QMol_CI_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sVelI),                           if strcmpi(obj.sVelI,'all')
            m           =   obj.CI.dim*obj.CI.N*m;      elseif strcmpi(obj.sVelI,{'dipole','dip'})
            m           =   m_dip;                      else
            m           =   0;                          end 
        elseif isempty(obj.sVelI)
            m           =   obj.CI.dim*obj.CI.N*m;
        else
            m           =   obj.CI.dim*(sum(obj.sVelI) )*m;
        end

        mem             =   mem + m;                
        if opt && m>0,      QMol_CI_profiler.showMemoryFootprint('Dipole velocity', m,2);  end
                                                                            end
                                                                            if obj.sAcc
        % Dipole acceleration
        if ischar(obj.sAccT) && any(strcmpi(obj.sAccT,{'dipole','dip'}))    %#ok<ALIGN> 
                                ind =   obj.getOutputIndex(obj.sDipT);
        else,                   ind =   obj.getOutputIndex(obj.sAccT);      end
        m               =   QMol_CI_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sAccI),                           if strcmpi(obj.sAccI,'all')
            m           =   obj.CI.dim*obj.CI.N*m;      elseif strcmpi(obj.sAccI,{'dipole','dip'})
            m           =   m_dip;                      else
            m           =   0;                          end 
        elseif isempty(obj.sAccI)
            m           =   obj.CI.dim*obj.CI.N*m;
        else
            m           =   obj.CI.dim*(sum(obj.sAccI) )*m;
        end

        mem             =   mem + m;                
        if opt && m>0,      QMol_CI_profiler.showMemoryFootprint('Dipole acceleration', m,2);  end
                                                                            end
    end

    if obj.sIon %----------------------------------------------------------
        if opt, QMol_DFT_profiler.showMemoryFootprint('Ionization statistics', 0,1); end

        % Ionization
        ind             =   obj.getOutputIndex(obj.sIonT);
        m               =   QMol_CI_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sIWfcnI),             if strcmpi(obj.sIWfcnI,'all')
            m           =   obj.CI.N*m;     else
            m           =   0;              end
        elseif isempty(obj.sIWfcnI)
            m           =   obj.CI.N*m;
        else
            m           =   (sum(obj.sIWfcnI) )*m;
        end

        mem             =   mem + m;                
        if opt && m>0,      QMol_CI_profiler.showMemoryFootprint('Ionization signal', m,2);  end
    end

    mem                 =   mem + obj.getMemoryProfileWaveFunction(opt);

    if isa(obj.sF,'function_handle')
        if opt, QMol_CI_profiler.showMemoryFootprint('Output function', 0,1); end

        % Output function of the wave functions
        ind             =   obj.getOutputIndex(obj.sFT);

        R               =   obj.sF(obj.CI,obj.tspan(1));                                            if isreal(R)
        m               =   QMol_CI_profiler.getMemoryFootprint(numel(ind)*(numel(R)+1),'real');    else
        m               =   QMol_CI_profiler.getMemoryFootprint(numel(ind)* numel(R)   ,'imag') +   ...
                            QMol_CI_profiler.getMemoryFootprint(numel(ind)             ,'real');    end

        mem             =   mem + m;                
        if opt,             QMol_CI_profiler.showMemoryFootprint('Output function', m,2);  end
    end


    % Children output
    mem                 =   mem + obj.getMemoryProfileChildren(opt);

    % Finalize ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('\n');
    obj.tspan           =   TSPAN;          obj.iref    =   IREF;

end