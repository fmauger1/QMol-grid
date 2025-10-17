function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint of the TDSE object with all its components initialized and
%   used.
    
    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nargin < 2,  opt     =   false;  end

    % SE model memory footprint ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ~obj.SE.isInit,      obj.SE.initialize;                              end

    fprintf('\n    ------------------------- SE model ------------------------\n');
    mem                 =   obj.SE.getMemoryProfile(opt);

    % TDSE propagator footprint ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('\n    --------------------- TDSE propagator ---------------------\n');
    mem                 =   mem + obj.getMemoryProfilePropagator(opt);

    if ~isempty(obj.EF)     % External field (returned to original state)
        disc            =   obj.EF.disc;        obj.EF.disc     =   obj.SE.disc;
        mem             =   mem + obj.EF.getMemoryProfilePropagator(opt);
        obj.EF.disc     =   disc;
    end

    if ~isempty(obj.ABC)    % Absorber (returned to original state)
        SE              =   obj.ABC.SE;         obj.ABC.SE      =   obj.SE;
        mem             =   mem + obj.ABC.getMemoryProfile(opt);
        obj.ABC.SE      =   SE;
    end
    
    fprintf('\n    ----------------------- TDSE output -----------------------\n');

    TSPAN               =   obj.tspan;      IREF        =   obj.iref;
    obj.tspan           =   uniquetol([obj.T(1):obj.dt:obj.T(end), obj.T(end)],1e-10,'DataScale',1);
    obj.iref            =   obj.getOutputIndex(obj.T);
    
    if obj.sESE   ||   obj.sEWfcn %----------------------------------------
        if opt, QMol_DFT_profiler.showMemoryFootprint('Schrodinger-equation and wave-function energies', 0,1); end

        % SE energy
        ind             =   obj.getOutputIndex(obj.sESET);
        m               =   5*QMol_SE_profiler.getMemoryFootprint(numel(ind),'real');

        mem             =   mem + m;
        if opt, QMol_DFT_profiler.showMemoryFootprint('Schrodinger-equation energy', m,2); end

        % Wave function energy
        ind             =   obj.getOutputIndex(obj.sEWfcnT);
        m               =   QMol_SE_profiler.getMemoryFootprint(numel(ind),'real');
        m               =   sum(obj.SE.N)*m;

        mem             =   mem + m;
        if opt, QMol_SE_profiler.showMemoryFootprint('Wave function energy', m,2); end
    end

    if obj.sDip   ||   obj.sVel   ||   obj.sAcc %--------------------------
        if opt, QMol_DFT_profiler.showMemoryFootprint('Dipole, dipole velocity, and dipole acceleration', 0,1); end
                                                                            if obj.sDip
        % Dipole
        ind             =   obj.getOutputIndex(obj.sDipT);
        m               =   QMol_SE_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sDipI),                           if strcmpi(obj.sDipI,'all')
            m           =   obj.SE.dim*obj.SE.N*m;      else
            m           =   0;                          end 
        elseif isempty(obj.sDipI)
            m           =   obj.SE.dim*obj.SE.N*m;
        else
            m           =   obj.SE.dim*(sum(obj.sDipI) )*m;
        end
        m_dip           =   m;

        mem             =   mem + m;                
        if opt && m>0,      QMol_SE_profiler.showMemoryFootprint('Dipole', m,2);  end
                                                                            else
        m_dip           =   0;
                                                                            end
                                                                            if obj.sVel
        % Dipole velocity
        if ischar(obj.sVelT) && any(strcmpi(obj.sVelT,{'dipole','dip'}))    %#ok<ALIGN> 
                                ind =   obj.getOutputIndex(obj.sDipT);
        else,                   ind =   obj.getOutputIndex(obj.sVelT);      end
        m               =   QMol_SE_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sVelI),                           if strcmpi(obj.sVelI,'all')
            m           =   obj.SE.dim*obj.SE.N*m;      elseif strcmpi(obj.sVelI,{'dipole','dip'})
            m           =   m_dip;                      else
            m           =   0;                          end 
        elseif isempty(obj.sVelI)
            m           =   obj.SE.dim*obj.SE.N*m;
        else
            m           =   obj.SE.dim*(sum(obj.sVelI) )*m;
        end

        mem             =   mem + m;                
        if opt && m>0,      QMol_SE_profiler.showMemoryFootprint('Dipole velocity', m,2);  end
                                                                            end
                                                                            if obj.sAcc
        % Dipole acceleration
        if ischar(obj.sAccT) && any(strcmpi(obj.sAccT,{'dipole','dip'}))    %#ok<ALIGN> 
                                ind =   obj.getOutputIndex(obj.sDipT);
        else,                   ind =   obj.getOutputIndex(obj.sAccT);      end
        m               =   QMol_SE_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sAccI),                           if strcmpi(obj.sAccI,'all')
            m           =   obj.SE.dim*obj.SE.N*m;      elseif strcmpi(obj.sAccI,{'dipole','dip'})
            m           =   m_dip;                      else
            m           =   0;                          end 
        elseif isempty(obj.sAccI)
            m           =   obj.SE.dim*obj.SE.N*m;
        else
            m           =   obj.SE.dim*(sum(obj.sAccI) )*m;
        end

        mem             =   mem + m;                
        if opt && m>0,      QMol_SE_profiler.showMemoryFootprint('Dipole acceleration', m,2);  end
                                                                            end
    end

    if obj.sIon %----------------------------------------------------------
        if opt, QMol_DFT_profiler.showMemoryFootprint('Ionization statistics', 0,1); end

        % Ionization
        ind             =   obj.getOutputIndex(obj.sIonT);
        m               =   QMol_SE_profiler.getMemoryFootprint(numel(ind),'real');
        if ischar(obj.sIWfcnI),             if strcmpi(obj.sIWfcnI,'all')
            m           =   obj.SE.N*m;     else
            m           =   0;              end
        elseif isempty(obj.sIWfcnI)
            m           =   obj.SE.N*m;
        else
            m           =   (sum(obj.sIWfcnI) )*m;
        end

        mem             =   mem + m;                
        if opt && m>0,      QMol_SE_profiler.showMemoryFootprint('Ionization signal', m,2);  end
    end

    mem                 =   mem + obj.getMemoryProfileWaveFunction(opt);

    if isa(obj.sF,'function_handle')
        if opt, QMol_SE_profiler.showMemoryFootprint('Output function', 0,1); end

        % Output function of the wave functions
        ind             =   obj.getOutputIndex(obj.sFT);

        R               =   obj.sF(obj.SE.wfcn,obj.tspan(1));                                      if isreal(R)
        m               =   QMol_SE_profiler.getMemoryFootprint(numel(ind)*(numel(R)+1),'real');   else
        m               =   QMol_SE_profiler.getMemoryFootprint(numel(ind)* numel(R)   ,'imag') + ...
                            QMol_SE_profiler.getMemoryFootprint(numel(ind)             ,'real');   end

        mem             =   mem + m;                
        if opt,             QMol_SE_profiler.showMemoryFootprint('Output function', m,2);  end
    end


    % Children output
    mem                 =   mem + obj.getMemoryProfileChildren(opt);

    % Finalize ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('\n');
    obj.tspan           =   TSPAN;          obj.iref    =   IREF;

end