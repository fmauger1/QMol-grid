function mem = getMemoryProfileWaveFunction(obj,opt)
%getMemoryProfileOrbitalDensity
    
    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nargin < 2,  opt     =   false;  end

    if (obj.sWfcn || obj.sWfcnP)   &&   opt
        QMol_SE_profiler.showMemoryFootprint('Wave functions', 0,1);
    end

    mem                 =   0;
    
    % Wave functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sWfcn
        % How much memory
        ind             =   obj.getOutputIndex(obj.sWfcn);
        if ischar(obj.sWfcnI),              if strcmpi(obj.sWfcnI,'all')
            m           =   obj.SE.N;       else
            m           =   0;              end
        else
            m           =   numel(obj.sWfcnI);
        end
        m               =   QMol_SE_profiler.getMemoryFootprint(m*prod(obj.SE.disc.SE_sizeWaveFunction)*numel(ind),'imag');
        mem             =   mem + m;

        % Show results
        if opt,     QMol_SE_profiler.showMemoryFootprint('Wave functions',m,2);  end
    end

    % Projection of the wave functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sWfcnP
        % Number of time saved
        ind             =   obj.getOutputIndex(obj.sWfcnPT);                        if isempty(obj.sWfcnB)
        % Projection basis
        d           =   QMol_disc_basis('x',obj.SE.disc.x,'v',obj.SE.wfcn.wfcn);    else
        d           =   obj.sWfcnPB;                                                end
                                                                if ischar(obj.sWfcnI),  if strcmpi(obj.sWfcnI,'all') %#ok<ALIGN> 
        % Number of orbitals
        m           =   obj.SE.N;               else,   m   =   0;                      end,    else
        m           =   numel(obj.sWfcnI);                                                      end
        
        % Memory footprint
        m           =   prod(d.SE_sizeWaveFunction) * m * numel(ind);
        m           =   QMol_SE_profiler.getMemoryFootprint(m,'imag') + ...
                        d.getMemoryProfile(false);
        mem         =   mem + m;

        % Show result
        if opt,     QMol_SE_profiler.showMemoryFootprint('Wave function projection',m,2);  end
    end

end