function mem = getMemoryProfileWaveFunction(obj,opt)
%getMemoryProfileOrbitalDensity
    
    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nargin < 2,  opt     =   false;  end

    if obj.sWfcn   &&   opt
        QMol_SE_profiler.showMemoryFootprint('Wave functions', 0,1);
    end

    mem                 =   0;
    
    % Wave functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sWfcn
        % How much memory
        ind             =   obj.getOutputIndex(obj.sWfcn);
        if ischar(obj.sWfcnI),                      if strcmpi(obj.sWfcnI,'all')
            m           =   size(obj.CI.wfcn,2);    else
            m           =   0;                      end
        else
            m           =   numel(obj.sWfcnI);
        end
        m               =   QMol_SE_profiler.getMemoryFootprint(m*size(obj.CI.wfcn,1)*numel(ind),'imag');
        mem             =   mem + m;

        % Show results
        if opt,     QMol_SE_profiler.showMemoryFootprint('Wave functions',m,2);  end
    end

end