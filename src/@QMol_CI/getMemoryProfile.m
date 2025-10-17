function mem = getMemoryProfile(obj,opt)
%getMemoryProfile estimate of the total memory footprint of CI object
%   Use getMemoryProfile to calculate an estimate for the total memory
%   footprint of the CI object with all its components initialized and
%   used. The methods will try to return the various objects in the same
%   configuration they were passed, but it's not 100% guarantied. 
%   getMemoryProfile should be used in isolation and not, e.g., preceding 
%   CI matrix or time propagation simulations.
%
%   mem = CI.getMemoryProfile and CI.getMemoryProfile(false) returns the
%   memory footprint estimate in bytes.
%
%   mem = CI.getMemoryProfile(true) also displays the details of the memory
%   footprint.
    
    % Initialization
    if nargin < 2,  opt     =   false;  end

    % Domain discretization
    disc                =   obj.getDiscCopy;
    mem                 =   disc.getMemoryProfile(opt);

    inDisc              =   obj.disc;
    obj.disc            =   disc;

    % External potential
    if ~isempty(obj.Vext)
        if ~obj.Vext.isInit,    obj.Vext.DFT    =   obj;                    end
        mem             =   mem + obj.Vext.getMemoryProfile(opt);
        if ~obj.Vext.isInit,    obj.Vext.reset();                           end
    end

    % Electron-interaction potential
    m                   =   QMol_DFT_profiler.getMemoryFootprint(2*numel(obj.disc.x)-1,'real');
    mem                 =   mem + m;
    if opt
        QMol_DFT_profiler.showMemoryFootprint('Electron-electron interaction',m,1);
    end

    % Spatial-orbital basis
    if ~isempty(obj.SOB),   if isreal(obj.SOB)                              %#ok<ALIGN>
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.SOB),'real');
    else
        m               =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.SOB),'imag');
    end, end
    mem                 =   mem + m;
    if opt
        QMol_DFT_profiler.showMemoryFootprint('Spatial orbital basis',m,1);
    end

    % Configuration basis, CI and dipole matrices
    if isempty(obj.CSB)
        % Assume the configuration basis will be built with setConfigurationBasis
        try
            % Configuration basis parameters
            [algo,locN,locRef,locAct,locFrz] = obj.getConfigurationBasisParameters(false);
            nbRef       =   size(locRef,1);
    
            % CIS(D) parameters
            if algo > 0,    nCSB = nbRef;                                       for k = 1:nbRef %#ok<ALIGN>
                % Allowed spin-orbital excitations
                iExP    =   locAct{k}(locAct{k} > 0);   nExP    =   numel(iExP);
                iExN    =   locAct{k}(locAct{k} < 0);   nExN    =   numel(iExN);
                % Allowed spin-orbital holes
                [~,iiHo]=   setdiff(locRef(k,:),locFrz);
                nHoP    =   sum(locRef(k,iiHo) > 0);    nHoN    =   sum(locRef(k,iiHo) < 0);
                % Number of basis vectors
                nCSB    =   nCSB + nHoP*nExP + nHoN*nExN;                      % reference plus number of single excitations
                if algo > 1                                                     % double excitations
                    nCSB=   nCSB + nHoP*(nHoP-1)/2 * nExP*(nExP-1)/2;       %  > both in up spin channel
                    nCSB=   nCSB + nHoN*(nHoN-1)/2 * nExN*(nExN-1)/2;       %  > both in down spin channel
                    nCSB=   nCSB + nHoP*nExP * nHoN*nExN;                   %  > in different spin channels
                end, end
            % RAS parameters
            else
                % Active orbitals
                locAct  =   locAct{1};
                iActP   =   locAct(locAct > 0);       nActP   =   sum(locRef > 0) - sum(locFrz > 0);
                iActN   =   locAct(locAct < 0);       nActN   =   sum(locRef < 0) - sum(locFrz < 0);
    
                % Number of configuration states
                nRasP   =   factorial(numel(iActP))/(factorial(numel(iActP)-nActP) * factorial(nActP));
                nRasN   =   factorial(numel(iActN))/(factorial(numel(iActN)-nActN) * factorial(nActN));
                nCSB    =   nRasP*nRasN;
            end
    
            % Display memory footprint
            m           =   QMol_DFT_profiler.getMemoryFootprint(locN*nCSB,'real');
            mem         =   mem + m;
            if opt
                QMol_DFT_profiler.showMemoryFootprint('Configuration state basis',m,1);
            end
    
            m           =   QMol_DFT_profiler.getMemoryFootprint(nCSB^2,'real');
            mem         =   mem + m;
            if opt
                QMol_DFT_profiler.showMemoryFootprint('CI matrix',m,1);
            end
    
            if obj.isDip
                mem     =   mem + m;        % dipole-coupling matrix has the same size as the CI one
                if opt
                    QMol_DFT_profiler.showMemoryFootprint('Dipole-coupling matrix',m,1);
                end
            end
        catch 
            % Nothing to catch
        end
    else
        % Assume the member configuration basis
        nCSB            =   size(obj.CSB,1);
        locN            =   size(obj.CSB,2);
    
        % Display memory footprint
        m               =   QMol_DFT_profiler.getMemoryFootprint(locN*nCSB,'real');
        mem             =   mem + m;
        if opt
            QMol_DFT_profiler.showMemoryFootprint('Configuration state basis',m,1);
        end

        m               =   QMol_DFT_profiler.getMemoryFootprint(nCSB^2,'real');
        mem             =   mem + m;
        if opt
            QMol_DFT_profiler.showMemoryFootprint('CI matrix',m,1);
        end

        if obj.isDip
            mem         =   mem + m;        % dipole-coupling matrix has the same size as the CI one
            if opt
                QMol_DFT_profiler.showMemoryFootprint('Dipole-coupling matrix',m,1);
            end
        end
    end

    % Wave function(s)
    if ~isempty(obj.wfcn)
        if isreal(obj.wfcn),    m   =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.wfcn),'real');
        else,                   m   =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.wfcn),'imag');   end
        mem         =   mem + m;
        if opt
            QMol_DFT_profiler.showMemoryFootprint('Wave function(s)',m,1);
        end
    end

    % Clean up
    obj.disc            =   inDisc;
    clear disc
end