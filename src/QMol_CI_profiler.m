classdef QMol_CI_profiler < QMol_profiler
%QMol_CI_profiler time/memory profiler for configuration-interaction models
%
%   QMol_CI_profiler(CI) or QMol_CI_profiler(CI,'all') displays an estimate
%   of the memory footprint of the configuration interaction object CI. If
%   the object defines a CI or dipole-coupling matrices, the profiling also
%   estimates their averate product time with a wave function.
%
%   QMol_CI_profiler(TDCI) or QMol_CI_profiler(TDCI,'all') performs the
%   profiling on the time-dependent configuration interaction propagator
%   object TDCI.
%
%   QMol_CI_profiler(___,'memory') only performs the memory profiling.
%
%   QMol_CI_profiler(___,'time') only performs the time profiling.
%
%   See also QMol_CI_conv
    
%   Version     Date        Author
%   01.23.000   05/25/2025  F. Mauger
%       Creation (from QMol_SE_profiler version 01.23.000)
%   01.23.001   06/05/2025  F. Mauger
%       Handle time profiling when no CI or dipole-coupling matrices

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.001','06/05/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_profiler})
function showInfo
    fprintf('  * QMol_CI_profiler:\n');
    fprintf('      > Memory/execution time profiling for CI\n'); 
    QMol_CI_profiler.version;
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties
    mode                =   'all'
    nbIter              =   []
end

%% Run profiling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
    function obj = QMol_CI_profiler(QM,opt,N)
    
    % CI profiler header ==================================================
    fprintf(['########################### QMol-grid package ############################\n'...
             'Quantum simulation methods for atomic and molecular systems [Mauger 2024b]\n\n' ...
             '  * Configuration interaction (CI) profiler (release) [Visentin 2025]\n']);
    QMol_CI_profiler.version;       fprintf('\n');
    QMol_doc.showFooter;            fprintf('\n');

    QMol_doc.showSection('External components');
    fprintf('  * convertUnit\n'); QMol_doc.showVersion(convertUnit.version,convertUnit.lastMod,convertUnit.author);
    fprintf('  * fourierTool\n'); QMol_doc.showVersion(convertUnit.version,convertUnit.lastMod,convertUnit.author);
    fprintf('\n')

    QMol_doc.showSection('License');
    QMol_doc.showLicense;
    fprintf('\n')

    QMol_doc.showBibliography({'Visentin 2025'});

    QMol_doc.showFunding;

    % Initialization ======================================================
    if nargin > 1,  obj.mode    =   opt;
    if nargin > 2,  obj.nbIter  =   N;          end, end

    % Memory ==============================================================
if any(strcmpi(obj.mode,{'all','memory','size'}))
    % Initialization
    QMol_doc.showSection('Memory footprint');

    fprintf('  Components                                                    -- Size --\n');

    % Get memory footprint
    mem                 =   QM.getMemoryProfile(true);

    % Show total memory footprint
    fprintf('                                                                ----------\n');
    fprintf('                                                        TOTAL =  ');
    QMol_profiler.showMemoryFormatted(mem);
    fprintf('\n\n');

end
    % Execution time ======================================================
if any(strcmpi(obj.mode,{'all','time','execution'}))
    QMol_doc.showSection('Execution time');

    % Initialization
    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility

    if isa(QM,'QMol_TDCI') % TDCI object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fprintf('*** No execution-time estimate for TDCI objects.                       ***\n\n');
    elseif isa(QM,'QMol_CI') % CI object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Initialization
        QM.initialize;      inWfcn  =   QM.wfcn;

        obj.nbIter      =   500;

        if isempty(QM.wfcn),    N   =   1;
        else,                   N   =   size(QM.wfcn,2);                    end
        p               =   rand(randStr,[size(QM.CSB,1),N],'like',1i)-.5-.5i;
        p               =   p / sqrt(sum(abs(p).^2));

        fprintf('  Operator components                                           -- Time --\n');

        % CI matrix
        if ~isempty(QM.CI)
            tic, for k = 1:obj.nbIter,  Hp  =   QM.CI*p;                    end,    tH = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-59s %8.2e s\n','CI matrix',tH);
        else,                                                                       tH = 0;
            fprintf('  * %-59s %8.2e s\n','No CI matrix to time',tH);
        end

        % Dipole-coupling matrix
        if ~isempty(QM.DX)
            tic, for k = 1:obj.nbIter,  Hp  =   QM.DX*p;                    end,    tDX = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-59s %8.2e s\n','Dipole-coupling matrix',tDX);
        else,                                                                       tDX = 0;
            fprintf('  * %-59s %8.2e s\n','No dipole-coupling matrix to test',tDX);
        end

        fprintf('                                                                ----------\n');
        fprintf('                                                        TOTAL = %8.2e s\n',tH+tDX);

        fprintf('\n  (average times for %u wave function(s) over %u iterations)\n\n',N,obj.nbIter);

        % Clean up
        QM.wfcn         =   inWfcn;
        
    end


end
    % Suppress the output object ==========================================
    QMol_doc.showFooter;
    clear obj

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

