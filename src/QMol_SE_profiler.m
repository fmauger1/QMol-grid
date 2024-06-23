classdef QMol_SE_profiler < QMol_suite
%QMol_SE_profiler memory and/or execution time profiling for Schrodinger-
%   equation model or functional
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_SE_profiler})
function showInfo
    fprintf('  * QMol_SE_profiler:\n');
    fprintf('      > memory/execution time profiling for SE\n'); 
    QMol_SE_profiler.version;
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties
    mode                =   'all'
    nbIter              =   []
end

%% Run profiling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
    function obj = QMol_SE_profiler(QM,opt,N)
    
    % SE profiler header ==================================================
    fprintf(['########################### QMol-grid package ############################\n'...
             'Quantum simulation methods for atomic and molecular systems [Mauger XXXX].\n\n' ...
             '  * Schrodinger-equation (SE) profiler (release)\n']);
    QMol_SE_profiler.version;       fprintf('\n');
    QMol_doc.showFooter;            fprintf('\n');

    QMol_doc.showSection('External components');
    fprintf('  * convertUnit\n'); QMol_doc.showVersion(convertUnit.version,convertUnit.lastMod,convertUnit.author);
    fprintf('  * fourierTool\n'); QMol_doc.showVersion(convertUnit.version,convertUnit.lastMod,convertUnit.author);
    fprintf('\n')

    QMol_doc.showSection('License');
    QMol_doc.showLicense;
    fprintf('\n')

    QMol_doc.showBibliography([]);

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
    QMol_SE_profiler.showMemoryFormatted(mem);
    fprintf('\n\n');

end
    % Execution time ======================================================
if any(strcmpi(obj.mode,{'all','time','execution'}))
    QMol_doc.showSection('Execution time');

    % Initialization
    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility

    if isa(QM,'QMol_TDSE') % TDSE object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fprintf('*** No execution-time estimate for TDSE objects.                       ***\n');
    elseif isa(QM,'QMol_SEq') % DFT object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Initialization
        QM.initialize;      inWfcn  =   QM.wfcn;

        if isempty(obj.nbIter)
            if     QM.dim == 1, obj.nbIter  =   500;
            elseif QM.dim == 2, obj.nbIter  =   50;
            elseif QM.dim == 3, obj.nbIter  =   10;
            else,  error('QMol:QMol_SE_profiler:dimSE',['Unexpected dimension ' num2str(QM.dim) ' (1, 2, or 3 allowed)']);
            end
        end

        N           =   QM.N;
        QM.wfcn     =   QM.disc.SE_allocateWaveFunction(N,[],randStr,false);
        p           =   QM.disc.SE_randomWaveFunction(randStr) + 1i*QM.disc.SE_randomWaveFunction(randStr);
        p           =   QM.disc.SE_normalizeWaveFunction(p);

        V           =   QM.V;

        fprintf('  Operator components                                           -- Time --\n');

        % Whole Hamiltonian operator
        tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.SE_operatorHamiltonian(V,p,0);      end, t = toc/obj.nbIter; %#ok<NASGU> 
        fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Hamiltonian'  ,N,t,N*t);

        fprintf('\n');

        % Kinetic & potential operators
        tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.SE_operatorKinetic(p,0);            end, tK = toc/obj.nbIter; %#ok<NASGU> 
        fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Kinetic',N,tK,N*tK);

        tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.SE_operatorPotential(V,p,0);        end, tV = toc/obj.nbIter; %#ok<NASGU> 
        fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Potential',N,tV,N*tV);

        fprintf('                                                                ----------\n');
        fprintf('                                                        TOTAL = %8.2e s\n',N*(tK+tV));

        % Clean up
        QM.wfcn         =   inWfcn;
        
    end

    fprintf('\n  (average times over %u iterations)\n\n',obj.nbIter);

end
    % Suppress the output object ==========================================
    QMol_doc.showFooter;
    clear obj

end
end
%% Memory profiling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_suite)
function m = getMemoryFootprint(N,opt)
%getMemoryFootprint obtain the memory footprint of an array with N elements
%   of the type opt. The memory footpring m is expressed in bites

    % Compute size
    switch lower(opt)
        case {'real','double'}
            m           =   8 * N;
        case {'imag','complex'}
            m           =   16 * N;
        otherwise
            warning('QMol:SE_profiler:arrayType',['Unknown array type: ' opt '. No memory evaluation performed']);
            m           =   0;
    end

end
function showMemoryFootprint(msg,mem,lvl)
%showMemoryFootprint displays the memory foortrint 
    
    % Initialization
    if nargin < 3,  lvl = 1;        end

    % Display memory footprint
    switch lvl
        case 1
            fprintf('  * %-60s ',msg);      QMol_SE_profiler.showMemoryFormatted(mem);
        case 2
            fprintf('    > %-58s ',msg);    QMol_SE_profiler.showMemoryFormatted(mem);
        otherwise
            fprintf('    > %-59s ',msg);    QMol_SE_profiler.showMemoryFormatted(mem);
            warning('QMol:SE_profiler:showLevel','Maximum show level is 2; Entry downgraded to level 2');
    end
    fprintf('\n');
end
end
methods (Static,Access=private)
function showMemoryFormatted(mem)
%showMemoryFormatted display the formatted memory footprint
    
    if mem == 0   ||   isnan(mem)   || isinf(mem)                           % Do not display
    elseif mem <= 1024,       fprintf('%6u  B',mem);
    elseif mem <= 1024^2,   fprintf('%6.1f KB',mem/1024);
    elseif mem <= 1024^3,   fprintf('%6.1f MB',mem/1024^2);
    elseif mem <= 1024^4,   fprintf('%6.1f GB',mem/1024^3);
    else,                   fprintf('%6.1f TB',mem/1024^4);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

