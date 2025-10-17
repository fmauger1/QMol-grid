classdef QMol_DFT_profiler < QMol_profiler
%QMol_DFT_profiler time/memory profiler for DFT systems
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release
%   01.23.000   05/25/2025  F. Mauger
%       Derive from QMol_profiler

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.000','05/25/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_profiler})
function showInfo
    fprintf('  * QMol_DFT_profiler:\n');
    fprintf('      > Memory/execution time profiling for DFT\n'); 
    QMol_DFT_profiler.version;
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties
    mode                =   'all'
    nbIter              =   []
end

%% Run profiling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function obj = QMol_DFT_profiler(QM,opt,N)
%QMol_DFT_profiler performs the profiling for the input DFT component
    
    % DFT profiler header =================================================
    fprintf(['########################### QMol-grid package ############################\n'...
             'Quantum simulation methods for atomic and molecular systems [Mauger 2024b]\n\n' ...
             '  * DFT profiler (release)\n']);
    QMol_DFT_profiler.version;      fprintf('\n');
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
    QMol_profiler.showMemoryFormatted(mem);
    fprintf('\n\n');

end
    % Execution time ======================================================
if any(strcmpi(obj.mode,{'all','time','execution'}))
    QMol_doc.showSection('Execution time');

    % Initialization
    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility

    if isa(QM,'QMol_TDDFT') % TDDFT object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fprintf('*** No execution-time estimate for TDDFT objects.                      ***\n');
    elseif isa(QM,'QMol_DFT') % DFT object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Initialization
        QM.initialize;      inKSO   =   QM.KSO;

        if isempty(obj.nbIter)
            if     QM.dim == 1, obj.nbIter  =   500;
            elseif QM.dim == 2, obj.nbIter  =   50;
            elseif QM.dim == 3, obj.nbIter  =   10;
            else,  error('QMol:QMol_DFT_profiler:dimDFT',['Unexpected dimension ' num2str(QM.dim) ' (1, 2, or 3 allowed)']);
            end
        end


        % Spin polarized  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if QM.isSpinPol
            % Initialization
            N           =   [numel(QM.occ{1}) numel(QM.occ{2})];
            QM.KSO      =   QM.disc.DFT_allocateOrbital(N,[],randStr,false);
            p           =   QM.disc.DFT_randomOrbital(randStr) + 1i*QM.disc.DFT_randomOrbital(randStr);
            p           =   QM.disc.DFT_normalizeOrbital(p);

            QM.setPotentialKernel;
            Vks         =   QM.getPotential();

            fprintf('  Operator components                                           -- Time --\n');

            % Whole Hamiltonian operator
            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorHamiltonian(Vks,p,true,0);  end, tUp = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Hamiltonian (spin up)'  ,N(1),tUp,N(1)*tUp);

            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorHamiltonian(Vks,p,false,0);  end, tDw = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Hamiltonian (spin down)',N(2),tDw,N(2)*tDw);

            fprintf('                                                                ----------\n');
            fprintf('                                                        TOTAL = %8.2e s\n\n',N(1)*tUp+N(2)*tDw);

            % Kinetic & potential operators
            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorKinetic(p,0);           end, tK = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Kinetic',N(1)+N(2),tK,sum(N)*tK);

            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,true,0);    end, tVup = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Potential (spin up)',N(1),tVup,N(1)*tVup);

            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,false,0);   end, tVdw = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Potential (spin down)',N(2),tVdw,N(2)*tVdw);

            fprintf('                                                                ----------\n');
            fprintf('                                                        TOTAL = %8.2e s\n\n',sum(N)*tK+N(1)*tVup+N(2)*tVdw);

            % DFT-functional operators
            QM.Vext.getPotential(Vks);
            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,true,0);    end, tExt = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','External potential',N(1)+N(2),tExt,sum(N)*tExt);

            QM.Vh.getPotential([],Vks);
            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,true,0);    end, tH = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Hartree potential',N(1)+N(2),tH,sum(N)*tH);

            if iscell(QM.Vxc), for l = 1:numel(QM.Vxc) %#ok<ALIGN> 
                QM.Vxc{l}.getPotential([],Vks);
                tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,true,0);    end, txc = toc/obj.nbIter; %#ok<NASGU> 
                fprintf('  * %-42s %3u x %8.2e = %8.2e s\n',['Exch.-corr. potential #' num2str(l) ' (spin up)'],N(1),txc,N(1)*txc);
                tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,true,0);    end, txc = toc/obj.nbIter; %#ok<NASGU> 
                fprintf('  * %-42s %3u x %8.2e = %8.2e s\n',['Exch.-corr. potential #' num2str(l) ' (spin down)'],N(2),txc,N(2)*txc);

            end, else
                QM.Vxc.getPotential([],Vks);
                tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,true,0);    end, txc = toc/obj.nbIter; %#ok<NASGU> 
                fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Exchange-correlation potential (spin up)',N(1),txc,N(1)*txc);
                tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,false,0);    end, txc = toc/obj.nbIter; %#ok<NASGU> 
                fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Exchange-correlation potential (spin down)',N(2),txc,N(2)*txc);
            end

            % Building potentials
            fprintf('\n  Building potential components                                 -- Time --\n');
            rho         =   QM.getDensity;
            
            tic, for k = 1:obj.nbIter,  QM.getDensity(rho);                     end, t = toc/obj.nbIter;
            fprintf('  * %-59s %8.2e s\n','One-body density',t);
            
            tic, for k = 1:obj.nbIter,  QM.Vext.getPotential(Vks);              end, t = toc/obj.nbIter;
            fprintf('  * %-59s %8.2e s\n','External functional',t);
            
            tic, for k = 1:obj.nbIter,  QM.Vh.getPotential(rho,Vks);            end, t = toc/obj.nbIter;
            fprintf('  * %-59s %8.2e s\n','Hartree functional',t);

            if iscell(QM.Vxc), for l = 1:numel(QM.Vxc) %#ok<ALIGN> 
                tic, for k = 1:obj.nbIter,  QM.Vxc{l}.getPotential(rho,Vks);    end, t = toc/obj.nbIter;
                fprintf('  * %-59s %8.2e s\n',['Exchange-correlation functional #' num2str(l)],t);

            end, else
                tic, for k = 1:obj.nbIter,  QM.Vxc.getPotential(rho,Vks);       end, t = toc/obj.nbIter;
                fprintf('  * %-59s %8.2e s\n','Exchange-correlation functional',t);
            end

        % Spin restricted ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else
            % Initialization
            N           =   numel(QM.occ);
            QM.KSO      =   QM.disc.DFT_allocateOrbital(N,[],randStr,false);
            p           =   QM.disc.DFT_randomOrbital(randStr) + 1i*QM.disc.DFT_randomOrbital(randStr);
            p           =   QM.disc.DFT_normalizeOrbital(p);

            QM.setPotentialKernel;
            Vks         =   QM.getPotential();

            fprintf('  Operator components                                           -- Time --\n');

            % Whole Hamiltonian operator
            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorHamiltonian(Vks,p,0);  end, t = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Hamiltonian'  ,N,t,N*t);

            fprintf('\n');

            % Kinetic & potential operators
            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorKinetic(p,0);           end, tK = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Kinetic',N,tK,N*tK);

            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,0);     end, tV = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Potential',N,tV,N*tV);

            fprintf('                                                                ----------\n');
            fprintf('                                                        TOTAL = %8.2e s\n\n',N*(tK+tV));

            % DFT-functional operators
            QM.Vext.getPotential(Vks);
            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,0);     end, tExt = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','External potential',N,tExt,N*tExt);

            QM.Vh.getPotential([],Vks);
            tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,0);     end, tH = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Hartree potential',N,tH,N*tH);

            if iscell(QM.Vxc), for l = 1:numel(QM.Vxc) %#ok<ALIGN> 
                QM.Vxc{l}.getPotential([],Vks);
                tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,0); end, txc = toc/obj.nbIter; %#ok<NASGU> 
                fprintf('  * %-42s %3u x %8.2e = %8.2e s\n',['Exchange-correlation potential #' num2str(l)],N,txc,N*txc);

            end, else
                QM.Vxc.getPotential([],Vks);
                tic, for k = 1:obj.nbIter,  Hp  =   QM.disc.DFT_operatorPotential(Vks,p,0); end, txc = toc/obj.nbIter; %#ok<NASGU> 
                fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Exchange-correlation potential',N,txc,N*txc);
            end

            % Building potentials
            fprintf('\n  Building potential components                                 -- Time --\n');
            rho         =   QM.getDensity;
            
            tic, for k = 1:obj.nbIter,  QM.getDensity(rho);                     end, t = toc/obj.nbIter;
            fprintf('  * %-59s %8.2e s\n','One-body density',t);
            
            tic, for k = 1:obj.nbIter,  QM.Vext.getPotential(Vks);              end, t = toc/obj.nbIter;
            fprintf('  * %-59s %8.2e s\n','External functional',t);
            
            tic, for k = 1:obj.nbIter,  QM.Vh.getPotential(rho,Vks);            end, t = toc/obj.nbIter;
            fprintf('  * %-59s %8.2e s\n','Hartree functional',t);

            if iscell(QM.Vxc), for l = 1:numel(QM.Vxc) %#ok<ALIGN> 
                tic, for k = 1:obj.nbIter,  QM.Vxc{l}.getPotential(rho,Vks);    end, t = toc/obj.nbIter;
                fprintf('  * %-59s %8.2e s\n',['Exchange-correlation functional #' num2str(l)],t);

            end, else
                tic, for k = 1:obj.nbIter,  QM.Vxc.getPotential(rho,Vks);       end, t = toc/obj.nbIter;
                fprintf('  * %-59s %8.2e s\n','Exchange-correlation functional',t);
            end

        end

        % Clean up
        QM.KSO          =   inKSO;
        
    elseif isa(QM,'QMol_Vmol') % External functional object ~~~~~~~~~~~~~~~
        % Initialization
        p           =   QM.DFT.disc.DFT_randomOrbital(randStr) + 1i*QM.DFT.disc.DFT_randomOrbital(randStr);
        p           =   QM.DFT.disc.DFT_normalizeOrbital(p);

        if QM.DFT.isSpinPol,    N   =  [numel(QM.DFT.occ{1}) numel(QM.DFT.occ{2})];
        else,                   N   =   numel(QM.DFT.occ);                          end

        if isempty(obj.nbIter)
            if     QM.DFT.dim == 1, obj.nbIter  =   500;
            elseif QM.DFT.dim == 2, obj.nbIter  =   50;
            elseif QM.DFT.dim == 3, obj.nbIter  =   10;
            else,  error('QMol:QMol_DFT_profiler:dimDFT',['Unexpected dimension ' num2str(QM.dim) ' (1, 2, or 3 allowed)']);
            end
        end

        % Time profiling
        Vks             =   QM.getPotential();

        fprintf('  Components                                                    -- Time --\n');

        if QM.DFT.isSpinPol, tic, for k = 1:obj.nbIter,  Hp  =   QM.DFT.disc.DFT_operatorPotential(Vks,p,true,0);    end        %#ok<NASGU> 
        else,                tic, for k = 1:obj.nbIter,  Hp  =   QM.DFT.disc.DFT_operatorPotential(Vks,p,     0);    end, end   %#ok<NASGU> 
        t               =   toc/obj.nbIter;
        fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','External potential operator',sum(N),t,sum(N)*t);

        tic, for k = 1:obj.nbIter,  QM.getPotential(Vks);              end, t = toc/obj.nbIter;
        fprintf('  * %-59s %8.2e s\n','Build external potential',t);

        
    elseif ~isprop(QM,'type') % Hartree functional object ~~~~~~~~~~~~~~~~~
        % Initialization
        QM.DFT.initialize;      inKSO       =   QM.DFT.KSO;

        if QM.DFT.isSpinPol,    N   =  [numel(QM.DFT.occ{1}) numel(QM.DFT.occ{2})];
        else,                   N   =   numel(QM.DFT.occ);                          end

        QM.DFT.KSO      =   QM.DFT.disc.DFT_allocateOrbital(N,[],randStr,false);
        QM.DFT.setPotentialKernel;

        p           =   QM.DFT.disc.DFT_randomOrbital(randStr) + 1i*QM.DFT.disc.DFT_randomOrbital(randStr);
        p           =   QM.DFT.disc.DFT_normalizeOrbital(p);

        if isempty(obj.nbIter)
            if     QM.DFT.dim == 1, obj.nbIter  =   500;
            elseif QM.DFT.dim == 2, obj.nbIter  =   50;
            elseif QM.DFT.dim == 3, obj.nbIter  =   10;
            else,  error('QMol:QMol_DFT_profiler:dimDFT',['Unexpected dimension ' num2str(QM.dim) ' (1, 2, or 3 allowed)']);
            end
        end

        % Time profiling
        Vks             =   QM.getPotential();

        fprintf('  Components                                                    -- Time --\n');

        if QM.DFT.isSpinPol, tic, for k = 1:obj.nbIter,  Hp  =   QM.DFT.disc.DFT_operatorPotential(Vks,p,true,0);    end        %#ok<NASGU> 
        else,                tic, for k = 1:obj.nbIter,  Hp  =   QM.DFT.disc.DFT_operatorPotential(Vks,p,     0);    end, end   %#ok<NASGU> 
        t               =   toc/obj.nbIter;
        fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Hartree potential operator',sum(N),t,sum(N)*t);

        rho             =   QM.DFT.getDensity();
        tic, for k = 1:obj.nbIter,  QM.getPotential(rho,Vks);               end, t = toc/obj.nbIter;
        fprintf('  * %-59s %8.2e s\n','Build Hartree potential',t);

        % Clean up
        QM.DFT.KSO      =   inKSO;

    else % Exchange-correlation functional object ~~~~~~~~~~~~~~~~~~~~~~~~~
        % Initialization
        QM.DFT.initialize;      inKSO       =   QM.DFT.KSO;

        if QM.DFT.isSpinPol,    N   =  [numel(QM.DFT.occ{1}) numel(QM.DFT.occ{2})];
        else,                   N   =   numel(QM.DFT.occ);                          end

        QM.DFT.KSO      =   QM.DFT.disc.DFT_allocateOrbital(N,[],randStr,false);
        QM.DFT.setPotentialKernel;

        p           =   QM.DFT.disc.DFT_randomOrbital(randStr) + 1i*QM.DFT.disc.DFT_randomOrbital(randStr);
        p           =   QM.DFT.disc.DFT_normalizeOrbital(p);

        if isempty(obj.nbIter)
            if     QM.DFT.dim == 1, obj.nbIter  =   500;
            elseif QM.DFT.dim == 2, obj.nbIter  =   50;
            elseif QM.DFT.dim == 3, obj.nbIter  =   10;
            else,  error('QMol:QMol_DFT_profiler:dimDFT',['Unexpected dimension ' num2str(QM.dim) ' (1, 2, or 3 allowed)']);
            end
        end

        % Time profiling
        Vks             =   QM.getPotential();

        fprintf('  Components                                                    -- Time --\n');

        if QM.DFT.isSpinPol
            % Spin up
            tic, for k = 1:obj.nbIter,  Hp  =   QM.DFT.disc.DFT_operatorPotential(Vks,p,true,0);    end,    t = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Exch.-corr. potential operator (spin up)',N(1),t,N(1)*t);

            % Spin down
            tic, for k = 1:obj.nbIter,  Hp  =   QM.DFT.disc.DFT_operatorPotential(Vks,p,false,0);   end,    t = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Exch.-corr. potential operator (spin down)',N(2),t,N(2)*t);
        else
            tic, for k = 1:obj.nbIter,  Hp  =   QM.DFT.disc.DFT_operatorPotential(Vks,p,     0);    end,    t = toc/obj.nbIter; %#ok<NASGU> 
            fprintf('  * %-42s %3u x %8.2e = %8.2e s\n','Exch.-corr. potential operator',N,t,N*t);
        end   
        
        rho             =   QM.DFT.getDensity();
        tic, for k = 1:obj.nbIter,  QM.getPotential(rho,Vks);               end, t = toc/obj.nbIter;
        fprintf('  * %-59s %8.2e s\n','Build Exchange-correlation potential',t);

        % Clean up
        QM.DFT.KSO      =   inKSO;
    
    end

    fprintf('\n  (average times over %u iterations)\n\n',obj.nbIter);

end
    % Suppress the output object ==========================================
    QMol_doc.showFooter;
    clear obj

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

