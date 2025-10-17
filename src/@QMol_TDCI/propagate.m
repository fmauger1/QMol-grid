function propagate(obj,opt)
%propagate propagates the TDCI equations using the (child object) scheme.

    % Initialization ======================================================
    obj.reset;          obj.isRun   =   true;                               % House keeping
    if obj.disp,        QMol_doc.showHeader;                                end


    if isa(opt,'QMol_CI')                                                   % TDCI propagation from scratch
        % Initialize the DFT object
        if obj.disp,    QMol_doc.showSection('Build the configuration interaction (CI) model'); %#ok<ALIGN> 
                        opt.initialize(2);
        else,           opt.initialize(0);                                  end

        % Initialize the propagator
        obj.initialize(opt,false);

        if obj.disp,    obj.showDocumentation;                              end

        % Progress bar and misc.
        kstart          =   1;
        ndisp           =   2;
        pBar            =   QMol_progressBar('N',31,'showZero',true);

        if obj.disp
            % Progress header
            QMol_doc.showSection('Time propagation')
    
            fprintf('    Iteration time (a.u.)    Iteration         Iteration progress (%%)\n');
            fprintf('      start      finish        number      0    20    40    60    80   100\n'); % 30 steps + 0
            fprintf('    ---------   ---------    ---------     -------------------------------\n');

            % Track mark
            fprintf('    %9.3f   %9.3f    %9u     ',obj.tspan(1),obj.tspan(obj.iref(2)),1);
            pBar.set('vMin',kstart,'vMax',obj.iref(2));
            pBar.start(kstart);
        end

        % Set time step and else
        obj.setTimeStep(obj.dt,obj.tspan(kstart));

    elseif  ischar(opt) && any(strcmpi(opt,{'restart','resume'}))           % Restart mode
        % Initialize the CI object
        if obj.disp,    QMol_doc.showSection('Build configuration interaction (CI) model'); %#ok<ALIGN> 
                        obj.CI.initialize(2);
        else,           obj.CI.initialize(0);                               end

        % Initialize the propagator
        obj.initialize([],true);

        if obj.disp,    obj.showDocumentation;                              end

        % Progress bar and misc.
        kstart          =   obj.restart.kstart;
        ndisp           =   obj.restart.ndisp;
        pBar            =   QMol_progressBar('N',31,'showZero',true,'motif','*');

        if obj.disp
            % Progress header
            QMol_doc.showSection('Time propagation')
    
            fprintf('    Iteration time (a.u.)    Iteration     Iteration progress (%%)\n');
            fprintf('      start      finish        number      0    20    40    60    80   100\n'); % 30 steps + 0
            fprintf('    ---------   ---------    ---------     -------------------------------\n');

            % Track mark
            fprintf('    %9.3f   %9.3f    %9u     ',obj.tspan(obj.iref(ndisp-1)),obj.tspan(obj.iref(ndisp)),obj.iref(ndisp-1));
            pBar.set('vMin',obj.iref(ndisp-1),'vMax',obj.iref(ndisp));
            pBar.start(kstart);
            pBar.set('motif','|')
            
        end

    else                                                                    % Unexpected case
        error('QMol:TDCI:propagate','Unexpected TDCI propagation case');
    end

    % Time propagation ====================================================
    for k = kstart:numel(obj.tspan)-2
        % New progress bar
        if k == obj.iref(ndisp)
            % Update display
            if obj.disp
                pBar.finish(true);
                fprintf('    %9.3f   %9.3f    %9u     ',obj.tspan(k),obj.tspan(obj.iref(ndisp+1)),k);
                pBar.set('vMin',k,'vMax',obj.iref(ndisp+1));
                pBar.start(k);
            end

            % Update display counter
            ndisp       =   ndisp + 1;
        end
        
        % Save data
        obj.saveOutputResults(k,obj.tspan(k));

        % Advance by one time step
        obj.applyTimeStep(obj.tspan(k));                                    % Unlike TDDFT/SE, must apply any mask damping

        % Update progress bar
        if obj.disp,                pBar.update(k);                         end

        % Restart file
        if k == obj.oRest.ind(obj.oRest.n)
            % Update restart structure
            obj.restart.kstart      =   k+1;
            obj.restart.ndisp       =   ndisp;

            % Children class add needed info to restart structure
            obj.saveRestartChildren(k+1,obj.tspan(k+1));

            % Update Counter
            obj.oRest.n             =   obj.oRest.n + 1;
    
            % Save copy of the object
            TDCI            =   obj;
            CI              =   obj.CI;
            save(obj.sRestF,'TDCI','CI');
        end
    end

    % Final time step =====================================================
        % Save data
        k               =   numel(obj.tspan)-1;                             % No restart in final step
        obj.saveOutputResults(k,obj.tspan(k));

        % Final step
        if abs(obj.tspan(end)-obj.tspan(end-1) - obj.dt) > 1e-10
            % Update time step in propagator
            obj.setTimeStep(obj.tspan(end)-obj.tspan(end-1),obj.tspan(end-1));
        end
        obj.applyTimeStep(obj.tspan(end-1));

        % Save data
        k               =   numel(obj.tspan);
        obj.saveOutputResults(k,obj.tspan(k));
        
        % Update progress bar
        if obj.disp,                pBar.finish(true);                      end
        
    % Finalization ========================================================
        % Cleanup output structures
        obj.finalize;
        obj.isRun       =   false;
        
        % Footer 
        if obj.disp
            fprintf('    ---------   ---------    ---------     -------------------------------\n');
            fprintf('    Time propagation finished without error. All results are saved in the\n');
            fprintf('    TDCI object and specified output MATLAB file(s).\n\n');

            QMol_doc.showFooter;
        end
end