function propagate(obj,opt)
%propagate propagates the TDDFT equations using the (child object) scheme.

    % Initialization ======================================================
    obj.reset;          obj.isRun   =   true;                               % House keeping
    if obj.disp,        QMol_doc.showHeader;                                end


    if isa(opt,'QMol_DFT')                                                  % TDDFT propagation from scratch
        % Initialize the DFT object
        if obj.disp,    QMol_doc.showSection('Build density-functional-theory (DFT) model'); %#ok<ALIGN> 
                        opt.initialize(2);
        else,           opt.initialize(0);                                  end

        % Initialize the propagator
        obj.initialize(opt,false);

        if obj.disp,    obj.showDocumentation;                              end

        % Progress bar and misc.
        kstart          =   1;
        ndisp           =   2;
        idp             =   0;

        if obj.disp
            % Progress header
            QMol_doc.showSection('Time propagation')
    
            fprintf('    Iteration time (a.u.)    Iteration         Iteration progress (%%)\n');
            fprintf('      start      finish        number      0    20    40    60    80   100\n'); % 30 steps + 0
            fprintf('    ---------   ---------    ---------     -------------------------------');

            % Track mark
            fprintf('\n    %9.3f   %9.3f    %9u     |',obj.tspan(1),obj.tspan(obj.iref(2)),1);
        end

        % Set time step and else
        obj.setTimeStep(obj.dt,obj.tspan(kstart));

    elseif  ischar(opt) && any(strcmpi(opt,{'restart','resume'}))           % Restart mode                                                 % TDDFT propagation from scratch
        % Initialize the DFT object
        if obj.disp,    QMol_doc.showSection('Build density-functional-theory (DFT) model'); %#ok<ALIGN> 
                        obj.DFT.initialize(2);
        else,           obj.DFT.initialize(0);                              end

        % Initialize the propagator
        obj.initialize([],true);

        if obj.disp,    obj.showDocumentation;                              end

        % Progress bar and misc.
        kstart          =   obj.restart.kstart;
        ndisp           =   obj.restart.ndisp;
        idp             =   obj.restart.idp;

        if obj.disp
            % Progress header
            QMol_doc.showSection('Time propagation')
    
            fprintf('    Iteration time (a.u.)    Iteration     Iteration progress (%%)\n');
            fprintf('      start      finish        number      0    20    40    60    80   100\n'); % 30 steps + 0
            fprintf('    ---------   ---------    ---------     -------------------------------');

            % Track mark
            fprintf('\n    %9.3f   %9.3f    %9u     *',obj.tspan(obj.iref(ndisp-1)),obj.tspan(obj.iref(ndisp)),obj.iref(ndisp-1));

            for k = 0:idp-1,    fprintf('*');                               end

        end

    else                                                                    % Unexpected case
        error('QMol:TDDFT:propagate','Unexpected TDDFT propagation case');
    end

    % Time propagation ====================================================
    for k = kstart:numel(obj.tspan)-2
        % New progress bar
        if k == obj.iref(ndisp)
            % Update display
            if obj.disp
                fprintf('\n    %9.3f   %9.3f    %9u     |',obj.tspan(k),obj.tspan(obj.iref(ndisp+1)),k);
                idp     =   0;
            end

            % Update display counter
            ndisp       =   ndisp + 1;
        end
        
        % Save data
        obj.saveOutputResults(k,obj.tspan(k));

        % Advance by one time step
        obj.applyTimeStep(obj.tspan(k));
        if ~isempty(obj.ABC),       obj.ABC.applyMask(obj.DFT.KSO);         end

        % Update progress bar
        if obj.disp
            ind         =   round(30*(k+1-obj.iref(ndisp-1))/(obj.iref(ndisp)-obj.iref(ndisp-1))); % 30 progress segments
            while ind > idp,    fprintf('|');   idp = idp + 1;  end
        end

        % Restart file
        if k == obj.oRest.ind(obj.oRest.n)
            % Update restart structure
            obj.restart.kstart      =   k+1;
            obj.restart.ndisp       =   ndisp;
            obj.restart.idp         =   idp;

            % Children class add needed info to restart structure
            obj.saveRestartChildren(k+1,obj.tspan(k+1));

            % Update Counter
            obj.oRest.n             =   obj.oRest.n + 1;
    
            % Save copy of the object
            TDDFT           =   obj;
            DFT             =   obj.DFT;
            save(obj.sRestF,'TDDFT','DFT');
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
        if ~isempty(obj.ABC),       obj.ABC.applyMask(obj.DFT.KSO);         end

        % Save data
        k               =   numel(obj.tspan);
        obj.saveOutputResults(k,obj.tspan(k));
        
        % Update progress bar
        if obj.disp,    while 30 > idp,    fprintf('|');   idp = idp + 1;  end, end % No matter what, we are done with the iteration
        
    % Finalization ========================================================
        % Cleanup output structures
        obj.finalize;
        obj.isRun       =   false;
        
        % Footer 
        if obj.disp
            fprintf('\n    ---------   ---------    ---------     -------------------------------\n');
            fprintf('    Time propagation finished without error. All results are saved in the\n');
            fprintf('    TDDFT object and specified output MATLAB file(s).\n\n');

            QMol_doc.showFooter;
        end
end