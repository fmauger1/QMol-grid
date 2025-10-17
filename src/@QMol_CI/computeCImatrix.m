function [CI, DX] = computeCImatrix(obj)
%computeCImatrix computes the CI matrix associted with the CI model.
%   Use computeCImatrix to compute the CI matrix associated with the CI
%   model. The calculation first precomputes the core Hamiltonian matrix
%   and two-electron 4-tensor, which are then used to build the CI matrix.
%
%   obj.computeCImatrix performs the CI-matrix calculation and stores the
%   result in the CImatrix member property. If obj.isDipole is true, the
%   dipole-coupling matrix is also calculated and stored in the DXmatrix
%   member property.
%
%   CI = obj.computeCImatrix returns a copy of the matrix
%
%   [CI,DX] = obj.computeCImatrix, when obj.isDipole is true, also returns 
%   a copy of the dipole-coupling matrix.

    % Initialization ======================================================
        % (Re)initialize everythings
        isCSB           =   obj.isSetConfigBasis;           obj.reset;
        
        if obj.disp,    QMol_doc.showHeader;                                %#ok<ALIGN> 
                        QMol_doc.showSection('Build the configuration-interaction (CI) model [Visentin 2025]');
                        obj.initialize(2);
                        fprintf('\n');
                        QMol_doc.showBibliography({'Visentin 2025'});
                        QMol_doc.showFunding;
                        QMol_doc.showFooter;
    
                        QMol_doc.showSection('CI-matrix calculation');
                        fprintf('                                       Calculation progress (%%)\n');
                        fprintf('  Calculation step          0       20       40       60       80      100\n'); % 30 steps + 0
                        fprintf('  ----------------------    ----------------------------------------------\n');

        else,           obj.initialize(0);                                  end
        obj.isSetConfigBasis =  isCSB;
    
        % Miscellaneous
        ISO             =   unique(abs(obj.CSB(:)));                        % list of used spatial orbitals
        N               =   numel(ISO);                                     % number of used spatial orbitals
        IISO            =   NaN(1,ISO(end));                                % map spatial-orbital index to core matrix ...
        IISO(ISO)       =   1:numel(ISO);                                   % and 4-tensor coefficients

    % Precalculate core integrals =========================================
        % Initialization
        if obj.disp,        fprintf('  Core integrals            ');        end
        ndp             =   N^2;
        idp             =   0;
        ikl             =   0;

        H               =   zeros(N,N);

        for k = 1:N, for l = k:N                                            %#ok<ALIGN>
            % Calculate component
            H(k,l)      =   obj.getCoreIntegral(ISO(k),ISO(l));             if l ~=k   % <k|h|l>
            H(l,k)      =   conj(H(k,l));           ikl     =   ikl + 1;    end        % <l|h|k> = <k|h|l>^*
            ikl         =   ikl + 1;
            
            % Update progress bar
            if obj.disp
                ind     =   round(46*ikl/ndp);                              % 46 progress segments
                while ind > idp,    fprintf('|');   idp     =   idp + 1;    end
            end
        end,end

        % Finalize
        H(abs(H) < obj.tol) =   0;

        if obj.disp
            while 46 > idp,         fprintf('|');   idp     =   idp + 1;    end
            fprintf('\n');
        end

    % Precalculate dipole-coupling elements ===============================
        % Initialization
        if obj.isDip
            if obj.disp,        fprintf('  Dipole coupling elements  ');        end
            ndp             =   N^2;
            idp             =   0;
            ikl             =   0;
    
            D               =   zeros(N,N);
    
            for k = 1:N, for l = k:N                                            %#ok<ALIGN>
                % Calculate component
                D(k,l)      =   obj.getDipoleIntegral(ISO(k),ISO(l));           if l ~=k   % <k|r|l>
                D(l,k)      =   conj(D(k,l));           ikl     =   ikl + 1;    end        % <l|r|k> = <k|r|l>^*
                ikl         =   ikl + 1;
                
                % Update progress bar
                if obj.disp
                    ind     =   round(46*ikl/ndp);                              % 46 progress segments
                    while ind > idp,    fprintf('|');   idp     =   idp + 1;    end
                end
            end,end
    
            % Finalize
            H(abs(H) < obj.tol) =   0;
    
            if obj.disp
                while 46 > idp,         fprintf('|');   idp     =   idp + 1;    end
                fprintf('\n');
            end
        end
    
    % Precalculate the 4-tensor ===========================================
        % Initialization
        if obj.disp,        fprintf('  Two-electron integrals    ');        end
        ndp             =   N^4;
        idp             =   0;
        iklmn           =   0;

        ERI             =   NaN(N,N,N,N);
        for k = 1:N, for l=1:N, for m = 1:N, for n = 1:N                    %#ok<ALIGN>
            % Check if already calculated
            if ~isnan(ERI(k,l,m,n)),    continue,   end                     % move to the next iteration

            % Calculate component
            ERI(k,l,m,n)=   obj.getTwoElectronIntegral(ISO(k),ISO(l),ISO(m),ISO(n));    if isnan(ERI(l,k,n,m)) && ( l ~= k || n ~= m )  % <kl|mn> 
            ERI(l,k,n,m)=   ERI(k,l,m,n);           iklmn   =   iklmn + 1;  end,        if isnan(ERI(m,n,k,l)) && ( m ~= k || n ~= l )  % <lk|nm> = <kl|mn>
            ERI(m,n,k,l)=   conj(ERI(k,l,m,n));     iklmn   =   iklmn + 1;  end,        if isnan(ERI(n,m,l,k)) && ( n ~= k || m ~= l )  % <mn|kl> = <kl|mn>^*
            ERI(n,m,l,k)=   conj(ERI(k,l,m,n));     iklmn   =   iklmn + 1;  end                                                         % <nm|lk> = <kl|mn>^*
            iklmn       =   iklmn + 1;
            
            % Update progress bar
            if obj.disp
                ind     =   round(46*iklmn/ndp);                            % 46 progress segments
                while ind > idp,    fprintf('|');   idp     =   idp + 1;    end
            end

        end, end, end, end

        % Finalize
        ERI(abs(ERI) < obj.tol) =   0;

        if obj.disp
            while 46 > idp,         fprintf('|');   idp     =   idp + 1;    end
            fprintf('\n');
        end

    % Calculate the CI matrix =============================================
        % Initialization
        if obj.disp,        fprintf('  CI matrix                 ');        end

        N               =   size(obj.CSB,1);
        ndp             =   N^2;
        idp             =   0;
        ikl             =   0;

        CI              =   zeros(N,N);                                     if obj.isDip
        DX              =   zeros(N,N);                                     end         % dipole-coupling matrix

        % Calculate the CI matrix
        for k = 1:N
            % map indexes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            iOk         =   IISO(abs(obj.CSB(k,:)));
            iOk         =   iOk .* sign(obj.CSB(k,:));                      % spin-orbital indexes in the active space

            % Diagonal elements
            for m = iOk                                                                 %#ok<ALIGN>
                CI(k,k) =   CI(k,k) + H(abs(m),abs(m));         if obj.isDip            % dipole-coupling matrix
                DX(k,k) =   DX(k,k) + D(abs(m),abs(m));         end,        for n = iOk %#ok<ALIGN> <m|h|m>
                CI(k,k) =   CI(k,k) + 0.5*ERI(abs(m),abs(n),abs(m),abs(n)); if m*n >0   % <mn|mn> (spins always match)
                CI(k,k) =   CI(k,k) - 0.5*ERI(abs(m),abs(n),abs(n),abs(m)); end         % <mn|nm> (only for matching spins)
            end, end

            % Update progress bar
            ikl         =   ikl + 1;

            if obj.disp
                ind     =   round(46*ikl/ndp);                             % 46 progress segments
                while ind > idp,    fprintf('|');   idp     =   idp + 1;    end
            end

            % Off-diagonal elements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for l = k+1:N
                % Relative excitation indexes
                [a,r,s] =   obj.getExcitationIndexes(obj.CSB(k,:),obj.CSB(l,:));
                ia      =   IISO(abs(a));                                   % spatial-orbital indexes in the active space
                ir      =   IISO(abs(r));

                if isscalar(a),                                                             if a*r > 0  % in all cases, a and r must have matching spin
                    CI(k,l) =   H(ia,ir);                       if obj.isDip                            % dipole-coupling matrix
                    DX(k,l) =   D(ia,ir);                       end,                for n = iOk         % <an|rn> (spins always match)
                    CI(k,l) =   CI(k,l) + ERI(ia,abs(n),ir,abs(n));         if a*n > 0                  % <an|nr> (only for matching spins)
                    CI(k,l) =   CI(k,l) - ERI(ia,abs(n),abs(n),ir);         end,    end,    end
    
                elseif numel(a) == 2,                                       if all(a.*r > 0)
                    CI(k,l) =   ERI(ia(1),ia(2),ir(1),ir(2));               end,    if all(a.*r([2 1]) > 0) % < a b | r s >
                    CI(k,l) =   CI(k,l)-ERI(ia(1),ia(2),ir(2),ir(1));       end                             % < a b | s r >   
                    
                end

                % Signature and symmetry
                if ~isscalar(s) || (abs(CI(k,l)) < obj.tol),    CI(k,l) =   0;  else    % remove round-off errors
                    CI(k,l) =   s*CI(k,l);                                              % permutation signature
                end
                if obj.isDip,                                                   if (~isscalar(s) || (abs(DX(k,l)) < obj.tol))
                    DX(k,l) =   0;                                              else    % remove round-off errors
                    DX(k,l) =   s*DX(k,l);                                      end     % permutation signature
                end
    
                CI(l,k)     =   conj(CI(k,l));                  if obj.isDip            % dipole-coupling matrix
                DX(l,k)     =   conj(DX(k,l));                  end
                ikl         =   ikl + 2;
    
                % Update progress bar
                if obj.disp
                    ind     =   round(46*ikl/ndp);                             % 46 progress segments
                    while ind > idp,    fprintf('|');   idp     =   idp + 1;    end
                end
            end
        end

        % Finalize
        if obj.disp
            while 46 > idp,         fprintf('|');   idp     =   idp + 1;    end
            fprintf('\n');
        end

    % Finalization ========================================================
        % Copy CI matrix
        obj.CI          =   CI;                                             if obj.isDip
        obj.DX          =   DX;                                             end

        % Footer
        if obj.disp
            fprintf('  ------------------------------------------------------------------------\n');
            fprintf('                                            CI-matrix calculation complete\n\n');
            QMol_doc.showFooter;                
        end
end

