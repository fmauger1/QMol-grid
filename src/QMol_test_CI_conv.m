classdef QMol_test_CI_conv < QMol_test
%QMol_test_CI_conv suite of unit tests for QMol_CI_conv

%   Version     Date        Author
%   01.22.000   05/20/2025  F. Mauger
%       Creation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_test)
function version
    QMol_doc.showVersion('01.22.000','05/20/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test_CI_conv\n'); 
    QMol_test_CI_conv.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class

    % CI matrix calculation
    obj.test_computeCImatrix
    
end
end
%% Test components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=private)
    function test_computeCImatrix(obj)

%% Initialization =========================================================
    obj.showSection('Compute the CI matrix (computeCImatrix)');

%% HF model ===============================================================
    % Molecular model (artificially increased bare charges to help converge the HF model)
    V_1                 =   QMol_Va_softCoulomb('atom','X_1','charge',5,'position',-2,'softeningParameter',1);
    V_2                 =   QMol_Va_softCoulomb('atom','X_2','charge',3,'position', 2,'softeningParameter',1);
    
    % Spin-restricted HF model
    Vext                =   QMol_DFT_Vext('atom',{V_1,V_2});
    Vh                  =   QMol_DFT_Vh_conv('interactionPotential',@(x) 1./sqrt(x.^2+1.7));  % increase strength of electron-electron repulsion
    Vxc                 =   QMol_DFT_Vx_XX_conv('interactionPotential',@(x) 1./sqrt(x.^2+1.7));
    
    HF                  =   QMol_DFT_spinRes(                                   ...
                                'xspan',                       -22:.2:22,       ...
                                'occupation',                  [2 2 2 0 0 0],   ...
                                'externalPotential',            Vext,           ...
                                'HartreePotential',             Vh,             ...
                                'exchangeCorrelationPotential', Vxc             );
    
    % Compute HF ground state
    SCF                 =   QMol_DFT_SCF_Anderson('display',false);
    SCF.solveSCF(HF);

%% CISD model =============================================================
    % Create CI model
    CI                  =   QMol_CI_conv(HF,'display',false,'tolerance',1e-7,'type','CIS');
    CI.setConfigurationBasis;

    % Compute the CI matrix
    CI.computeCImatrix;

%% Test CI matrix =========================================================
    % Initialization
    dx                  =   HF.discretization.dx;                           % discretization step
    T                   =   HF.discretization.T;                            % kinetic operator (in Fourier)
    V                   =   HF.externalPotential.getPotential;              % external (electron-nuclear) potential
    V                   =   V.potential;
    V2e                 =   HF.HartreePotential.interactionPotential(   ...   electron-electron interaction kernel
                                (0:2*numel(HF.xspan)-2).' * dx - (HF.xspan(end)-HF.xspan(1)) );

    CSB                 =   CI.configurationBasis;

    Ref                 =   CSB(1,:);
    chi                 =   CI.orbitalBasis;

    function h = getCoreIntegral(chi_k,chi_l)
    %calculates and returns the core integral < k | h | l >
        
        % Kinetic component
        h   =   sum( conj(chi_k) .* ifft(T .* fft(chi_l)) );
    
        % Potential component
        h   =   h + sum(conj(chi_k) .* V .* chi_l);
    
        % Finalize
        h   =   h * dx;
    end

    function v = getTwoElectronIntegral(chi_k,chi_l,chi_m,chi_n)
        %calculates and returns the two-electron integral < k l | m n>
        v   =   sum(conj(chi_k).*chi_m .* conv(conj(chi_l).*chi_n,V2e,'same')) * dx^2;
    end

    % Check HF reference energy
    E_HF                =   HF.getEnergy('DFT');
    obj.showResult('HF reference energy',abs(E_HF-CI.CImatrix(1,1)) < CI.tolerance);

    % Check Brillioun's theorem for HF orbitals
    No                  =   sum(HF.occupation > 0);                         % number of occupied spatial orbitals
    Nv                  =   numel(HF.occupation)-No;                        % number of virtual orbitals
    N_CIS               =   2*No*Nv;

    obj.showResult('Brillioun''s theorem for HF orbitals',all(CI.CImatrix(1,2:N_CIS+1) == 0));

    % Check single-excitation diagonal components
    for n = 1:N_CIS
        % Identify excitation
        C               =   intersect(Ref,CSB(n+1,:));
        a               =   setdiff(Ref,C);
        r               =   setdiff(CSB(n+1,:),C);
    
        % Calculate CI matrix element
        E               =   CI.CImatrix(1,1) + ...
                                getCoreIntegral(chi(:,abs(r)),chi(:,abs(r))) - ...
                                getCoreIntegral(chi(:,abs(a)),chi(:,abs(a))) + ...
                                getTwoElectronIntegral(chi(:,abs(r)),chi(:,abs(a)),chi(:,abs(a)),chi(:,abs(r))) - ...
                                getTwoElectronIntegral(chi(:,abs(r)),chi(:,abs(a)),chi(:,abs(r)),chi(:,abs(a)));
        for k = Ref
            E           =   E + getTwoElectronIntegral(chi(:,abs(r)),chi(:,abs(k)),chi(:,abs(r)),chi(:,abs(k))) - ...
                                getTwoElectronIntegral(chi(:,abs(r)),chi(:,abs(k)),chi(:,abs(k)),chi(:,abs(r))) * (r*k > 0) - ...
                                getTwoElectronIntegral(chi(:,abs(a)),chi(:,abs(k)),chi(:,abs(a)),chi(:,abs(k))) + ...
                                getTwoElectronIntegral(chi(:,abs(a)),chi(:,abs(k)),chi(:,abs(k)),chi(:,abs(a))) * (a*k > 0);
        end
    
        % Check result
        isOK            =   abs(CI.CImatrix(n+1,n+1)-E) < CI.tolerance;
        if ~isOK,           break;                                          end
    end

    obj.showResult('Single-excitation (diagonal elements) for HF orbitals',isOK);

    % Check single-excitation off-diagonal components
    for n = 1:N_CIS, for m = n+1:N_CIS                                                      %#ok<ALIGN>
        % Identify excitation
        C               =   intersect(Ref,CSB(n+1,:));
        a               =   setdiff(Ref,C);
        r               =   setdiff(CSB(n+1,:),C);
    
        C               =   intersect(Ref,CSB(m+1,:));
        b               =   setdiff(Ref,C);
        s               =   setdiff(CSB(m+1,:),C);
    
        % Calculate CI matrix element
        E               =   getTwoElectronIntegral(chi(:,abs(r)),chi(:,abs(b)),chi(:,abs(a)),chi(:,abs(s))) - ...
                            getTwoElectronIntegral(chi(:,abs(r)),chi(:,abs(b)),chi(:,abs(s)),chi(:,abs(a))) * (a*b > 0);
        if (a == b) && (r*s > 0)                                                            %#ok<ALIGN>
            E           =   E + getCoreIntegral(chi(:,abs(r)),chi(:,abs(s)));    for k = Ref %#ok<ALIGN>
            E           =   E + getTwoElectronIntegral(chi(:,abs(r)),chi(:,abs(k)),chi(:,abs(s)),chi(:,abs(k))) - ...
                                getTwoElectronIntegral(chi(:,abs(r)),chi(:,abs(k)),chi(:,abs(k)),chi(:,abs(s))) * (r*k > 0);
        end, end
        if (r == s) && (a*b > 0)                                                            %#ok<ALIGN>
            E           =   E - getCoreIntegral(chi(:,abs(a)),chi(:,abs(b)));    for k = Ref %#ok<ALIGN>
            E           =   E - getTwoElectronIntegral(chi(:,abs(a)),chi(:,abs(k)),chi(:,abs(b)),chi(:,abs(k))) + ...
                                getTwoElectronIntegral(chi(:,abs(a)),chi(:,abs(k)),chi(:,abs(k)),chi(:,abs(b))) * (a*k > 0);
        end, end
    
        % Check result
        isOK            =   abs(CI.CImatrix(n+1,m+1)-E) < CI.tolerance;
        if ~isOK,           break;                                          end
    
    end, end

    obj.showResult('Single-excitation (off-diagonal elements) for HF orbitals',isOK);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

