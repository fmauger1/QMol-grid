classdef QMol_TDCI_SSO_4BM < QMol_TDCI_SSO
%QMol_TDDFT_SSO_4BM 4th order Blanes and Moan optimized symplectic split-
%   operator scheme

%   Version     Date        Author
%   01.23.000   06/07/2025  F. Mauger
%       Creation (from QMol_TDSE_SSO_4BM)

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.000','06/07/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDCI})
function showInfo
    fprintf('  * QMol_TDCI_SSO_4BM:\n');
    fprintf('      > Symplectic split-operator TDCI propagator\n'); 
    fprintf('      > Optimized Blanes and Moan scheme (4th order)\n'); 
    QMol_TDCI_SSO_4BM.version;
end
end
methods (Access=protected)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the TDCI solver

    % Call parent documentation
    ref                 =   showDoc@QMol_TDCI_sympSplitOp(obj);
    
    % Implementation of the TDDFT model
    fprintf('  * Optimized Blanes and Moan scheme [Blanes 2002]               4th order\n');
    obj.version;

    ref                 =   [ref,{'Blanes 2002'}];

end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    % Parent-class components
    [~,PropNames]       =   QMol_TDCI_sympSplitOp.propertyNames;

    ClassName           =   'QMol_TDCI_SSO_4BM';
end
end
%% Time propagation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=protected)
function setTimeStep(obj,dt,t)
%setTimeStep sets the time step, and performs associated initializations,
%   for the Strang splitting propagator
    
    % Copy time step
    obj.dt_             =   dt;

    % Initialize split-operator components
    obj.t_              =   t;

end
function applyTimeStep(obj,t)
%applyTimeStep applies one iteration of the Strang splitting operator
    
    % Update t
    obj.t_              =   t;

    % Time propagation
    if obj.isHDH,       obj.applyExpH(QMol_cte.symp_4BM2_sC1);
                        obj.applyExpD(QMol_cte.symp_4BM2_sD1);
                        obj.applyExpH(QMol_cte.symp_4BM2_sC2);
                        obj.applyExpD(QMol_cte.symp_4BM2_sD2);
                        obj.applyExpH(QMol_cte.symp_4BM2_sC3);
                        obj.applyExpD(QMol_cte.symp_4BM2_sD3);
                        obj.applyExpH(QMol_cte.symp_4BM2_sC4);
                        obj.applyExpD(QMol_cte.symp_4BM2_sD4);
                        obj.applyExpH(QMol_cte.symp_4BM2_sC5);
                        obj.applyExpD(QMol_cte.symp_4BM2_sD5);
                        obj.applyExpH(QMol_cte.symp_4BM2_sC6);
                        obj.applyExpD(QMol_cte.symp_4BM2_sD6);
                        obj.applyExpH(QMol_cte.symp_4BM2_sC7);
    else,               obj.applyExpD(QMol_cte.symp_4BM2_sC1);
                        obj.applyExpH(QMol_cte.symp_4BM2_sD1);
                        obj.applyExpD(QMol_cte.symp_4BM2_sC2);
                        obj.applyExpH(QMol_cte.symp_4BM2_sD2);
                        obj.applyExpD(QMol_cte.symp_4BM2_sC3);
                        obj.applyExpH(QMol_cte.symp_4BM2_sD3);
                        obj.applyExpD(QMol_cte.symp_4BM2_sC4);
                        obj.applyExpH(QMol_cte.symp_4BM2_sD4);
                        obj.applyExpD(QMol_cte.symp_4BM2_sC5);
                        obj.applyExpH(QMol_cte.symp_4BM2_sD5);
                        obj.applyExpD(QMol_cte.symp_4BM2_sC6);
                        obj.applyExpH(QMol_cte.symp_4BM2_sD6);
                        obj.applyExpD(QMol_cte.symp_4BM2_sC7);
    end

    % Mask damping
    if ~isempty(obj.D),     obj.CI.wfcn     =   obj.D*obj.CI.wfcn;          end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

