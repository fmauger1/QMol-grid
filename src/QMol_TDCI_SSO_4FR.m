classdef QMol_TDCI_SSO_4FR < QMol_TDCI_SSO
%QMol_TDCI_SSO_4FR symplectic split-operator TDCI propagator using a
%   Forest-Ruth scheme (4th order)

%   Version     Date        Author
%   01.23.000   06/07/2025  F. Mauger
%       Creation (from QMol_TDSE_SSO_4FR)

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.000','06/07/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDCI})
function showInfo
    fprintf('  * QMol_TDCI_SSO_4FR:\n');
    fprintf('      > Symplectic split-operator TDDFT propagator\n'); 
    fprintf('      > Forest-Ruth scheme (4th order)\n'); 
    QMol_TDCI_SSO_4FR.version;
end
end
methods (Access=protected)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the TDCI solver

    % Call parent documentation
    ref                 =   showDoc@QMol_TDCI_sympSplitOp(obj);
    
    % Implementation of the TDDFT model
    fprintf('  * Forest-Ruth propagation scheme [Forest 1990]                 4th order\n');
    obj.version;

    ref                 =   [ref,{'Forest 1990'}];

end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    % Parent-class components
    [~,PropNames]       =   QMol_TDCI_sympSplitOp.propertyNames;

    ClassName           =   'QMol_TDCI_SSO_4FR';
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
    if obj.isHDH,       obj.applyExpH(QMol_cte.symp_4FR2_sC1);
                        obj.applyExpD(QMol_cte.symp_4FR2_sD1);
                        obj.applyExpH(QMol_cte.symp_4FR2_sC2);
                        obj.applyExpD(QMol_cte.symp_4FR2_sD2);
                        obj.applyExpH(QMol_cte.symp_4FR2_sC3);
                        obj.applyExpD(QMol_cte.symp_4FR2_sD3);
                        obj.applyExpH(QMol_cte.symp_4FR2_sC4);
    else,               obj.applyExpD(QMol_cte.symp_4FR2_sC1);
                        obj.applyExpH(QMol_cte.symp_4FR2_sD1);
                        obj.applyExpD(QMol_cte.symp_4FR2_sC2);
                        obj.applyExpH(QMol_cte.symp_4FR2_sD2);
                        obj.applyExpD(QMol_cte.symp_4FR2_sC3);
                        obj.applyExpH(QMol_cte.symp_4FR2_sD3);
                        obj.applyExpD(QMol_cte.symp_4FR2_sC4);
    end

    % Mask damping
    if ~isempty(obj.D),     obj.CI.wfcn     =   obj.D*obj.CI.wfcn;          end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

