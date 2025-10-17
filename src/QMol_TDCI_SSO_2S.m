classdef QMol_TDCI_SSO_2S < QMol_TDCI_SSO
%QMol_TDCI_SSO_2S symplectic split-operator TDDFT propagator using a
%   Strang-splitting scheme (2nd order)

%   Version     Date        Author
%   01.23.000   06/04/2025  F. Mauger
%       Creation (from QMol_TDSE_SSO_2S)
%   01.23.001   06/07/2025  F. Mauger
%       Use QMol_cte split coefficients (instead of hard coded)

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.001','06/07/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDCI})
function showInfo
    fprintf('  * QMol_TDCI_SSO_2S:\n');
    fprintf('      > Symplectic split-operator TDCI propagator\n'); 
    fprintf('      > Strang-splitting scheme (2nd order)\n'); 
    QMol_TDCI_SSO_2S.version;
end
end
methods (Access=protected)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the TDCI solver

    % Call parent documentation
    ref                 =   showDoc@QMol_TDCI_sympSplitOp(obj);
    
    % Implementation of the TDDFT model
    fprintf('  * Strang-splitting (a.k.a. Verlet) propagation scheme          2nd order\n');
    obj.version;

end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    % Parent-class components
    [~,PropNames]       =   QMol_TDCI_sympSplitOp.propertyNames;

    ClassName           =   'QMol_TDCI_SSO_2S';
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
    if obj.isHDH,       obj.applyExpH(QMol_cte.symp_2S2_sC1);
                        obj.applyExpD(QMol_cte.symp_2S2_sD1);
                        obj.applyExpH(QMol_cte.symp_2S2_sC2);
    else,               obj.applyExpD(QMol_cte.symp_2S2_sC1);
                        obj.applyExpH(QMol_cte.symp_2S2_sD1);
                        obj.applyExpD(QMol_cte.symp_2S2_sC2);
    end

    % Mask damping
    if ~isempty(obj.D),     obj.CI.wfcn     =   obj.D*obj.CI.wfcn;          end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

