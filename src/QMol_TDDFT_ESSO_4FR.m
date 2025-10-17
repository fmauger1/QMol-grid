classdef QMol_TDDFT_ESSO_4FR < QMol_TDDFT_ESSO
%QMol_TDDFT_ESSO_4FR symplectic split-operator TDDFT propagator using a
%   Forest-Ruth scheme (4th order)

%   Version     Date        Author
%   01.23.000   07/23/2025  F. Mauger
%       Creation (from QMol_TDDFT_SSO_4FR version 01.21.000)

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.000','07/23/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDDFT})
function showInfo
    fprintf('  * QMol_TDDFT_ESSO_4FR:\n');
    fprintf('      > Extended symplectic split-operator TDDFT propagator\n'); 
    fprintf('      > Forest-Ruth scheme (4th order)\n'); 
    QMol_TDDFT_ESSO_4FR.version;
end
end
methods (Access=protected)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the TDDFT solver

    % Call parent documentation
    ref                 =   showDoc@QMol_TDDFT_extSymp(obj);
    
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
    [~,PropNames]       =   QMol_TDDFT_extSymp.propertyNames;

    ClassName           =   'QMol_TDDFT_ESSO_4FR';
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
    obj.setExpT(QMol_cte.symp_4FRn_a1);
end
function applyTimeStep(obj,t)
%applyTimeStep applies one iteration of propagator
    
    % Update t
    obj.t_              =   t;

    % Time propagation
    obj.chiFchiB([QMol_cte.symp_4FRn_a1 QMol_cte.symp_4FRn_a2 QMol_cte.symp_4FRn_a3],-1);
    obj.chiFchiB([QMol_cte.symp_4FRn_a3 QMol_cte.symp_4FRn_a4 QMol_cte.symp_4FRn_a5], 0);
    obj.chiFchiB([QMol_cte.symp_4FRn_a5 QMol_cte.symp_4FRn_a6                      ], 1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

