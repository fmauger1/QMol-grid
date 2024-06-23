classdef QMol_TDSE_SSO_2S < QMol_TDSE_SSO
%QMol_TDSE_SSO_2S symplectic split-operator TDDFT propagator using a
%   Strang-splitting scheme (2nd order)

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDSE})
function showInfo
    fprintf('  * QMol_TDSE_SSO_2S:\n');
    fprintf('      > Symplectic split-operator TDSE propagator\n'); 
    fprintf('      > Strang-splitting scheme (2nd order)\n'); 
    QMol_TDSE_SSO_2S.version;
end
end
methods (Access=protected)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the TDSE solver

    % Call parent documentation
    ref                 =   showDoc@QMol_TDSE_sympSplitOp(obj);
    
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
    [~,PropNames]       =   QMol_TDSE_sympSplitOp.propertyNames;

    ClassName           =   'QMol_TDSE_SSO_2S';
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
    
    if obj.isVTV,   obj.setExpV(.5);    if obj.FG ~= 2, obj.setExpT(1 );    end
    else,           obj.setExpT(.5);                                        end

end
function applyTimeStep(obj,t)
%applyTimeStep applies one iteration of the Strang splitting operator
    
    % Update t
    obj.t_              =   t;

    % Time propagation
    if obj.isVTV,                                   obj.applyExpV;          % expV is still good
        if obj.FG == 2, obj.setExpT(1 );    end;    obj.applyExpT;
                        obj.setExpV(.5);            obj.applyExpV;
    else,                                           obj.applyExpT;          % expT is still good
                        obj.setExpV(1 );            obj.applyExpV;
        if obj.FG == 2, obj.setExpT(.5);    end;    obj.applyExpT;
    end

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

