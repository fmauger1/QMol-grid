classdef QMol_TDSE_SSO_4FR < QMol_TDSE_SSO
%QMol_TDSE_SSO_4FR symplectic split-operator TDSE propagator using a
%   Forest-Ruth scheme (4th order)

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
    fprintf('  * QMol_TDSE_SSO_4FR:\n');
    fprintf('      > Symplectic split-operator TDSE propagator\n'); 
    fprintf('      > Forest-Ruth scheme (4th order)\n'); 
    QMol_TDSE_SSO_4FR.version;
end
end
methods (Access=protected)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the TDSE solver

    % Call parent documentation
    ref                 =   showDoc@QMol_TDSE_sympSplitOp(obj);
    
    % Implementation of the TDDFT model
    fprintf('  * Forest-Ruth propagation scheme [Forest 1990]                 4th order\n');
    obj.version;

    ref                 =   [ref,{'Forest 1990'}];

end
end
%% Properties (scheme coefficients) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Constant,Access=private)
    sC1                 =   .5               / (2-2^(1/3))
    sC2                 =   .5 * (1-2^(1/3)) / (2-2^(1/3))
    sC3                 =   .5 * (1-2^(1/3)) / (2-2^(1/3))
    sC4                 =   .5               / (2-2^(1/3))
    sD1                 =   1                / (2-2^(1/3))
    sD2                 =          -2^(1/3)  / (2-2^(1/3))
    sD3                 =   1                / (2-2^(1/3))
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    % Parent-class components
    [~,PropNames]       =   QMol_TDSE_sympSplitOp.propertyNames;

    ClassName           =   'QMol_TDSE_SSO_4FR';
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
    
    if obj.isVTV,   obj.setExpV(obj.sC1);   if obj.FG ~= 2, obj.setExpT(obj.sD1);   end
    else,           obj.setExpT(obj.sC1);                                           end

end
function applyTimeStep(obj,t)
%applyTimeStep applies one iteration of the Strang splitting operator
    
    % Update t
    obj.t_              =   t;

    % Time propagation
    if obj.isVTV,                                       obj.applyExpV;      % expV is still good
        if obj.FG == 2, obj.setExpT(obj.sD1);   end;    obj.applyExpT;      % expT might still be good
                        obj.setExpV(obj.sC2);           obj.applyExpV;
                        obj.setExpT(obj.sD2);           obj.applyExpT;
                        obj.setExpV(obj.sC3);           obj.applyExpV;
                        obj.setExpT(obj.sD3);           obj.applyExpT;
                        obj.setExpV(obj.sC4);           obj.applyExpV;
    else,                                               obj.applyExpT;      % expT is still good
                        obj.setExpV(obj.sD1);           obj.applyExpV;
                        obj.setExpT(obj.sC2);           obj.applyExpT;
                        obj.setExpV(obj.sD2);           obj.applyExpV;
                        obj.setExpT(obj.sC3);           obj.applyExpT;
                        obj.setExpV(obj.sD3);           obj.applyExpV;
                        obj.setExpT(obj.sC4);           obj.applyExpT;
    end

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

