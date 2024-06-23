classdef QMol_TDDFT_SSO_4BM < QMol_TDDFT_SSO
%QMol_TDDFT_SSO_4BM 4th order Blanes and Moan optimized symplectic split-
%   operator scheme

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDDFT})
function showInfo
    fprintf('  * QMol_TDDFT_SSO_4BM:\n');
    fprintf('      > Symplectic split-operator TDDFT propagator\n'); 
    fprintf('      > Optimized Blanes and Moan scheme (4th order)\n'); 
    QMol_TDDFT_SSO_4BM.version;
end
end
methods (Access=protected)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the TDDFT solver

    % Call parent documentation
    ref                 =   showDoc@QMol_TDDFT_sympSplitOp(obj);
    
    % Implementation of the TDDFT model
    fprintf('  * Optimized Blanes and Moan scheme [Blanes 2002]               4th order\n');
    obj.version;

    ref                 =   [ref,{'Blanes 2002'}];

end
end
%% Properties (scheme coefficients) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Constant,Access=private)
    sC1                 =   0.0792036964311957
    sD1                 =   0.209515106613362
    sC2                 =   0.353172906049774
    sD2                 =  -0.143851773179818
    sC3                 =  -0.0420650803577195
    sD3                 =   0.434336666566456
    sC4                 =   0.2193769557534996
    sD4                 =   0.434336666566456
    sC5                 =  -0.0420650803577195
    sD5                 =  -0.143851773179818
    sC6                 =   0.353172906049774
    sD6                 =   0.209515106613362
    sC7                 =   0.0792036964311957
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    % Parent-class components
    [~,PropNames]       =   QMol_TDDFT_sympSplitOp.propertyNames;

    ClassName           =   'QMol_TDDFT_SSO_4BM';
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
    
    if obj.isVTV,   obj.setExpV(obj.sC1);   if obj.FG ~= 2, obj.setExpT(obj.sD1);   end %#ok<ALIGN> 
    else,           obj.setExpT(obj.sC1);
                    obj.DFT.Vks =   obj.DFT.getPotential([],obj.DFT.Vks);           end % Initialize V_KS

    if ~isempty(obj.DFT.Vks.Vimp)
        % Implicit potential component
        error('QMol:TDDFT:SSOBlanesMoan64:implicitPotential', ...
            'The Blanes and Moan 6 (4) propagator does not support implicit Kohn-Sham potentials');
    end

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
                        obj.setExpT(obj.sD4);           obj.applyExpT;
                        obj.setExpV(obj.sC5);           obj.applyExpV;
                        obj.setExpT(obj.sD5);           obj.applyExpT;
                        obj.setExpV(obj.sC6);           obj.applyExpV;
                        obj.setExpT(obj.sD6);           obj.applyExpT;
                        obj.setExpV(obj.sC7);           obj.applyExpV;
    else,                                               obj.applyExpT;      % expT is still good
                        obj.setExpV(obj.sD1);           obj.applyExpV;
                        obj.setExpT(obj.sC2);           obj.applyExpT;
                        obj.setExpV(obj.sD2);           obj.applyExpV;
                        obj.setExpT(obj.sC3);           obj.applyExpT;
                        obj.setExpV(obj.sD3);           obj.applyExpV;
                        obj.setExpT(obj.sC4);           obj.applyExpT;
                        obj.setExpV(obj.sD4);           obj.applyExpV;
                        obj.setExpT(obj.sC5);           obj.applyExpT;
                        obj.setExpV(obj.sD5);           obj.applyExpV;
                        obj.setExpT(obj.sC6);           obj.applyExpT;
                        obj.setExpV(obj.sD6);           obj.applyExpV;
                        obj.setExpT(obj.sC7);           obj.applyExpT;
    end

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

