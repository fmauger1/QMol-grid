classdef QMol_TDDFT_ESSO_2S < QMol_TDDFT_ESSO
%QMol_TDDFT_ESSO_2S extended symplectic split-operator TDDFT propagator 
%   using a Strang-splitting scheme (2nd order)

%   Version     Date        Author
%   01.23.000   07/22/2025  F. Mauger
%       Creation (from QMol_TDDFT_SSO_2S version 01.21.000)
%               07/23/2050  F. Mauger
%       Add HRH, TVR, and TRV split motifs
%   01.23.001   07/25/2024  F. Mauger
%       Fix inclusion of external field and CAP without potential split

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.001','07/30/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDDFT})
function showInfo
    fprintf('  * QMol_TDDFT_ESSO_2S:\n');
    fprintf('      > Extended symplectic split-operator TDDFT propagator\n'); 
    fprintf('      > Strang-splitting scheme (2nd order)\n'); 
    QMol_TDDFT_ESSO_2S.version;
end
end
methods (Access=protected)
function ref = showDoc(obj)
%showDoc displays the documentation reflecting the specific implementation
%   of the TDDFT solver

    % Call parent documentation
    ref                 =   showDoc@QMol_TDDFT_extSymp(obj);
    
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
    [~,PropNames]       =   QMol_TDDFT_extSymp.propertyNames;

    ClassName           =   'QMol_TDDFT_ESSO_2S';
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
    obj.setExpT(QMol_cte.symp_2Sn_a1);
end
function applyTimeStep(obj,t)
%applyTimeStep applies one iteration of the Strang splitting operator
    
    % Update t
    obj.t_              =   t;

    % Time propagation
    switch obj.algo
        case 1 % chi_f = T1 V1 T2 V2 R      chi_b = R V2 T2 V1 T1 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a1;
            obj.applyExpT1(c);          obj.applyExpV1(c,true);
            obj.applyExpT2(c);          obj.applyExpV2(c,true);
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpR(QMol_cte.symp_2Sn_a1+QMol_cte.symp_2Sn_a2);
                if obj.FG == 2,     obj.setExpT(QMol_cte.symp_2Sn_a2);      end
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a2;
            obj.applyExpV2(c,false);    obj.applyExpT2(c);
            obj.applyExpV1(c,false);    obj.applyExpT1(c);
        case 2 % chi_f = T1 V1 R T2 V2      chi_b = V2 T2 R V1 T1 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a1;
            obj.applyExpT1(c);          obj.applyExpV1(c,true);  
            obj.applyExpR(c);   if obj.FG == 2,     obj.setExpT(c);         end
            obj.applyExpT2(c);
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpV22(QMol_cte.symp_2Sn_a1,QMol_cte.symp_2Sn_a2)
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a2;
            obj.applyExpT2(c);
            obj.applyExpR(c);   if obj.FG == 2,     obj.setExpT(c);         end
            obj.applyExpV1(c,false);    obj.applyExpT1(c);
        case 3 % chi_f = T1 T2 V1 V2 R      chi_b = R V2 V1 T2 T1 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a1;
            obj.applyExpT1(c);          obj.applyExpT2(c);
            obj.applyExpV1(c,true);     obj.applyExpV2(c,true);
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpR(QMol_cte.symp_2Sn_a1+QMol_cte.symp_2Sn_a2)
                if obj.FG == 2,     obj.setExpT(QMol_cte.symp_2Sn_a2);      end
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a2;
            obj.applyExpV2(c,false);    obj.applyExpV1(c,false);
            obj.applyExpT2(c);          obj.applyExpT1(c);
        case 4 % chi_f = T1 T2 R V1 V2      chi_b = V2 V1 R T2 T1 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a1;
            obj.applyExpT1(c);          obj.applyExpT2(c);
            obj.applyExpR(c);   if obj.FG == 2,     obj.setExpT(c);         end
            obj.applyExpV1(c,true);
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpV22(QMol_cte.symp_2Sn_a1,QMol_cte.symp_2Sn_a2)
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a2;
            obj.applyExpV1(c,false);
            obj.applyExpR(c);   if obj.FG == 2,     obj.setExpT(c);         end
            obj.applyExpT2(c);          obj.applyExpT1(c);
        case 5 % chi_f = T1 V2 T2 V1 R      chi_b = R V1 T2 V2 T1 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a1;
            obj.applyExpT1(c);          obj.applyExpV2(c,true);
            obj.applyExpT2(c);          obj.applyExpV1(c,true);
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpR(QMol_cte.symp_2Sn_a1+QMol_cte.symp_2Sn_a2);
                if obj.FG == 2,     obj.setExpT(QMol_cte.symp_2Sn_a2);      end
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a2;
            obj.applyExpV1(c,false);    obj.applyExpT2(c);
            obj.applyExpV2(c,false);    obj.applyExpT1(c);
        case 6 % chi_f = T2 V1 R T1 V2      chi_b = V2 T1 R V1 T2 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a1;
            obj.applyExpT2(c);          obj.applyExpV1(c,true);  
            obj.applyExpR(c);   if obj.FG == 2,     obj.setExpT(c);         end
            obj.applyExpT1(c);
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpV22(QMol_cte.symp_2Sn_a1,QMol_cte.symp_2Sn_a2)
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c   =   QMol_cte.symp_2Sn_a2;
            obj.applyExpT1(c);
            obj.applyExpR(c);   if obj.FG == 2,     obj.setExpT(c);         end
            obj.applyExpV1(c,false);    obj.applyExpT2(c);
        otherwise  % ======================================================
            error('QMol:QMol_TDDFT_ESSO_2S:applyTimeStep','Unexpected error; contact a developer');
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

