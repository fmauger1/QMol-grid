classdef QMol_TDSE_SSO < QMol_TDSE_sympSplitOp
%QMol_TDSE_SSO dimension-specific methods for symplectic split-operator
%   TDSE propagators in one dimension

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
    fprintf('  * QMol_TDSE_SSO:\n');
    fprintf('      > TDSE propagator\n'); 
    fprintf('      > Symplectic split-operator\n'); 
    fprintf('      > 1D-specific methods\n'); 
    QMol_TDSE_SSO.version;
end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N/A

%% Propagate TDDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=protected)
function saveOutputIonization(obj,K,t) %=====================================
%saveOutputIonization

    % Ionization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ~isempty(obj.oIon.indexWaveFunction)
        I               =   sum(abs(obj.SE.wfcn.wfcn(:,obj.oIon.indexWaveFunction)).^2,1)*obj.dv;
        obj.oIon.waveFunction(:,obj.oIon.n)      =   max(1-I,0);
    end

    % Add external field ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.sEF,         obj.addOutputExternalField('oIon',K,t);             end

    % Update counter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.oIon.n          =   obj.oIon.n + 1;

end
function saveOutputDipole(obj,K,t) %=======================================
%saveOutputDipole

    % Dipole signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if K == obj.oDip.ind(obj.oDip.n)
        % dipole
        for l = 1:numel(obj.oDip.indexWaveFunction)
            % Spin restricted
            obj.oDip.waveFunction_x(l,obj.oDip.n)   =   sum(obj.X .* abs(obj.SE.wfcn.wfcn(  :,obj.oDip.indexWaveFunction(l))).^2,'all')*obj.dv;
        end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oDip',K,t);        end
        
        % Update counter
        obj.oDip.n      =   obj.oDip.n + 1;

    end
    % Dipole velocity signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if K == obj.oVel.ind(obj.oVel.n)
        % dipole
        for l = 1:numel(obj.oDip.indexWaveFunction)
            % Spin restricted
            obj.oVel.waveFunction_x(l,obj.oVel.n)   =   real(-1i* sum(conj(obj.SE.wfcn.wfcn(:,obj.oVel.indexWaveFunction(l))).*ifft(obj.SE.disc.D.*fft(obj.SE.wfcn.wfcn(:,obj.oVel.indexWaveFunction(l)))),'all') )*obj.dv;
        end

        if obj.FG == 2
            obj.oVel.waveFunction_x(:,obj.oVel.n)   =   obj.oVel.waveFunction_x(:,obj.oVel.n)+ obj.FA;
        end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oVel',K,t);        end
        
        % Update counter
        obj.oVel.n      =   obj.oVel.n + 1;

    end

    % Dipole acceleration signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if K == obj.oAcc.ind(obj.oAcc.n)
        % External field
        if obj.FG > 0,  E   =   obj.FE;
        else,           E   =   0;                                          end

        % Dipole acceleration
        for l = 1:numel(obj.oAcc.indexWaveFunction)
            obj.oAcc.waveFunction_x(l,obj.oAcc.n)   =  -sum(obj.SE.V.DV .* abs(obj.SE.wfcn.wfcn(:,obj.oAcc.indexWaveFunction(l))).^2,'all')*obj.dv - E;
        end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oAcc',K,t);        end

        % Update counter
        obj.oAcc.n      =   obj.oAcc.n + 1;
    end
end
function [E,DE] = getExternalFieldEnergy(obj,~)
%getExternalFieldEnergy

    % External-field energy
    if obj.FG == 1  % Length gauge
        E               =   sum(obj.X.*abs(obj.SE.wfcn.wfcn).^2,  'all');
        E               =   E * obj.FE * obj.dv;
    else            % No field / velocity Gauge (energy included in kinetic term
        E               =   0;
    end
    
    % Autonomization
    DE                  =   obj.xi;

end
function applyExpT(obj) %==================================================
%applyExpT apply the Liouville operator for the kinetic term
    
    % Update other components (not affecting the orbital propagation next)
    if obj.FG == 1
        % Update t (time variable)
        obj.t_          =   obj.t_ + obj.c_(1)*obj.dt_;
        obj.FE          =   obj.getFE(obj.t_);

    elseif obj.FG == 2   &&   obj.sESE
        % Update xi (input-energy variable)
        E               =   0;
        for k = 1:obj.nWfcn
            E           =   E + sum(conj(obj.SE.wfcn.wfcn(:,k)).*ifft(obj.SE.disc.D.*fft(obj.SE.wfcn.wfcn(:,k))),'all'); 
        end

        E               =   obj.FE*(real(-1i*E*obj.SE.disc.dx) + obj.FA*sum(obj.SE.N));

        obj.xi          =   obj.xi + obj.c_(1)*obj.dt_*E;                   %end
    end
    
    % Propagate orbitals
    for k = 1:obj.nWfcn
        obj.SE.wfcn.wfcn(:,k)   =   ifft(obj.expT.*fft(obj.SE.wfcn.wfcn(:,k))); 
    end
end
function applyExpV(obj) %==================================================
%applyExpV apply the Liouville operator for the potential term
    
    % Propagate orbitals
    for k = 1:obj.nWfcn
        obj.SE.wfcn.wfcn(:,k)   =   obj.expV  .*obj.SE.wfcn.wfcn(:,k);
    end

    % Update other components
    if obj.FG == 1   &&   obj.sESE
        % Electric field derivative
        obj.FDE         =   obj.getFDE(obj.t_);

        % Update xi (input-energy variable)
        E               =   sum(obj.X.*abs(obj.SE.wfcn.wfcn).^2,'all');
    
        obj.xi      =   obj.xi - obj.c_(2)*obj.dt_*obj.FDE*E * obj.SE.disc.dx;

    elseif obj.FG == 2,                                                             if obj.uA
        % Update A
        obj.FA          =   obj.EF.potentialVector(obj.t_ + obj.c_(2)*obj.dt_);     else
        obj.FA          =   obj.FA - obj.c_(2)*obj.dt_/6*(...
                                obj.EF.electricField(obj.t_) + ...
                                4*obj.EF.electricField(obj.t_+.5*obj.c_(2)*obj.dt_) + ...
                                obj.EF.electricField(obj.t_+obj.c_(2)*obj.dt_));    end
    
        % Update t (time variable)
        obj.t_          =   obj.t_ + obj.c_(2)*obj.dt_;                     if (obj.sESE || obj.sAcc)
        obj.FE          =   obj.getFE(obj.t_);                              end
        
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

