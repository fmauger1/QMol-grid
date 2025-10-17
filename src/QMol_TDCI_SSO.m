classdef QMol_TDCI_SSO < QMol_TDCI_sympSplitOp
%QMol_TDCI_SSO dimension-specific methods for symplectic split-operator
%   TDCI propagators in one dimension

%   Version     Date        Author
%   01.23.000   06/04/2025  F. Mauger
%       Creation (from QMol_TDCI_SSO)
%   01.23.001   06/10/2025  F. Mauger
%       Fix dipole dignal for field-free propagation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.001','06/10/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDCI})
function showInfo
    fprintf('  * QMol_TDCI_SSO:\n');
    fprintf('      > TDCI propagator\n'); 
    fprintf('      > Symplectic split-operator\n'); 
    fprintf('      > 1D-specific methods\n'); 
    QMol_TDCI_SSO.version;
end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N/A

%% Propagate TDCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=protected)
function saveOutputIonization(obj,K,t) %=====================================
%saveOutputIonization

    % Ionization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ~isempty(obj.oIon.indexWaveFunction)
        I               =   sum(abs(obj.CI.wfcn(:,obj.oIon.indexWaveFunction)).^2,1);
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
        obj.oDip.waveFunction_x(:,obj.oDip.n)   =   real(sum(conj(obj.CI.wfcn(:,obj.oDip.indexWaveFunction)).*(obj.CI.DX*obj.CI.wfcn(:,obj.oDip.indexWaveFunction)),1));
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oDip',K,t);        end
        
        % Update counter
        obj.oDip.n      =   obj.oDip.n + 1;

    end
    % Dipole velocity signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if K == obj.oVel.ind(obj.oVel.n)
        if obj.FG == 1, H   =   obj.CI.CI + obj.FE*obj.CI.DX;
        else,           H   =   obj.CI.CI;                                  end
        H               =   H*obj.CI.DX - obj.CI.DX*H;

        obj.oVel.waveFunction_x(:,obj.oVel.n)   =  -imag(sum(conj(obj.CI.wfcn(:,obj.oDip.indexWaveFunction)).*(H*obj.CI.wfcn(:,obj.oDip.indexWaveFunction)),1));
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oVel',K,t);        end
        
        % Update counter
        obj.oVel.n      =   obj.oVel.n + 1;

    end

    % Dipole acceleration signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if K == obj.oAcc.ind(obj.oAcc.n)
        if obj.FG == 1, H   =   obj.CI.CI + obj.FE*obj.CI.DX;
        else,           H   =   obj.CI.CI;                                  end
        HD              =   H*obj.CI.DX - obj.CI.DX*H;
        HD              =   H*HD - HD*H;

        obj.oAcc.waveFunction_x(:,obj.oAcc.n)   =   -real(sum(conj(obj.CI.wfcn(:,obj.oDip.indexWaveFunction)).*(HD*obj.CI.wfcn(:,obj.oDip.indexWaveFunction)),1));
        
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
        E               =   real(sum(conj(obj.CI.wfcn).*(obj.CI.DX*obj.CI.wfcn),'all'));
        E               =   E * obj.FE;
    else            % No field
        E               =   0;
    end
    
    % Autonomization
    DE                  =   obj.xi;

end
function applyExpH(obj,c) %================================================
%applyExpH apply the Liouville operator for the field-free Hamiltonian
    
    % Update other components (not affecting the orbital propagation next)
    if obj.FG == 1
        % Update t (time variable)
        obj.t_          =   obj.t_ + c*obj.dt_;
        obj.FE          =   obj.getFE(obj.t_);
    end
    
    % Propagate wave functions
    for k = 1:obj.nWfcn
        obj.CI.wfcn(:,k)   =   expmv(obj.H0,obj.CI.wfcn(:,k),c*obj.dt_); 
    end
end
function applyExpD(obj,c) %================================================
%applyExpD apply the Liouville operator for the dipole coupling

    % Apply dipole-coupling term only if there is a field
    if obj.FG == 1
        % CI energy
        if obj.sECI
            % Electric field derivative
            obj.FDE         =   obj.getFDE(obj.t_);
    
            % Update xi (input-energy variable)
            E               =   real(sum(conj(obj.CI.wfcn).*(obj.CI.DX*obj.CI.wfcn),'all'));
        
            obj.xi          =   obj.xi - c*obj.dt_*obj.FDE*E;
        end

        % Propagate wave functions
        D                   =   -1i*obj.FE*obj.CI.DX;
        for k = 1:obj.nWfcn
            obj.CI.wfcn(:,k)=   expmv(D,obj.CI.wfcn(:,k),c*obj.dt_);
        end
    end
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

