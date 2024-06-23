classdef QMol_TDDFT_SSO < QMol_TDDFT_sympSplitOp
%QMol_TDDFT_SSO dimension-specific methods for symplectic split-operator
%   TDDFT propagators in one dimension

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
    fprintf('  * QMol_TDDFT_SSO:\n');
    fprintf('      > TDDFT propagator\n'); 
    fprintf('      > Symplectic split-operator\n'); 
    fprintf('      > 1D-specific methods\n'); 
    QMol_TDDFT_SSO.version;
end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N/A

%% Propagate TDDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=protected)
function saveOutputIonization(obj,K,t) %=====================================
%saveOutputIonization

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.DFT.rho         =   obj.DFT.getDensity(obj.DFT.rho);
                                                                            if obj.DFT.isSpinPol
    % Total ionization signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.oIon.totalUp(obj.oIon.n)    =   max(obj.DFT.Ntot(1) - sum(obj.DFT.rho.rhoUp,'all')*obj.dv,0);
    obj.oIon.totalDown(obj.oIon.n)  =   max(obj.DFT.Ntot(2) - sum(obj.DFT.rho.rhoDw,'all')*obj.dv,0);
    obj.oIon.total(obj.oIon.n)      =   obj.oIon.totalUp(obj.oIon.n) + obj.oIon.totalDown(obj.oIon.n);
                                                                            else
    obj.oIon.total(obj.oIon.n)      =   max(obj.DFT.Ntot - sum(obj.DFT.rho.rho,'all')*obj.dv,0);
                                                                            end,    if obj.DFT.isSpinPol
    % Orbital-resolved ionization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ~isempty(obj.oIon.indexOrbitalUp)
        I               =   sum(abs(obj.DFT.KSO.KSOup(:,obj.oIon.indexOrbitalUp  )).^2,1)*obj.dv;
        obj.oIon.orbitalUp(:,obj.oIon.n)    =   max(1-I,0);
    end
    if ~isempty(obj.oIon.indexOrbitalDown)
        I               =   sum(abs(obj.DFT.KSO.KSOdw(:,obj.oIon.indexOrbitalDown)).^2,1)*obj.dv;
        obj.oIon.orbitalDown(:,obj.oIon.n)  =   max(1-I,0);
    end,                                                                            else
    if ~isempty(obj.oIon.indexOrbital)
        I               =   sum(abs(obj.DFT.KSO.KSO(:,  obj.oIon.indexOrbital    )).^2,1)*obj.dv;
        obj.oIon.orbital(:,obj.oIon.n)      =   max(1-I,0);
    end
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
        % Total dipole
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);                if obj.DFT.isSpinPol
        
        obj.oDip.totalUp(1,obj.oDip.n)      =   sum(obj.X .* obj.DFT.rho.rhoUp,'all') * obj.dv;
        obj.oDip.totalDown(1,obj.oDip.n)    =   sum(obj.X .* obj.DFT.rho.rhoDw,'all') * obj.dv;
        obj.oDip.total(:,obj.oDip.n)        =   obj.oDip.totalUp(:,obj.oDip.n) + obj.oDip.totalDown(:,obj.oDip.n);
                                                                            else
        obj.oDip.total(1,obj.oDip.n)        =   sum(obj.X .* obj.DFT.rho.rho,  'all') * obj.dv;
                                                                            end 
        % Orbital-resolved dipole
        if obj.DFT.isSpinPol,                                               for l = 1:numel(obj.oDip.indexOrbitalUp)                        %#ok<ALIGN> 
            % Spin up
            obj.oDip.orbitalUp_x(l,obj.oDip.n)  =   sum(obj.X .* abs(obj.DFT.KSO.KSOup(:,obj.oDip.indexOrbitalUp(  l))).^2,'all')*obj.dv;
                                                                            end, for l = 1:numel(obj.oDip.indexOrbitalDown)
            % Spin down
            obj.oDip.orbitalDown_x(l,obj.oDip.n)=   sum(obj.X .* abs(obj.DFT.KSO.KSOdw(:,obj.oDip.indexOrbitalDown(l))).^2,'all')*obj.dv;
                                                                                 end, else, for l = 1:numel(obj.oDip.indexOrbital)          %#ok<ALIGN> 
            % Spin restricted
            obj.oDip.orbital_x(l,obj.oDip.n)    =   sum(obj.X .* abs(obj.DFT.KSO.KSO(  :,obj.oDip.indexOrbital(    l))).^2,'all')*obj.dv;
        end, end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oDip',K,t);        end
        
        % Update counter
        obj.oDip.n      =   obj.oDip.n + 1;

    end
    % Dipole velocity signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if K == obj.oVel.ind(obj.oVel.n)
        % Total dipole velocity
        if obj.DFT.isSpinPol
            % Spin up
            V           =   NaN(obj.nKSO(1),1);                                                                                 for l = 1:obj.nKSO(1)
            V(l)        =   sum(conj(obj.DFT.KSO.KSOup(:,l)).*ifft(obj.DFT.disc.D.*fft(obj.DFT.KSO.KSOup(:,l))),'all');         end
            V           =   real(-1i*V*obj.dv);                                         if obj.FG == 2,     V   =   V + obj.FA; end
            obj.oVel.totalUp(      1,obj.oVel.n)    =   sum(obj.DFT.occ{1}(:).* V);     if ~isempty(obj.oVel.indexOrbitalUp)
            obj.oVel.orbitalUp_x(  :,obj.oVel.n)    =   V(obj.oVel.indexOrbitalUp);     end

            % Spin down
            V           =   NaN(obj.nKSO(2),1);                                                                                 for l = 1:obj.nKSO(2)
            V(l)        =   sum(conj(obj.DFT.KSO.KSOdw(:,l)).*ifft(obj.DFT.disc.D.*fft(obj.DFT.KSO.KSOdw(:,l))),'all');         end
            V           =   real(-1i*V*obj.dv);                                         if obj.FG == 2,     V   =   V + obj.FA; end
            obj.oVel.totalDown(    1,obj.oVel.n)    =   sum(obj.DFT.occ{2}(:).* V);     if ~isempty(obj.oVel.indexOrbitalDown)
            obj.oVel.orbitalDown_x(:,obj.oVel.n)    =   V(obj.oVel.indexOrbitalDown);   end
            
            % Total signal
            obj.oVel.total(        :,obj.oVel.n)    =   obj.oVel.totalUp(1,obj.oVel.n) + obj.oVel.totalDown(1,obj.oVel.n);
        else
            % Spin restricted
            V           =   NaN(obj.nKSO,1);                                                                                    for l = 1:obj.nKSO
            V(l)        =   sum(conj(obj.DFT.KSO.KSO(  :,l)).*ifft(obj.DFT.disc.D.*fft(obj.DFT.KSO.KSO(  :,l))),'all');         end
            V           =   real(-1i*V*obj.dv);                                         if obj.FG == 2,     V   =   V + obj.FA; end
            obj.oVel.total(        1,obj.oVel.n)    =   sum(obj.DFT.occ(:).* V);        if ~isempty(obj.oVel.indexOrbital)
            obj.oVel.orbital_x(    :,obj.oVel.n)    =   V(obj.oVel.indexOrbital);       end
        end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oVel',K,t);        end

        % Update counter
        obj.oVel.n      =   obj.oVel.n + 1;
    end

    % Dipole acceleration signal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if K == obj.oAcc.ind(obj.oAcc.n)
        % Kohn-Sham potential derivative
        if obj.sDip,    obj.DFT.DVks        =   obj.DFT.getPotentialGradient(obj.DFT.rho,obj.DFT.DVks,1);   % One-body density already computed for the dipole
        else,           obj.DFT.DVks        =   obj.DFT.getPotentialGradient([]         ,obj.DFT.DVks,1);   end
        if obj.FG > 0,  E                   =   obj.FE;                     else,   E       =   0;          end
                                                                            if obj.DFT.isSpinPol
        % Total dipole acceleration
        obj.oAcc.totalUp(  1,obj.oAcc.n)    =  -sum(obj.DFT.DVks.DVup .* obj.DFT.rho.rhoUp,'all') * obj.dv - E*obj.DFT.Ntot(1);
        obj.oAcc.totalDown(1,obj.oAcc.n)    =  -sum(obj.DFT.DVks.DVdw .* obj.DFT.rho.rhoDw,'all') * obj.dv - E*obj.DFT.Ntot(2);
        obj.oAcc.total(    1,obj.oAcc.n)    =   obj.oAcc.totalUp(1,obj.oAcc.n) + obj.oAcc.totalDown(1,obj.oAcc.n);
                                                                            else
        obj.oAcc.total(    1,obj.oAcc.n)    =  -sum(obj.DFT.DVks.DV   .* obj.DFT.rho.rho,  'all') * obj.dv - E*obj.DFT.Ntot;
                                                                            end
        % Orbital-resolved dipole acceleration
        if obj.DFT.isSpinPol,                                               for l = 1:numel(obj.oAcc.indexOrbitalUp)                        %#ok<ALIGN> 
            % Spin up
            obj.oAcc.orbitalUp_x(  l,obj.oAcc.n)=  -sum(obj.DFT.DVks.DVup .* abs(obj.DFT.KSO.KSOup(:,obj.oAcc.indexOrbitalUp(  l))).^2,'all')*obj.dv - E;
                                                                            end, for l = 1:numel(obj.oAcc.indexOrbitalDown)
            % Spin down
            obj.oAcc.orbitalDown_x(l,obj.oAcc.n)=  -sum(obj.DFT.DVks.DVdw .* abs(obj.DFT.KSO.KSOdw(:,obj.oAcc.indexOrbitalDown(l))).^2,'all')*obj.dv - E;
                                                                                 end, else, for l = 1:numel(obj.oAcc.indexOrbital)          %#ok<ALIGN> 
            % Spin restricted
            obj.oAcc.orbital_x(    l,obj.oAcc.n)=  -sum(obj.DFT.DVks.DV   .* abs(obj.DFT.KSO.KSO(  :,obj.oAcc.indexOrbital(    l))).^2,'all')*obj.dv - E;
        end, end
        
        % Add external field
        if obj.sEF,         obj.addOutputExternalField('oAcc',K,t);        end

        % Update counter
        obj.oAcc.n      =   obj.oAcc.n + 1;
    end
end
function [E,DE] = getExternalFieldEnergy(obj,~)
%getExternalFieldEnergy

    % External-field energy
    if obj.FG == 1  % Length gauge (rho was initialized in previous DFT-energy-component calculations)
        if obj.DFT.isSpinPol,   E   =  [sum(obj.X.*obj.DFT.rho.rhoUp,'all'); sum(obj.X.*obj.DFT.rho.rhoDw,'all')];
        else,                   E   =   sum(obj.X.*obj.DFT.rho.rho,  'all');    end
        E               =   E * obj.FE * obj.dv;
    else            % No field / velocity Gauge (energy included in kinetic term
        if obj.DFT.isSpinPol,   E   =   [0;0];
        else,                   E   =   0;                                  end
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

    elseif obj.FG == 2   &&   obj.sEDFT
        % Update xi (input-energy variable)
        E               =   0;
        if obj.DFT.isSpinPol
            % Spin polarized
            for k = 1:obj.nKSO(1),  E   =   E + obj.DFT.occ{1}(k)*sum(conj(obj.DFT.KSO.KSOup(:,k)).*ifft(obj.DFT.disc.D.*fft(obj.DFT.KSO.KSOup(:,k))),'all'); end
            for k = 1:obj.nKSO(2),  E   =   E + obj.DFT.occ{2}(k)*sum(conj(obj.DFT.KSO.KSOdw(:,k)).*ifft(obj.DFT.disc.D.*fft(obj.DFT.KSO.KSOdw(:,k))),'all'); end
        else
            % Spin restricted
            for k = 1:obj.nKSO,     E   =   E + obj.DFT.occ(k)   *sum(conj(obj.DFT.KSO.KSO(  :,k)).*ifft(obj.DFT.disc.D.*fft(obj.DFT.KSO.KSO(  :,k))),'all'); end
        end

        E               =   obj.FE*(real(-1i*E*obj.DFT.disc.dx) + obj.FA*sum(obj.DFT.Ntot));

        obj.xi          =   obj.xi + obj.c_(1)*obj.dt_*E;                   %end
    end
    
    % Propagate orbitals
    if obj.DFT.isSpinPol
        % Spin polarized
        for k = 1:obj.nKSO(1),  obj.DFT.KSO.KSOup(:,k)  =   ifft(obj.expT.*fft(obj.DFT.KSO.KSOup(:,k))); end
        for k = 1:obj.nKSO(2),  obj.DFT.KSO.KSOdw(:,k)  =   ifft(obj.expT.*fft(obj.DFT.KSO.KSOdw(:,k))); end
    else
        % Spin restricted
        for k = 1:obj.nKSO,     obj.DFT.KSO.KSO(  :,k)  =   ifft(obj.expT.*fft(obj.DFT.KSO.KSO(  :,k))); end
    end
end
function applyExpV(obj) %==================================================
%applyExpV apply the Liouville operator for the potential term
    
    % Propagate orbitals
    if obj.DFT.isSpinPol
        % Spin polarized
        for k = 1:obj.nKSO(1),  obj.DFT.KSO.KSOup(:,k)  =   obj.expVup.*obj.DFT.KSO.KSOup(:,k); end
        for k = 1:obj.nKSO(2),  obj.DFT.KSO.KSOdw(:,k)  =   obj.expVdw.*obj.DFT.KSO.KSOdw(:,k); end
    else
        % Spin restricted
        for k = 1:obj.nKSO,     obj.DFT.KSO.KSO(  :,k)  =   obj.expV  .*obj.DFT.KSO.KSO(  :,k); end
    end

    % Update other components
    if obj.FG == 1   &&   obj.sEDFT
        % Electric field derivative
        obj.FDE         =   obj.getFDE(obj.t_);

        % Update xi (input-energy variable)
        if obj.DFT.isSpinPol,   E   =   sum(obj.X.*obj.DFT.rho.rhoUp,'all')+sum(obj.X.*obj.DFT.rho.rhoDw,'all');
        else,                   E   =   sum(obj.X.*obj.DFT.rho.rho,  'all');    end
    
        obj.xi      =   obj.xi - obj.c_(2)*obj.dt_*obj.FDE*E * obj.DFT.disc.dx;

    elseif obj.FG == 2,                                                             if obj.uA
        % Update A
        obj.FA          =   obj.EF.potentialVector(obj.t_ + obj.c_(2)*obj.dt_);     else
        obj.FA          =   obj.FA - obj.c_(2)*obj.dt_/6*(...
                                obj.EF.electricField(obj.t_) + ...
                                4*obj.EF.electricField(obj.t_+.5*obj.c_(2)*obj.dt_) + ...
                                obj.EF.electricField(obj.t_+obj.c_(2)*obj.dt_));    end
    
        % Update t (time variable)
        obj.t_          =   obj.t_ + obj.c_(2)*obj.dt_;                     if (obj.sEDFT || obj.sAcc)
        obj.FE          =   obj.getFE(obj.t_);                              end
        
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

