classdef QMol_TDDFT_ESSO < QMol_TDDFT_extSymp
%QMol_TDDFT_ESSO dimension-specific methods for extended symplectic split-
%   operator TDDFT propagators in one dimension

%   Version     Date        Author
%   01.23.000   07/22/2025  F. Mauger
%       Creation (from QMol_TDDFT_SSO version 01.21.000)
%   01.23.001   07/25/2024  F. Mauger
%       Fix inclusion of external field and CAP without potential split
%   01.23.002   08/06/2025  F. Mauger
%       Skip mixing (restrain) for omega = 0

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.002','08/06/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_TDDFT})
function showInfo
    fprintf('  * QMol_TDDFT_ESSO:\n');
    fprintf('      > TDDFT propagator\n'); 
    fprintf('      > Extended symplectic split-operator\n'); 
    fprintf('      > 1D-specific methods\n'); 
    QMol_TDDFT_ESSO.version;
end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N/A

%% Propagate TDDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=protected)
function saveOutputIonization(obj,K,t) %===================================
%saveOutputIonization

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.setOutputOrbital;
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

    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %saveOutputResults already limits call of the various methods to when
    %saving is required
    obj.setOutputOrbital;

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

    % Initialization
    obj.setOutputOrbital;

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
    DE                  =   .5*obj.xi;                                      % double counting the energy in the extended phase space

end
% Factorized chi_f+chi_b propagator 
function chiFchiB(obj,c,opt)
%chiFchiB performs one step of the chi_f + chi_p propagation for the
%   specified pair of parameters. The parameter opt specifies
%     -1  <=>  initial step of the propagation
%      0  <=>  middle step in the propagation
%      1  <=>  final step in the propagation

    switch obj.algo
        case 1 % chi_f = T1 V1 T2 V2 R      chi_b = R V2 T2 V1 T1 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if opt > -1,                    obj.setExpT(c(1));              else
            obj.applyExpT1(c(1));   end,    obj.applyExpV1(c(1),true);
            obj.applyExpT2(c(1));           obj.applyExpV2(c(1),true);
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpR(c(1)+c(2));
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                            obj.setExpT(c(2));
            obj.applyExpV2(c(2),false);     obj.applyExpT2(c(2));
            obj.applyExpV1(c(2),false);
            if opt == 1,                    obj.applyExpT1(c(2));           else %#ok<ALIGN>
            % chi_b + chi_f refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.setExpT(c(2)+c(3));         obj.applyExpT1(c(2)+c(3));      end

        case 2 % chi_f = T1 V1 R T2 V2      chi_b = V2 T2 R V1 T1 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if opt > -1,                    obj.setExpT(c(1));              else
            obj.applyExpT1(c(1));   end,    obj.applyExpV1(c(1),true); 
            obj.applyExpR(c(1));    if obj.FG == 2,     obj.setExpT(c(1));  end
            obj.applyExpT2(c(1));   
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpV22(c(1),c(2))
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                            obj.setExpT(c(2));
            obj.applyExpT2(c(2));
            obj.applyExpR(c(2));    if obj.FG == 2,     obj.setExpT(c(2));  end
            obj.applyExpV1(c(2),false);
            if opt == 1,                    obj.applyExpT1(c(2));           else %#ok<ALIGN>
            % chi_b + chi_f refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.setExpT(c(2)+c(3));         obj.applyExpT1(c(2)+c(3));      end

        case 3 % chi_f = T1 T2 V1 V2 R      chi_b = R V2 V1 T2 T1 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if opt > -1,                    obj.setExpT(c(1));              else
            obj.applyExpT1(c(1));   end,    obj.applyExpT2(c(1));
            obj.applyExpV1(c(1),true);      obj.applyExpV2(c(1),true);
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpR(c(1)+c(2));
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpV2(c(2),false);     obj.applyExpV1(c(2),false);     obj.setExpT(c(2));
            obj.applyExpT2(c(2));
            if opt == 1,                    obj.applyExpT1(c(2));           else %#ok<ALIGN>
            % chi_b + chi_f refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.setExpT(c(2)+c(3));         obj.applyExpT1(c(2)+c(3));      end

        case 4 % chi_f = T1 T2 R V1 V2      chi_b = V2 V1 R T2 T1 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if opt > -1,                    obj.setExpT(c(1));              else
            obj.applyExpT1(c(1));   end,    obj.applyExpT2(c(1));
            obj.applyExpR(c(1));            obj.applyExpV1(c(1),true);
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpV22(c(1),c(2))
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpV1(c(2),false);     obj.applyExpR(c(2));
            obj.setExpT(c(2));              obj.applyExpT2(c(2));
            if opt == 1,                    obj.applyExpT1(c(2));           else %#ok<ALIGN>
            % chi_b + chi_f refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.setExpT(c(2)+c(3));         obj.applyExpT1(c(2)+c(3));      end

        case 5 % chi_f = T1 V2 T2 V1 R      chi_b = R V1 T2 V2 T1 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if opt > -1,                    obj.setExpT(c(1));              else
            obj.applyExpT1(c(1));   end,    obj.applyExpV2(c(1),true);
            obj.applyExpT2(c(1));           obj.applyExpV1(c(1),true);
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpR(c(1)+c(2));
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                            obj.setExpT(c(2));
            obj.applyExpV1(c(2),false);     obj.applyExpT2(c(2));
            obj.applyExpV2(c(2),false);
            if opt == 1,                    obj.applyExpT1(c(2));           else %#ok<ALIGN>
            % chi_b + chi_f refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.setExpT(c(2)+c(3));         obj.applyExpT1(c(2)+c(3));      end

        case 6 % chi_f = T2 V1 R T1 V2      chi_b = V2 T1 R V1 T2 =========
            % chi_f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if opt > -1,                    obj.setExpT(c(1));              else
            obj.applyExpT2(c(1));   end,    obj.applyExpV1(c(1),true); 
            obj.applyExpR(c(1));    if obj.FG == 2,     obj.setExpT(c(1));  end
            obj.applyExpT1(c(1));   
            % chi_f + chi_b refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.applyExpV22(c(1),c(2))
            % chi_b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                            obj.setExpT(c(2));
            obj.applyExpT1(c(2));
            obj.applyExpR(c(2));    if obj.FG == 2,     obj.setExpT(c(2));  end
            obj.applyExpV1(c(2),false);
            if opt == 1,                    obj.applyExpT2(c(2));           else %#ok<ALIGN>
            % chi_b + chi_f refactorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            obj.setExpT(c(2)+c(3));         obj.applyExpT2(c(2)+c(3));      end

        otherwise  % ======================================================
            error('QMol:QMol_TDDFT_ESSO_4BM:applyTimeStep','Unexpected error; contact a developer');
    end
end
% Elemental propagation bits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function applyExpT1(obj,c),      obj.applyExpT(obj.p1,c);        end %=====
function applyExpT2(obj,c),      obj.applyExpT(obj.p2,c);        end %=====
function applyExpV1(obj,c,isFwd) %=========================================
    if obj.spltV,                                                                           if isFwd
        obj.applyExpVExplicit(obj.p1,       c);     obj.applyExpVImplicit(obj.p1,obj.p2,c); else
        obj.applyExpVImplicit(obj.p1,obj.p2,c);     obj.applyExpVExplicit(obj.p1,       c); end
    else
        obj.applyExpVFull(    obj.p1,obj.p2,c);
    end
end
function applyExpV2(obj,c,isFwd) %=========================================
    if obj.spltV,                                                                           if isFwd
        obj.applyExpVExplicit(obj.p2,       c);     obj.applyExpVImplicit(obj.p2,obj.p1,c); else
        obj.applyExpVImplicit(obj.p2,obj.p1,c);     obj.applyExpVExplicit(obj.p2,       c); end
    else
        obj.applyExpVFull(    obj.p2,obj.p1,c);
    end
end
function applyExpV22(obj,c1,c2) %==========================================
    if obj.spltV
        obj.applyExpVExplicit(obj.p2,       c1)
        obj.applyExpVImplicit(obj.p2,obj.p1,c1+c2)
        obj.applyExpVExplicit(obj.p2,       c2)
    else
        obj.applyExpVFull(    obj.p2,obj.p1,c1+c2);
    end
end
function applyExpT(obj,p,c) %==============================================
%applyExpT apply the Liouville operator for the kinetic term on input
%   orbitals

    % Energy conservation
    if obj.FG == 2   &&   obj.sEDFT
        % Update xi (input-energy variable)
        E               =   0;
        if obj.DFT.isSpinPol
            % Spin polarized
            for k = 1:obj.nKSO(1),  E   =   E + obj.DFT.occ{1}(k)*sum(conj(p.KSOup(:,k)).*ifft(obj.DFT.disc.D.*fft(p.KSOup(:,k))),'all'); end
            for k = 1:obj.nKSO(2),  E   =   E + obj.DFT.occ{2}(k)*sum(conj(p.KSOdw(:,k)).*ifft(obj.DFT.disc.D.*fft(p.KSOdw(:,k))),'all'); end
        else
            % Spin restricted
            for k = 1:obj.nKSO,     E   =   E + obj.DFT.occ(k)   *sum(conj(p.KSO(  :,k)).*ifft(obj.DFT.disc.D.*fft(p.KSO(  :,k))),'all'); end
        end

        E               =   obj.FE*(real(-1i*E*obj.DFT.disc.dx) + obj.FA*sum(obj.DFT.Ntot));

        obj.xi          =   obj.xi + c*obj.dt_*E;
    end

    % Propagate orbitals
    if obj.DFT.isSpinPol
        % Spin polarized
        for k = 1:obj.nKSO(1),  p.KSOup(:,k)  =   ifft(obj.expT.*fft(p.KSOup(:,k))); end
        for k = 1:obj.nKSO(2),  p.KSOdw(:,k)  =   ifft(obj.expT.*fft(p.KSOdw(:,k))); end
    else
        % Spin restricted
        for k = 1:obj.nKSO,     p.KSO(  :,k)  =   ifft(obj.expT.*fft(p.KSO(  :,k))); end
    end
end
function applyExpVFull(obj,p,p2,c) %=======================================
%applyExpVfull apply the Liouville operator for the entire potential term
%   on the input orbital, kernel, and step

    % Initialization
    h                   =   c*obj.dt_;

    obj.DFT.KSO         =   p2;
    obj.DFT.Vks         =   obj.DFT.getPotential([],obj.DFT.Vks);
    if ~isempty(obj.ABC),   obj.ABC.getPotential(obj.DFT.Vks);              end
    if obj.FG == 1,                                                         if obj.DFT.isSpinPol
        % Electric field
        obj.DFT.Vks.Vup =   obj.DFT.Vks.Vup + obj.FE*obj.X;
        obj.DFT.Vks.Vdw =   obj.DFT.Vks.Vdw + obj.FE*obj.X;                 else
        obj.DFT.Vks.V   =   obj.DFT.Vks.V   + obj.FE*obj.X;                 end
    end
    obj.DFT.setPotentialKernel;

    % Apply potential operator
    if obj.DFT.isSpinPol,                                                                               for k = 1:obj.nKSO(1)
        % Apply potential operator
        p.KSOup(:,k)    =  p.KSOup(:,k) - 1i*h*obj.DFT.Vks.applyPotential(p2.KSOup(:,k),true);           end,    for k = 1:obj.nKSO(2)
        p.KSOdw(:,k)    =  p.KSOdw(:,k) - 1i*h*obj.DFT.Vks.applyPotential(p2.KSOdw(:,k),false);          end
    else,                                                                                               for k = 1:obj.nKSO
        % Apply potential operator
        p.KSO(  :,k)    =  p.KSO(  :,k) - 1i*h*obj.DFT.Vks.applyPotential(p2.KSO(  :,k));                end
    end

    % Energy conservation
    if obj.FG == 1   &&   obj.sEDFT
        % Electric field derivative
        obj.FDE         =   obj.getFDE(obj.t_);

        % Update xi (input-energy variable ... density calculated for potential)
        if obj.DFT.isSpinPol,   E   =   sum(obj.X.*obj.DFT.rho.rhoUp,'all')+sum(obj.X.*obj.DFT.rho.rhoDw,'all');
        else,                   E   =   sum(obj.X.*obj.DFT.rho.rho,  'all');    end
    
        obj.xi      =   obj.xi - h*obj.FDE*E * obj.DFT.disc.dx;
    end

end
function applyExpVExplicit(obj,p,c) %======================================
%applyExpVexplicit apply the Liouville operator for the explicit part of 
%   the potential term on the input orbital and step

    % Initialization
    obj.setExpV(p,c);

    % Propagate orbitals
    if obj.DFT.isSpinPol
        % Spin polarized
        for k = 1:obj.nKSO(1),  p.KSOup(:,k)  =   obj.expVup.*p.KSOup(:,k); end
        for k = 1:obj.nKSO(2),  p.KSOdw(:,k)  =   obj.expVdw.*p.KSOdw(:,k); end
    else
        % Spin restricted
        for k = 1:obj.nKSO,     p.KSO(  :,k)  =   obj.expV  .*p.KSO(  :,k); end
    end

    % Energy conservation
    if obj.FG == 1   &&   obj.sEDFT
        % Electric field derivative
        obj.FDE         =   obj.getFDE(obj.t_);

        % Update xi (input-energy variable ... density calculated for potential)
        if obj.DFT.isSpinPol,   E   =   sum(obj.X.*obj.DFT.rho.rhoUp,'all')+sum(obj.X.*obj.DFT.rho.rhoDw,'all');
        else,                   E   =   sum(obj.X.*obj.DFT.rho.rho,  'all');    end
    
        obj.xi      =   obj.xi - c*obj.dt_*obj.FDE*E * obj.DFT.disc.dx;
    end

end
function applyExpVImplicit(obj,p,p2,c) %=======================================
%applyExpVfull apply the Liouville operator for the entire potential term
%   on the input orbital, kernel, and step

    % Initialization
    h                   =   c*obj.dt_;

    obj.DFT.KSO         =   p2;
    obj.DFT.Vks         =   obj.DFT.getPotential([],obj.DFT.Vks);
    obj.DFT.setPotentialKernel;

    % Apply the implicit part of the potential
    for l = 1:numel(obj.DFT.Vks.Vimp),                                      if obj.DFT.isSpinPol
        % Spin polarized
        for k = 1:obj.nKSO(1),  p.KSOup(:,k)    =  p.KSOup(:,k) - 1i*h*obj.DFT.Vks.Vimp{l}(p2.KSOup(:,k),true);     end
        for k = 1:obj.nKSO(2),  p.KSOdw(:,k)    =  p.KSOdw(:,k) - 1i*h*obj.DFT.Vks.Vimp{l}(p2.KSOdw(:,k),false);    end
                                                                            else
        % Spin restricted
        for k = 1:obj.nKSO,     p.KSO(  :,k)    =  p.KSO(  :,k) - 1i*h*obj.DFT.Vks.Vimp{l}(p2.KSO(  :,k));          end
                                                                            end
    end
    % end
end
function applyExpR(obj,c)
%applyRestrain apply the restrain operation on the wave functions

    % Initialization
    h                   =   c*obj.dt_;
    
    % Mixing wave functions
    if obj.w ~= 0, if obj.DFT.isSpinPol                                     %#ok<ALIGN>
        obj.pDFT.KSOup  =   .5*(obj.p1.KSOup+obj.p2.KSOup);                 obj.pDFT.KSOdw  =   .5*(obj.p1.KSOdw+obj.p2.KSOdw);
        obj.p1.KSOup    =   obj.pDFT.KSOup;                                 obj.p1.KSOdw    =   obj.pDFT.KSOdw;
        obj.pDFT.KSOup  =   obj.pDFT.KSOup-obj.p2.KSOup;                    obj.pDFT.KSOdw  =   obj.pDFT.KSOdw-obj.p2.KSOdw;
        obj.p2.KSOup    =   obj.p1.KSOup;                                   obj.p2.KSOdw    =   obj.p1.KSOdw;

        obj.p1.KSOup    =   obj.p1.KSOup + exp(2i*obj.w*h)*obj.pDFT.KSOup;  obj.p1.KSOdw    =   obj.p1.KSOdw + exp(2i*obj.w*h)*obj.pDFT.KSOdw;
        obj.p2.KSOup    =   obj.p2.KSOup - exp(2i*obj.w*h)*obj.pDFT.KSOup;  obj.p2.KSOdw    =   obj.p2.KSOdw - exp(2i*obj.w*h)*obj.pDFT.KSOdw;
    else
        obj.pDFT.KSO    =   .5*(obj.p1.KSO+obj.p2.KSO);
        obj.p1.KSO      =   obj.pDFT.KSO;
        obj.pDFT.KSO    =   obj.pDFT.KSO-obj.p2.KSO;
        obj.p2.KSO      =   obj.p1.KSO;

        obj.p1.KSO      =   obj.p1.KSO + exp(2i*obj.w*h)*obj.pDFT.KSO;
        obj.p2.KSO      =   obj.p2.KSO - exp(2i*obj.w*h)*obj.pDFT.KSO;
    end, end

    % Update time and external field
    obj.t_              =   obj.t_ + h;
    
    if obj.FG == 1
        obj.FE          =   obj.getFE(obj.t_);
    elseif obj.FG == 2,                                                     if obj.uA
        % Update A
        obj.FA          =   obj.EF.potentialVector(obj.t_);                 else
        obj.FA          =   obj.FA - h/6*(...
                                obj.EF.electricField(obj.t_-h) + ...
                                4*obj.EF.electricField(obj.t_-.5*h) + ...
                                obj.EF.electricField(obj.t_));              end
                                                                            if (obj.sEDFT || obj.sAcc)
        % Update E (time variable)
        obj.FE          =   obj.getFE(obj.t_);                              end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

