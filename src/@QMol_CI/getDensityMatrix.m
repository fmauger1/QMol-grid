function [rUp,rDw] = getDensityMatrix(obj,psi)
%getDensityMatrix computes the one-particle reduced density matrix (1-RDM)
%   Use getDensityMatrix to calculate the 1-RMD associated with the wave
%   function(s) of a CI model.
%
%   [rhoUp,rhoDown] = obj.getDensityMatrix returns the 1-RDM matrices for
%   the up (rhoUp) and down (rhoDown) spin channels of the waveFunction
%   property. For a single wave function (column vector waveFunction), the 
%   1-RDM are size(configurationBasis,1)-by-size(configurationBasis,1)
%   matrices. If waveFunction contains several wave functions, the outputs
%   are arrays with rhoUp(:,:,k) and rhoDown(:,:,k) containing the 1-RDM 
%   matrices for the kth wave function.
%
%   [rhoUp,rhoDown] = obj.getDensityMatrix(psi) returns the 1-RDM matrices
%   for the input wave function psi instead of waveFunction. Each column in
%   psi specifies the population coefficients in the configuration-state 
%   basis configurationBasis.

    % Initialization
    if nargin == 1,     psi     =   obj.waveFunction;                   end
    N                   =   size(obj.CSB,1);
    nMO                 =   size(obj.SOB,2);
    nPsi                =   size(psi,2);

    rUp                 =   zeros(nMO,nMO,nPsi);
    rDw                 =   zeros(nMO,nMO,nPsi);

    % Calculate 1-RDM
    for k = 1:N
        % Diagonal elements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if nPsi == 1,                                                       if abs(psi(k).^2) > obj.tol
            for m = obj.CSB(k,:)
                if m > 0,   rUp( m, m)      =   rUp( m, m) + abs(psi(k)).^2;   
                else,       rDw(-m,-m)      =   rDw(-m,-m) + abs(psi(k)).^2;    end
            end,                                                            end
        else,                                                               if any(abs(psi(k,:).^2) > obj.tol)
            for m = obj.CSB(k,:)
                if m > 0,   rUp( m, m,:)    =   squeeze(rUp( m, m,:)) + abs(psi(k,:).').^2;   
                else,       rDw(-m,-m,:)    =   squeeze(rDw(-m,-m,:)) + abs(psi(k,:).').^2;     end
            end,                                                            end
        end

        % Off-diagonal elements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if nPsi == 1,                                                       if abs(psi(k)) > obj.tol
            for l = k+1:N,                                                  if abs(psi(k)*psi(l)) > obj.tol
                [a,r]       =   obj.getExcitationIndexes(obj.CSB(k,:),obj.CSB(l,:));
                if isscalar(a),                                                 if a > 0 && r > 0
                    rUp( a, r)      =   rUp( a, r) + (conj(psi(k)).*psi(l));
                    rUp( r, a)      =   rUp( r, a) + (conj(psi(l)).*psi(k));    elseif a < 0 && r < 0
                    rDw(-a,-r)      =   rDw(-a,-r) + (conj(psi(k)).*psi(l));
                    rDw(-r,-a)      =   rDw(-r,-a) + (conj(psi(l)).*psi(k));    end
                end,                                                        end
            end,                                                            end
        else,                                                               if any(abs(psi(k,:)) > obj.tol)
            for l = k+1:N,                                                  if any(abs(psi(k,:).*psi(l,:)) > obj.tol)
                [a,r]       =   obj.getExcitationIndexes(obj.CSB(k,:),obj.CSB(l,:));
                if isscalar(a),                                                                 if a > 0 && r > 0
                    rUp( a, r,:)    =   squeeze(rUp( a, r,:)) + (conj(psi(k,:)).*psi(l,:)).';
                    rUp( r, a,:)    =   squeeze(rUp( r, a,:)) + (conj(psi(l,:)).*psi(k,:)).';   elseif a < 0 && r < 0
                    rDw(-a,-r,:)    =   squeeze(rDw(-a,-r,:)) + (conj(psi(k,:)).*psi(l,:)).';
                    rDw(-r,-a,:)    =   squeeze(rDw(-r,-a,:)) + (conj(psi(l,:)).*psi(k,:)).';   end
                end,                                                        end
            end,                                                            end
        end
    end
end