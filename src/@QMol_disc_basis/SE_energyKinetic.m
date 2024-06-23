function E = SE_energyKinetic(obj,wfcn)
%SE_kineticEnergy computes the kinetic energy of the member SE, and
%   associated wave function, object

    % Parse orbitals
    E                   =   0;

    if isempty(obj.Tv), for k = 1:size(wfcn.wfcn,2)                         %#ok<ALIGN>
        E           =   E + real( wfcn.wfcn(:,k)'*obj.mT *wfcn.wfcn(:,k) );
    end, else, for k = 1:size(wfcn.wfcn,2)                                  %#ok<ALIGN> 
        E           =   E + real( wfcn.wfcn(:,k)'*obj.mTv*wfcn.wfcn(:,k) );
    end, end

end