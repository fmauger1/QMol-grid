function proj = DFT_projectOrbital(obj,KSO,proj)
%DFT_projectOrbital computes the projection coefficients associated with
%    the (discretized) KS orbitals

    
    % Initialization
    if nargin == 2,     proj    =   [];                     end
    if isempty(proj),   proj    =   QMol_DFT_orbital_basis; end
    
    if isempty(KSO.KSO)     % Spin polarized
        % Initialize components
        proj.isSpinPol  =   true;
        proj.KSO        =   [];
        proj.KSOup      =   NaN(obj.nV,size(KSO.KSOup,2));
        proj.KSOdw      =   NaN(obj.nV,size(KSO.KSOdw,2));

        % Compute projections
        for k = 1:obj.nV, if isreal(obj.v)                                  %#ok<ALIGN> 
            for l = 1:size(KSO.KSOup,2),    proj.KSOup(k,l) =   sum(     obj.v(:,k) .*KSO.KSOup(:,l)) * obj.dx; end
            for l = 1:size(KSO.KSOdw,2),    proj.KSOdw(k,l) =   sum(     obj.v(:,k) .*KSO.KSOdw(:,l)) * obj.dx; end, else
            for l = 1:size(KSO.KSOup,2),    proj.KSOup(k,l) =   sum(conj(obj.v(:,k)).*KSO.KSOup(:,l)) * obj.dx; end
            for l = 1:size(KSO.KSOdw,2),    proj.KSOdw(k,l) =   sum(conj(obj.v(:,k)).*KSO.KSOdw(:,l)) * obj.dx; end, end
        end
    else                    % Spin restricted
        % Initialize components
        proj.isSpinPol  =   false;
        proj.KSO        =   NaN(obj.nV,size(KSO.KSO,2));
        proj.KSOup      =   [];             proj.KSOdw      =   [];

        % Compute projections
        for k = 1:obj.nV, if isreal(obj.v)                                  %#ok<ALIGN> 
            for l = 1:size(KSO.KSO,2),      proj.KSO(k,l)   =   sum(     obj.v(:,k) .*KSO.KSO(:,l)) * obj.dx;   end, else
            for l = 1:size(KSO.KSO,2),      proj.KSO(k,l)   =   sum(conj(obj.v(:,k)).*KSO.KSO(:,l)) * obj.dx;   end, end
        end
    end

    % Initialize (mirror KS orbital input)
    if KSO.isInit,      proj.initialize(obj);        end

end