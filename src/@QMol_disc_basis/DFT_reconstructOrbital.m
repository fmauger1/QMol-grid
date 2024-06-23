function KSO = DFT_reconstructOrbital(obj,proj,KSO)
%DFT_reconstructOrbital reconstructs the (discretization) of the KS
%   orbitals associated with the proj projection coefficients
    
    % Initialization
    if nargin == 2,     KSO     =   [];     end
    
    if obj.QM.isSpinPol     % Spin polarized
        if isempty(KSO)
            KSO         =   QMol_DFT_orbital('isSpinPol',true,      ...
                                'KSOup',zeros(numel(obj.x),size(proj.KSOup,2)),   ...
                                'KSOdw',zeros(numel(obj.x),size(proj.KSOdw,2)));
        else
            KSO.KSO     =   [];
            KSO.KSOup   =   zeros(numel(obj.x),size(proj.KSOup,2));
            KSO.KSOdw   =   zeros(numel(obj.x),size(proj.KSOdw,2));
        end
    else                    % Spin restricted
        if isempty(KSO)
            KSO         =   QMol_DFT_orbital('isSpinPol',false,      ...
                                'KSO',zeros(numel(obj.x),size(proj.KSO,2)));
        else
            KSO.KSO     =   zeros(numel(obj.x),size(proj.KSO,2));
            KSO.KSOup   =   [];         KSO.KSOdw   =   [];
        end
    end

    % Reconstruct wave function
    if obj.QM.isSpinPol, for k = 1:obj.nV, for l = 1:size(proj.KSOup,2)     %#ok<ALIGN> 
        KSO.KSOup(:,l)  =   KSO.KSOup(:,l) + proj.KSOup(k,l)*obj.v(:,k);   end, for l = 1:size(proj.KSOdw,2)
        KSO.KSOdw(:,l)  =   KSO.KSOdw(:,l) + proj.KSOdw(k,l)*obj.v(:,k);   end, end
    else, for k = 1:obj.nV, for l = 1:size(proj.KSO,2)
        KSO.KSO(:,l)    =   KSO.KSO(:,l)   + proj.KSO(k,l)*obj.v(:,k);     end, end
    end

    % Initialize (mirror projection input)
    if proj.isInit,     KSO.initialize(obj);         end

end