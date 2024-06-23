function Vks = DFT_allocatePotential(obj,Vks)
%DFT_allocatePotential allocates a DFT-potential operator object to store
%   the Kohn-Sham, external, Hartree, and/or exchange correlation
%   potentials
    
    % Initialization
    if nargin == 1, Vks =   [];     end
    
    % Allocate/reinitialized density object (allocate to zero)
    if isempty(Vks)
        Vks             =   QMol_DFT_Vks_basis('isSpinPol',obj.QM.isSpinPol);
    else
        Vks.clear;
        Vks.set('isSpinPol',obj.QM.isSpinPol);
    end

    % Do not initialize on allocate
end