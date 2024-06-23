function DVks = DFT_allocatePotentialGradient(obj,DVks)
%DFT_allocatePotentialGradient allocates a DFT-potential gradient operator 
%   object to store the Kohn-Sham, external, Hartree, and/or exchange 
%   correlation potential gradients
    
    % Initialization
    if nargin == 1, DVks=   [];     end
    
    % Allocate/reinitialized density object (allocate to zero)
    if isempty(DVks)
        DVks            =   QMol_DFT_Vks_grad('isSpinPol',obj.QM.isSpinPol);
    else
        DVks.clear;
        DVks.set('isSpinPol',obj.QM.isSpinPol);
    end

    % Do not initialize on allocate
end