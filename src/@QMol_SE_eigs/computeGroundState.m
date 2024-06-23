function computeGroundState(obj,SE)
%computeGroundState computes the ground state(s) for the input Schrodinger-
%   equation model
    
    % Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.initialize(SE);                                                     % resets the object first
    obj.isRun   =   true;
    
    if obj.dispGS,  QMol_doc.showHeader;        end
    
    % (Re)link everything ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.dispGS
        QMol_doc.showSection('Build Schrodinger-equation (SE) model');
        SE.initialize(2);
    else
        SE.initialize(0);
    end
    
    if obj.dispGS,  obj.showDocumentation;  fprintf('\n');  end

    % Compute eigen states ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.computeEigenstates;

    % Finalization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if obj.dispGS
        
        % Show wave function and Schrodinger-equation energies
        QMol_doc.showSection('Wave-function energies');
        SE.showEnergy('wfcn');      fprintf('\n');
        
        QMol_doc.showSection('Schrodinger-equation-component energies');
        SE.showEnergy('SE');        fprintf('\n');
        
        % Footer
        QMol_doc.showFooter;
    end
    
    % Clear member temporary variables
    obj.isRun           =   false;
end

