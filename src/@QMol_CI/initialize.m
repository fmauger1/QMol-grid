function initialize(obj,showBuild)
%initialize initializes the CI object and all its components
    
    % Initialization
    if nargin == 1,         showBuild   =   [];     end
    if isempty(showBuild),  showBuild   =   0;      end

    function buildSection(txt)
        if showBuild > 0,   fprintf('  * %-66s',txt);   end
    end
    function buildStatus(t)
        if showBuild > 0,   if t, fprintf('  OK\n'); else, fprintf('FAIL\n'); end, end
    end

    if showBuild == 1
        QMol_doc.showHeader;
        QMol_doc.showFunding;
        QMol_doc.showSection('Building:');
    end

    obj.reset;
    
    % Discretization
    buildSection('Discretization');
    if isempty(obj.disc)
        warning('QMol:CI:noDiscretization','Missing the discretization component; initialization skipped on it.')
        buildStatus(false);
    else
        % Check if discretization has changed
        if ~obj.disc.isInit,        obj.reset;      end

        % Initialize the discretization component
        obj.disc.initialize(obj);
        buildStatus(obj.disc.isInit);
    end
    
    % Atomic or molecular potential
    buildSection('External (atomic or molecular) potential');
    if isempty(obj.Vext)
        warning('QMol:CI:noAtomicMolecularPotential','Missing the atomic or molecular potential component; initialization skipped it.')
        buildStatus(false)
    else
        % Implicit definition of the potential
        if ~isa(obj.Vext,'QMol_DFT_Vext'),                                  if isnumeric(obj.Vext) || isa(obj.Vext, 'function_handle')
            obj.Vext    =   QMol_DFT_Vext('Vext',obj.Vext);                 elseif iscell(obj.Vext)
            obj.Vext    =   QMol_DFT_Vext('atom',obj.Vext);                 else
            obj.Vext    =   QMol_DFT_Vext('atom',{obj.Vext});               end
        end

        % Initialize the potential
        obj.Vext.initialize(obj);
        buildStatus(obj.Vext.isInit);
    end

    % Electron-electron interaction potential
    buildSection('Electron-electron interaction potential');
    buildStatus(obj.initializeVee);

    % Miscellaneous
    obj.isReal          =   isreal(obj.SOB);

    % House keeping
    obj.isInit          =   true;

    if showBuild == 1,  QMol_doc.showFooter;      end
end