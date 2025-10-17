function test_setConfigurationBasis_CIS(obj)

    % Initialization ======================================================
    obj.showSection('Configuration state basis with CIS (setConfigurationBasis)');

    % Default parameters ==================================================
    x                   =   1:5;
    SOB                 =   rand(8,5);
    Vext                =   QMol_DFT_Vext;
    CI                  =   QMol_CI_conv('display', false,  ...
                                'numberElectron',   6,      ...
                                'xspan',            x,      ...
                                'orbitalBasis',     SOB,    ...
                                'externalPotential',Vext    );
    CI.setConfigurationBasis;

    % Reference
    ref                 =   CI.configurationBasis(1,:);
    ref_th              =   [-1 -2 -3 1 2 3];
    obj.showResult('reference (default)',all(sort(ref) == sort(ref_th)));

    % (single) excitations
    CSB                 =   unique(sort(CI.configurationBasis(2:end,:),2),'rows');
    CSB_th              =   unique(sort([-1 -2 -4 1 2 3; -1 -2 -3 1 2 4; -1 -2 -5 1 2 3; -1 -2 -3 1 2 5; ...
                                         -1 -4 -3 1 2 3; -1 -2 -3 1 4 3; -1 -5 -3 1 2 3; -1 -2 -3 1 5 3; ...
                                         -4 -2 -3 1 2 3; -1 -2 -3 4 2 3; -5 -2 -3 1 2 3; -1 -2 -3 5 2 3],2),'rows');
    obj.showResult('configurationBasis (default)', all(size(CSB) == size(CSB_th)) && all(ref == ref_th,'all'));

    % User-defined parameters =============================================
    CI.set('reference',[-1 -2 -3 1 2],'active',[-3 -5 -6 1 3 4 5]);
    CI.setConfigurationBasis;

    % Number of electrons
    obj.showResult('numberElectron (from reference)',CI.numberElectron == 5);

    % Configuration basis
    CSB                 =   unique(sort(CI.configurationBasis,2),'rows');
    CSB_th              =   unique(sort([-1 -2 -3 1 2;                      ... reference
                                         -5 -2 -3 1 2; -6 -2 -3 1 2;        ... down spin
                                         -1 -5 -3 1 2; -1 -6 -3 1 2;        ...
                                         -1 -2 -5 1 2; -1 -2 -6 1 2;        ...
                                         -1 -2 -3 3 2; -1 -2 -3 4 2; -1 -2 -3 5 2; ... up spin
                                         -1 -2 -3 1 3; -1 -2 -3 1 4; -1 -2 -3 1 5],2),'rows');
    obj.showResult('configurationBasis (user defined)', all(size(CSB) == size(CSB_th)) && all(CSB == CSB_th,'all'));

    % Freezing orbitals 
    CI.set('frozen',[-2 1]);
    CI.setConfigurationBasis;

    CSB                 =   unique(sort(CI.configurationBasis,2),'rows');
    CSB_th              =   unique(sort([-1 -2 -3 1 2;                      ... reference
                                         -5 -2 -3 1 2; -6 -2 -3 1 2;        ... down spin
                                         -1 -2 -5 1 2; -1 -2 -6 1 2;        ...
                                         -1 -2 -3 1 3; -1 -2 -3 1 4; -1 -2 -3 1 5],2),'rows'); % spin down
    obj.showResult('configurationBasis (with frozen)', all(size(CSB) == size(CSB_th)) && all(CSB == CSB_th,'all'));

    % noDouble constraint 
    CI.set('frozen',[],'noDouble',[2 3]);
    CI.setConfigurationBasis;
    CSB                 =   unique(sort(CI.configurationBasis,2),'rows');
    CSB_th              =   unique(sort([-1 -2 -3 1 2;                      ... (noDouble does not apply to the reference)
                                         -1 -5 -3 1 2; -1 -6 -3 1 2;        ...
                                         -1 -2 -3 1 4; -1 -2 -3 1 5],2),'rows');
    obj.showResult('configurationBasis (with noDouble)', all(size(CSB) == size(CSB_th)) && all(CSB == CSB_th,'all'));

    % noEmpty constraint 
    CI.set('noEmpty',[3 5]);
    CI.setConfigurationBasis;
    CSB                 =   unique(sort(CI.configurationBasis,2),'rows');
    CSB_th              =   unique(sort([-1 -2 -3 1 2;                      ... (noEmpty does not apply to the reference)
                                         -1 -5 -3 1 2; -1 -2 -3 1 5],2),'rows');
    obj.showResult('configurationBasis (with noEmpty)', all(size(CSB) == size(CSB_th)) && all(CSB == CSB_th,'all'));

    % Multi-reference =====================================================
    CI.set('reference',[-1 -2 -3 1 2; -1 -2 -4 1 2],'active',[-3 -5 -6 1 3 4 5], ...
        'frozen',[-1 -2 1],'noDouble',[],'noEmpty',[]);
    CI.setConfigurationBasis;

    CSB                 =   unique(sort(CI.configurationBasis,2),'rows');
    CSB_th              =   unique(sort([-1 -2 -3 1 2; -1 -2 -4 1 2;        ... reference
                                         -1 -2 -5 1 2; -1 -2 -6 1 2;        ... down spin
                                         -1 -2 -3 1 3; -1 -2 -3 1 4; -1 -2 -3 1 5;        ... down spin
                                         -1 -2 -4 1 3; -1 -2 -4 1 4; -1 -2 -4 1 5],2),'rows');
    obj.showResult('configurationBasis (multi-reference)', all(size(CSB) == size(CSB_th)) && all(CSB == CSB_th,'all'));

end

