function test_initialization(obj)

    % Initialization
    obj.showSection('Initialization');

    % Common parameters
    x                   =   -20:.1:15;

    SOB                 =   rand(numel(x),5);
    CSB                 =   [-1 -2 1 2 3;-1 -3 1 2 3];

    % User-defined CI model
    Ve                  =   @(x) x+3;
    CI_user             =   QMol_CI_conv('display',     false,              ...
                                'numberElectron',       5,                  ...
                                'xspan',                x,                  ...
                                'orbitalBasis',         SOB,                ...
                                'configurationBasis',   CSB,                ...
                                'externalPotential',    QMol_DFT_Vext,      ...
                                'interactionPotential', Ve                  );
    CI_user.initialize;

    % CI model from DFT
    Va                  =  {QMol_Va_Gaussian('V0',1,'s',.5,'X0',-3), ...
                            QMol_Va_Gaussian('V0',sqrt(2),'s',exp(1),'X0',log(5)),...
                            QMol_Va_softCoulomb('Z',1/3,'a',pi,'X0',5)};
    V                   =   @(x) -3*x./sqrt((x-sqrt(2)).^2 + exp(2));
    DV                  =   @(x) -3./sqrt((x-sqrt(2)).^2 + exp(2)) + 3*x.*(x-sqrt(2)).*((x-sqrt(2)).^2 + exp(2)).^-1.5;
    Vext                =   QMol_DFT_Vext('atom',Va,'externalPotential',V,'externalPotentialDerivative',DV);
    Vh                  =   QMol_DFT_Vh_conv;
    Vx                  =   QMol_DFT_Vx_LDA_exp;
    DFT                 =   QMol_DFT_spinRes( ...
                                'xspan',                        x,          ...
                                'occupation',                   [2 2 1],    ...
                                'externalPotential',            Vext,       ...
                                'HartreePotential',             Vh,         ...
                                'exchangeCorrelationPotential', Vx          );
    DFT.initialize;
    CI_DFT              =   QMol_CI_conv(DFT,'display', false,  ...
                                'orbitalBasis',         SOB,    ...
                                'configurationBasis',   CSB);
    CI_DFT.initialize;

    % discretization 
    obj.showResult('discretization (user defined)',CI_user.discretization == QMol_disc('xspan',x));
    obj.showResult('discretization (from DFT/HF)', CI_DFT.discretization  == QMol_disc('xspan',x));

    % Total charge
    obj.showResult('numberElectron (user defined)',CI_user.numberElectron == 5);
    obj.showResult('numberElectron (from DFT/HF)', CI_DFT.numberElectron  == 5);

    % isInitialized
    obj.showResult('isInitialized',CI_user.isInitialized && CI_DFT.isInitialized);

    % orbitalBasis
    obj.showResult('orbitalBasis',all(CI_user.orbitalBasis == SOB & CI_DFT.orbitalBasis == SOB));

    % configurationBasis
    obj.showResult('configurationBasis',all(CI_user.configurationBasis == CSB & CI_DFT.configurationBasis == CSB));

    % externalPotential
    Vext.clear;             % CI should have made a local copy
    Vc                  =   CI_DFT.externalPotential;
    isOK                =   isequal(Vc.Vext,V) & isequal(Vc.DVext,DV) & numel(Vc.atom) == 3;
    isOK                =   isOK && Vc.atom{1}.V0 == Va{1}.V0 && Vc.atom{1}.s == Va{1}.s && Vc.atom{1}.X0 == Va{1}.X0;
    isOK                =   isOK && Vc.atom{2}.V0 == Va{2}.V0 && Vc.atom{2}.s == Va{2}.s && Vc.atom{2}.X0 == Va{2}.X0;
    isOK                =   isOK && Vc.atom{3}.Z  == Va{3}.Z  && Vc.atom{3}.a == Va{3}.a && Vc.atom{3}.X0 == Va{3}.X0;
    
    obj.showResult('externalPotential (from DFT/HF)',isOK);

    % interactionPotential
    obj.showResult('interactionPotential (user defined)',isequal(CI_user.interactionPotential,Ve));
    obj.showResult('interactionPotential (from DFT/HF)',isequal(CI_DFT.interactionPotential,Vh.interactionPotential));

end

