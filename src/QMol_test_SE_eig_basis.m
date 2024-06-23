classdef QMol_test_SE_eig_basis < QMol_test
%QMol_test_SE_eig_basis unit test for the QMol_SE_eig_basis class in 1D

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_test)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test_SE_eig_basis\n'); 
    QMol_test_SE_eig_basis.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Real basis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    x                   =   (-15:.1:20).';
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)];
    disc                =   QMol_disc_basis('xspan',x,'basis',v);   disc.orthonormalizeBasis;

    V                   =   QMol_SE_V('atom', ...
                               {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5), ...
                                QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3), ...
                                QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)});

    SE                  =   QMol_SE(                            ...
                                'disc',                 disc,   ...
                                'numberWaveFunction',   3,      ...
                                'potential',            V);

    ES                  =   QMol_SE_eig_basis('display',false);
    ES.computeGroundState(SE);
    E                   =   SE.getEnergy('wfcn');

    R_size              =   all(size(SE.wfcn.wfcn) == [disc.basisSize 3]);

    R_norm              =   abs( sum(abs(SE.wfcn.wfcn(:,1)).^2) - 1) < 1e-10  &&  ... % wfcn are normalized
                            abs( sum(abs(SE.wfcn.wfcn(:,2)).^2) - 1) < 1e-10  &&  ...
                            abs( sum(abs(SE.wfcn.wfcn(:,3)).^2) - 1) < 1e-10;

    Hp                  =   SE.disc.SE_operatorHamiltonian(V,SE.wfcn.wfcn(:,1));
    E1                  =   sum( SE.wfcn.wfcn(:,1) .* Hp );
    DE1                 =   sqrt(sum( abs(Hp - E(1)*SE.wfcn.wfcn(:,1)).^2 ));

    Hp                  =   SE.disc.SE_operatorHamiltonian(V,SE.wfcn.wfcn(:,2));
    E2                  =   sum( SE.wfcn.wfcn(:,2) .* Hp );
    DE2                 =   sqrt(sum( abs(Hp - E(2)*SE.wfcn.wfcn(:,2)).^2 ));

    Hp                  =   SE.disc.SE_operatorHamiltonian(V,SE.wfcn.wfcn(:,3));
    E3                  =   sum( SE.wfcn.wfcn(:,3) .* Hp );
    DE3                 =   sqrt(sum( abs(Hp - E(3)*SE.wfcn.wfcn(:,3)).^2 ));

    R_energy            =   max(abs(E - [E1;E2;E3])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3]) < 1e-10;
    
    obj.showResult('computeGroundState (real basis)',R_size && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of wave functions\n'); end
    if ~R_norm,     fprintf('      - Wave functions are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Wave functions are not converged\n'); end

    % Complex basis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    v                   =   [exp(-(x-2).^2),exp(-(x-1).^2/.7),exp(-x.^2),exp(-(x+1).^2/2)] .* ...
                            [exp( 2i*x),    exp( 3i*x),       exp(-4i*x),exp(-2i*x)];
    v                   =   [v, exp(-(x-4).^2/2.5)];
    disc.set('basis',v);    disc.orthonormalizeBasis;
    SE.reset();             %SE.initialize;

    ES.computeGroundState(SE);
    E                   =   SE.getEnergy('wfcn');

    R_size              =   all(size(SE.wfcn.wfcn) == [disc.basisSize 3]);

    R_norm              =   abs( sum(abs(SE.wfcn.wfcn(:,1)).^2) - 1) < 1e-10  &&  ... % wfcn are normalized
                            abs( sum(abs(SE.wfcn.wfcn(:,2)).^2) - 1) < 1e-10  &&  ...
                            abs( sum(abs(SE.wfcn.wfcn(:,3)).^2) - 1) < 1e-10;

    Hp                  =   SE.disc.SE_operatorHamiltonian(V,SE.wfcn.wfcn(:,1));
    E1                  =   sum( conj(SE.wfcn.wfcn(:,1)) .* Hp );
    DE1                 =   sqrt(sum( abs(Hp - E(1)*SE.wfcn.wfcn(:,1)).^2 ));

    Hp                  =   SE.disc.SE_operatorHamiltonian(V,SE.wfcn.wfcn(:,2));
    E2                  =   sum( conj(SE.wfcn.wfcn(:,2)) .* Hp );
    DE2                 =   sqrt(sum( abs(Hp - E(2)*SE.wfcn.wfcn(:,2)).^2 ));

    Hp                  =   SE.disc.SE_operatorHamiltonian(V,SE.wfcn.wfcn(:,3));
    E3                  =   sum( conj(SE.wfcn.wfcn(:,3)) .* Hp );
    DE3                 =   sqrt(sum( abs(Hp - E(3)*SE.wfcn.wfcn(:,3)).^2 ));

    R_energy            =   max(abs(E - [E1;E2;E3])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3]) < 1e-10;

    obj.showResult('computeEigenstates (complex basis)',R_size && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of Kohn-Sham orbitals\n'); end
    if ~R_norm,     fprintf('      - Kohn-Sham orbitals are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Kohn-Sham orbitals are not converged\n'); end
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

