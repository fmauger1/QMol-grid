classdef QMol_test_SE_eigs < QMol_test
%QMol_test_SE_eigs suite of unit tests for QMol_SE_eigs

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
    fprintf('  * QMol_test_SE_eigs\n'); 
    QMol_test_SE_eigs.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Object initialization
    obj.test_groundState;
end
end
methods (Access=private)
function test_groundState(obj)

    % No symmetry imposed =================================================
    Vse                 =   QMol_SE_V('atom', ...
                               {QMol_Va_softCoulomb('name','(1)','Z',1,'X0',-1.5), ...
                                QMol_Va_softCoulomb('name','(2)','Z',1,'X0', 1/3), ...
                                QMol_Va_softCoulomb('name','(3)','Z',2,'X0', pi)});

    SE                  =   QMol_SE(...
                                'xspan',                   -15:.1:20,       ...
                                'numberWaveFunction',       3,              ...
                                'potential',                Vse);

    ES                  =   QMol_SE_eigs('display','off');
    ES.computeGroundState(SE);
    E                   =   SE.getEnergy('wfcn');

    x                   =   SE.x(:);
    wfcn                =   SE.wfcn;

    R_size              =   all(size(wfcn.wfcn) == [numel(x) 3]);

    dx                  =   x(2)-x(1);
    R_norm              =   abs( sum(abs(wfcn.wfcn(:,1)).^2)*dx - 1) < 1e-10  &&  ... % wfcn are normalized
                            abs( sum(abs(wfcn.wfcn(:,2)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(wfcn.wfcn(:,3)).^2)*dx - 1) < 1e-10;

    Hp                  =   SE.disc.SE_operatorHamiltonian(SE.V,wfcn.wfcn(:,1));
    E1                  =   sum( wfcn.wfcn(:,1) .* Hp ) * dx;
    DE1                 =   sqrt(sum( abs(Hp - E(1)*wfcn.wfcn(:,1)).^2 ) *dx);

    Hp                  =   SE.disc.SE_operatorHamiltonian(SE.V,wfcn.wfcn(:,2));
    E2                  =   sum( wfcn.wfcn(:,2) .* Hp ) * dx;
    DE2                 =   sqrt(sum( abs(Hp - E(2)*wfcn.wfcn(:,2)).^2 ) *dx);

    Hp                  =   SE.disc.SE_operatorHamiltonian(SE.V,wfcn.wfcn(:,3));
    E3                  =   sum( wfcn.wfcn(:,3) .* Hp ) * dx;
    DE3                 =   sqrt(sum( abs(Hp - E(3)*wfcn.wfcn(:,3)).^2 ) *dx);

    R_energy            =   max(abs(E - [E1;E2;E3])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3]) < 1e-10;
    
    obj.showResult('computeGroundState (no symmetry)',R_size && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of wave functions\n'); end
    if ~R_norm,     fprintf('      - wave functions are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Wave functions are not converged\n'); end

    % Symmetry imposed ====================================================
    SE.set('xspan',-20:.1:20,'numberWaveFunction',4);
    Vse.set('atom',{QMol_Va_softCoulomb('name','(1)','Z',2,'a',.5,'X0',-1.5), ...
                    QMol_Va_softCoulomb('name','(2)','Z',2,'a',.5,'X0', 1.5)});

    ES.set('Symmetry','1 Sx + 3 Ax');
    ES.computeGroundState(SE);
    E                   =   SE.getEnergy('wfcn');

    x                   =   SE.x(:);
    wfcn                =   SE.wfcn;

    R_size              =   all(size(wfcn.wfcn) == [numel(x) 4]);
    
    R_norm              =   abs( sum(abs(wfcn.wfcn(:,1)).^2)*dx - 1) < 1e-10  &&  ...   wfcn are normalized
                            abs( sum(abs(wfcn.wfcn(:,2)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(wfcn.wfcn(:,3)).^2)*dx - 1) < 1e-10  &&  ...
                            abs( sum(abs(wfcn.wfcn(:,4)).^2)*dx - 1) < 1e-10;
    
    R_symmetry          =   max(abs(wfcn.wfcn(:,1)-flip(wfcn.wfcn(:,1)))) < 1e-10 && ... wfcn symmetry
                            max(abs(wfcn.wfcn(:,2)+flip(wfcn.wfcn(:,2)))) < 1e-10 && ...
                            max(abs(wfcn.wfcn(:,3)+flip(wfcn.wfcn(:,3)))) < 1e-10 && ...
                            max(abs(wfcn.wfcn(:,4)+flip(wfcn.wfcn(:,4)))) < 1e-10;

    Hp                  =   SE.disc.SE_operatorHamiltonian(SE.V,wfcn.wfcn(:,1));
    E1                  =   sum( wfcn.wfcn(:,1) .* Hp ) * dx;
    DE1                 =   sqrt(sum( abs(Hp - E(1)*wfcn.wfcn(:,1)).^2 ) *dx);

    Hp                  =   SE.disc.SE_operatorHamiltonian(SE.V,wfcn.wfcn(:,2));
    E2                  =   sum( wfcn.wfcn(:,2) .* Hp ) * dx;
    DE2                 =   sqrt(sum( abs(Hp - E(2)*wfcn.wfcn(:,2)).^2 ) *dx);

    Hp                  =   SE.disc.SE_operatorHamiltonian(SE.V,wfcn.wfcn(:,3));
    E3                  =   sum( wfcn.wfcn(:,3) .* Hp ) * dx;
    DE3                 =   sqrt(sum( abs(Hp - E(3)*wfcn.wfcn(:,3)).^2 ) *dx);

    Hp                  =   SE.disc.SE_operatorHamiltonian(SE.V,wfcn.wfcn(:,4));
    E4                  =   sum( wfcn.wfcn(:,4) .* Hp ) * dx;
    DE4                 =   sqrt(sum( abs(Hp - E(4)*wfcn.wfcn(:,4)).^2 ) *dx);

    R_energy            =   max(abs(E - [E1;E2;E3;E4])) < 1e-10;
    R_error             =   max([DE1 DE2 DE3 DE4]) < 1e-10;

    obj.showResult('computeGroundState (imposed symmetry)',R_size && R_symmetry && R_norm && R_energy && R_error);
    if ~R_size,     fprintf('      - Wrong number of wave functions\n'); end
    if ~R_size,     fprintf('      - Wrong wave function symmetries\n'); end
    if ~R_norm,     fprintf('      - Wave functions are not normalized\n'); end
    if ~R_energy,   fprintf('      - Wrong energy values\n'); end
    if ~R_error,    fprintf('      - Wave functions are not converged\n'); end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

