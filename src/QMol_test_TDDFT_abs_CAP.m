classdef QMol_test_TDDFT_abs_CAP < QMol_test
%QMol_test_TDDFT_abs_CAP suite of unit tests for QMol_TDDFT_abs_CAP

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
    fprintf('  * QMol_test_TDDFT_abs_CAP\n'); 
    QMol_test_TDDFT_abs_CAP.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test components
    obj.runRunTime;
    obj.runApplyGobbler;
    
end
end
%% Test components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=private)
function runRunTime(obj)
%runTest unit tests for the class

    % Initialization
    obj.showSection('Run time variables');
    x                   =   -15:.1:20;
    
    disc                =   QMol_disc('x',x);
    DFT                 =   QMol_DFT_spinRes('disc',disc);
    disc.initialize(disc);

    % Create mask
    abc                 =   QMol_TDDFT_abs_CAP('length',[exp(2) log(100)]);
    abc.initialize(DFT,true);

    R                   =   all(size(abc.V) == [numel(x) 1]);
    obj.showResult('CAP size' ,R);
    if ~R, fprintf('      SKIPPING THE OTHER TESTS\n'); return,   end

    % sin shape
    abc.set('shape','sin','amplitude',.8);
    
    abc.initialize(DFT,true);
    X                   =   (x(13)+15)/exp(2);
    R                   =   abs(abc.V(13)+1i*.8*(1-sin(.5*pi*X))) < 1e-10;

    abc.initialize(DFT,false);
    X                   =   (20-x(end-20))/log(100);
    R                   =   R && abs(abc.V(end-20)-1i*.8*(1-sin(.5*pi*X))) < 1e-10;
    obj.showResult('sin shape' ,R);

    % sin^2 shape
    abc.set('shape','sin^2','amplitude',.8);
    
    abc.initialize(DFT,true);
    X                   =   (x(13)+15)/exp(2);
    R                   =   abs(abc.V(13)+1i*.8*(1-sin(.5*pi*X)^2)) < 1e-10;

    abc.initialize(DFT,false);
    X                   =   (20-x(end-20))/log(100);
    R                   =   R && abs(abc.V(end-20)-1i*.8*(1-sin(.5*pi*X)^2)) < 1e-10;
    obj.showResult('sin^2 shape' ,R);

    % sin^1/8 shape
    abc.set('shape','sin^1/8','amplitude',.8);
    
    abc.initialize(DFT,true);
    X                   =   (x(13)+15)/exp(2);
    R                   =   abs(abc.V(13)+1i*.8*(1-sin(.5*pi*X)^.125)) < 1e-10;

    abc.initialize(DFT,false);
    X                   =   (20-x(end-20))/log(100);
    R                   =   R && abs(abc.V(end-20)-1i*.8*(1-sin(.5*pi*X)^.125)) < 1e-10;
    obj.showResult('sin^1/8 shape' ,R);

    % user defined shape
    abc.set('shape',@(x) sin(.5*pi*x).^3,'amplitude',.8);
    
    abc.initialize(DFT,true);
    X                   =   (x(13)+15)/exp(2);
    R                   =   abs(abc.V(13)+1i*.8*cos(.5*pi*X)^3) < 1e-10;

    abc.initialize(DFT,false);
    X                   =   (20-x(end-20))/log(100);
    R                   =   R && abs(abc.V(end-20)-1i*.8*cos(.5*pi*X)^3) < 1e-10;
    obj.showResult('user defined shape' ,R);

end
function runApplyGobbler(obj)
%runApplyGobbler unit tests for the application of the gobbler
    
    % Initialization
    obj.showSection('Absorbing boundaries');

    % Spin restricted
    disc                =   QMol_disc('xspan',-15:.1:20);
    DFT                 =   QMol_DFT_spinRes('disc',disc,'occupation',[1 2 3]);
    disc.initialize(DFT);

    Vks                 =   disc.DFT_allocatePotential;

    abc                 =   QMol_TDDFT_abs_CAP('length',[8 8]);
    abc.initialize(DFT,true);

    abc.getPotential(Vks);
    R                   =   max(abs(Vks.V-abc.V) < 1e-10);
    obj.showResult('getPotential (spin restricted)' ,R);

    % Spin polarized
    DFT                 =   QMol_DFT_spinPol('disc',disc,'occupation',{[1 2 3],[4 5]});
    disc.initialize(DFT);

    disc.DFT_allocatePotential(Vks);

    abc.initialize(DFT,false);

    abc.getPotential(Vks);
    R                   =   max(abs(Vks.Vup-abc.V) < 1e-10) && ...
                            max(abs(Vks.Vdw-abc.V) < 1e-10);                %#ok<LOGMAX>
    obj.showResult('getPotential (spin polarized)' ,R);

    % applyMask
    abc.applyMask([]);      % Will cause an error if does something
    obj.showResult('applyMask',true);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

