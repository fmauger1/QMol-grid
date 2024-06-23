classdef QMol_test_TDSE_abs_mask < QMol_test
%QMol_test_TDSE_abs_mask suite of unit tests for QMol_TDSE_abs_mask

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
    fprintf('  * QMol_test_TDSE_abs_mask\n'); 
    QMol_test_TDSE_abs_mask.version;
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
    SE                  =   QMol_SE('disc',disc);
    disc.initialize(disc);

    % Create mask
    abc                 =   QMol_TDSE_abs_mask('length',[exp(2) log(100)]);
    abc.initialize(SE);

    R                   =   all(size(abc.W) == [numel(x) 1]);
    obj.showResult('Mask size' ,R);
    if ~R, fprintf('      SKIPPING THE OTHER TESTS\n'); return,   end

    % cos^2 shape
    abc.set('shape','cos^2');   abc.initialize(SE);

    X                   =   (x(13)+15)/exp(2);
    R                   =   abs(abc.W(13)-sin(.5*pi*X)^2) < 1e-10;

    X                   =   (20-x(end-20))/log(100);
    R                   =   R && abs(abc.W(end-20)-sin(.5*pi*X)^2) < 1e-10;
    obj.showResult('cos^2 shape' ,R);

    % cos^1/8 shape
    abc.set('shape','cos^1/8');   abc.initialize(SE);

    X                   =   (x(13)+15)/exp(2);
    R                   =   abs(abc.W(13)-sin(.5*pi*X)^.125) < 1e-10;

    X                   =   (20-x(end-20))/log(100);
    R                   =   R && abs(abc.W(end-20)-sin(.5*pi*X)^.125) < 1e-10;
    obj.showResult('cos^1/8 shape' ,R);

    % user-defined shape
    abc.set('shape',@(x) cos(.5*pi*x));   abc.initialize(SE);

    X                   =   (x(13)+15)/exp(2);
    R                   =   abs(abc.W(13)-sin(.5*pi*X)) < 1e-10;

    X                   =   (20-x(end-20))/log(100);
    R                   =   R && abs(abc.W(end-20)-sin(.5*pi*X)) < 1e-10;
    obj.showResult('user-defined shape' ,R);

end
function runApplyGobbler(obj)
%runApplyGobbler unit tests for the application of the gobbler
    
    % Initialization
    obj.showSection('Absorbing boundaries');

    % Apply mask
    disc                =   QMol_disc('xspan',-15:.1:20);
    SE                  =   QMol_SE('disc',disc,'numberWaveFunction',3);
    disc.initialize(SE);
    
    Wfcn                =   disc.SE_allocateWaveFunction(3);
    Wfcn.wfcn(:,:)      =   1;
    SE.set('Wfcn',Wfcn);
    
    abc                 =   QMol_TDSE_abs_mask('length',[8 8]);
    abc.initialize(SE);
    abc.applyMask(Wfcn);

    R                   =   max(abs(Wfcn.wfcn(:,1)-abc.W) < 1e-10) && ...
                            max(abs(Wfcn.wfcn(:,2)-abc.W) < 1e-10) && ...
                            max(abs(Wfcn.wfcn(:,3)-abc.W) < 1e-10);         %#ok<LOGMAX>
    obj.showResult('applyMask' ,R);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

