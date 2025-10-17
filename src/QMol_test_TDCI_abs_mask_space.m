classdef QMol_test_TDCI_abs_mask_space < QMol_test
%QMol_test_TDCI_abs_CAP_space of unit tests for QMol_TDCI_abs_CAP_space

%   Version     Date        Author
%   01.23.000   06/07/2025  F. Mauger
%       Creation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_test)
function version
    QMol_doc.showVersion('01.23.000','06/07/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_test})
function showInfo
    fprintf('  * QMol_test_TDCI_abs_mask_space\n'); 
    QMol_test_TDCI_abs_mask_space.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class

    % Initialization
    obj.showSection('Mask-type damping matrix (from spatial mask)');

    x                   =   -15:.1:20;
    SOB                 =  [exp(-.5*x.^2);                  x.*exp(-.5*x.^2/2); ...
                            cos(x-2).*exp(-.5*(x-2).^2/3);  sin(x-4).*exp(-.5*(x-4).^2/4)].';
    CSB                 =  [-1 -2 1 2; -1 -3 1 2; -1 -4 1 2;  -1 -2 1 3; -1 -4 1 3];
    Vext                =   QMol_DFT_Vext('atom',QMol_Va_softCoulomb);

    CI                  =   QMol_CI_conv(...
                                'discretization',       QMol_disc('x',x),   ...
                                'numberElectron',       4,                  ...
                                'orbitalBasis',         SOB,                ...
                                'configurationBasis',   CSB,                ...
                                'externalPotential',    Vext                );
    CI.initialize;

    % Create CAP
    dx                  =   x(2)-x(1);

    function V = getMask(x,D,L,fcn)
        % Initialization
        V               =   ones(numel(x),1);

        % Negative positions
        ind             =   x < -D(1);          V(ind)  =   fcn(abs(x(ind)+D(1))/L(1));
        ind             =   x <=-(D(1)+L(1));   V(ind)  =   0;

        % Positive positions
        ind             =   x >  D(2);          V(ind)  =   fcn(abs(x(ind)-D(2))/L(2));
        ind             =   x >= D(2)+L(2);     V(ind)  =   0;
    end

    function M = getDampingMatrix(V)
        % Initialization
        M               =   eye(size(CSB,1),size(CSB,1));

        % Build damping matrix
        for k = 1:size(CSB,1)
            for l = 1:size(CSB,2)
                M(k,k)  =   M(k,k) * sum(V.* (SOB(:,abs(CSB(k,l))).^2) )*dx;
            end
        end
    end

    % Test sin^1/8
    D                   =   2;
    L                   =   5;
    abc                 =   QMol_TDCI_abs_mask_space('distance',D,'length',L);
    abc.initialize(CI);

    V                   =   getMask(x,[D D],[L L],@(x) cos(.5*pi*x).^.125);
    M                   =   abc.getDampingMatrix - getDampingMatrix(V);
    obj.showResult('sin^1/8 shape' ,all(abs(M)<1e-10,'all'));

    % Test sin^2
    D                   =   [1 3];
    L                   =   5;
    abc.set('distance',D,'length',L,'shape','sin^2');
    abc.initialize(CI);

    V                   =   getMask(x,D,[L L],@(x) 1-sin(.5*pi*x).^2);
    M                   =   abc.getDampingMatrix - getDampingMatrix(V);
    obj.showResult('sin^2 shape' ,all(abs(M)<1e-10,'all'));

    % Test sin
    D                   =   2.5;
    L                   =   [3 8];
    abc.set('distance',D,'length',L,'shape','sin');
    abc.initialize(CI);

    V                   =   getMask(x,[D D],L,@(x) cos(.5*pi*x));
    M                   =   abc.getDampingMatrix - getDampingMatrix(V);
    obj.showResult('sin shape' ,all(abs(M)<1e-10,'all'));

    % Test function handle
    D                   =   [3.5 5];
    L                   =   [3 8];
    abc.set('distance',D,'length',L,'shape',@(x) 1-x);
    abc.initialize(CI);

    V                   =   getMask(x,D,L,@(x) 1-x);
    M                   =   abc.getDampingMatrix - getDampingMatrix(V);
    obj.showResult('function handle (user-defined shape)' ,all(abs(M)<1e-10,'all'));
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

