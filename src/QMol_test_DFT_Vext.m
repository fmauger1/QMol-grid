classdef QMol_test_DFT_Vext < QMol_test
%QMol_test_DFT_Vext suite of unit tests for QMol_DFT_Vext

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
    fprintf('  * QMol_test_DFT_Vext\n'); 
    QMol_test_DFT_Vext.version;
end
end
%% Run tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=?QMol_test)
function testUnit(obj)
%testUnit run all unit tests on the class
    
    % Run test units
    obj.runDefinePotential;
    obj.runUsePotential;
    
end
end
methods (Access=private)
function runDefinePotential(obj) %=========================================
%runDefinePotential unit tests defining the external potential
    
    % Initialization
    obj.showSection('External potential discretization');

    disc                =   QMol_disc('xspan',-15:.1:20);
    DFT                 =   QMol_DFT_spinRes('discretization',disc);
    Vext                =   QMol_DFT_Vext;
    
    disc.initialize(DFT);

    % User-defined function
    V                   =   @(x) -3*x./sqrt((x-sqrt(2)).^2 + exp(2));
    
    Vext.set('externalPotential',V);
    Vext.initialize(DFT);

    obj.showResult('function handle' ,max(abs( Vext.V-V(disc.x(:)) )) < 1e-10);

    % User-defined discretization
    V                   =   V(disc.x(:));
    
    Vext.set('externalPotential',V.');
    Vext.initialize(DFT);
    
    obj.showResult('user-defined discretization' ,max(abs( Vext.V-V )) < 1e-10);
    
    % List of atomic centers
    Va                  =  {QMol_Va_Gaussian('V0',1,'s',.5,'X0',-3), ...
                            QMol_Va_Gaussian('V0',sqrt(2),'s',exp(1),'X0',log(5)),...
                            QMol_Va_softCoulomb('Z',1/3,'a',pi,'X0',5)};
    
    Vext.set('externalPotential',[],'atom',Va);
    Vext.initialize(DFT);
    
    V                   =   Va{1}.getPotential(disc.x(:)) + ...
                            Va{2}.getPotential(disc.x(:)) + ...
                            Va{3}.getPotential(disc.x(:));
    R                   =   max(abs( Vext.V-V )) < 1e-10;

    obj.showResult('list of atomic centers' ,R);
    if ~R, fprintf('      (check QMol_Va_Gaussian and QMol_Va_softCoulomb components)\n');    end

    % Both user-defined and list of atoms
    V                   =   @(x) -3*x./sqrt((x-sqrt(2)).^2 + exp(2));

    Vext.set('externalPotential',V,'atom',Va{1});
    Vext.initialize(DFT);

    V                   =   V(disc.x(:)) + Va{1}.getPotential(disc.x(:));
    R                   =   max(abs( Vext.V-V )) < 1e-10;

    obj.showResult('user-defined + list of atomic centers' ,R);
    if ~R, fprintf('      (check QMol_Va_Gaussian and QMol_Va_softCoulomb components)\n');    end

    % Potential derivative ~~~~~~~~~~~~~~~~~
    obj.showSection('External potential derivative discretization');

    V                   =   @(x) -3*x./sqrt((x-sqrt(2)).^2 + exp(2));
    DV                  =   @(x) -3./sqrt((x-sqrt(2)).^2 + exp(2)) + 3*x.*(x-sqrt(2)).*((x-sqrt(2)).^2 + exp(2)).^-1.5;

    Vext.set('externalPotential',V,'atom',[]);
    Vext.initialize(DFT);   Vext.getPotentialDerivative(1);                 % Store result in Vext.DV
    obj.showResult('finite difference approximation',max(abs(Vext.DV-DV(disc.x(:)))) < 5e-10);
    
    Vext.set('externalPotentialDerivative',DV);
    Vext.initialize(DFT);   Vext.getPotentialDerivative(1);                 % Store result in Vext.DV
    obj.showResult('function handle',max(abs(Vext.DV-DV(disc.x(:)))) < 1e-10);

    DV                  =   DV(disc.x(:));
    Vext.set('externalPotentialDerivative',DV);
    Vext.initialize(DFT);   Vext.getPotentialDerivative(1);                 % Store result in Vext.DV
    obj.showResult('user-defined discretization',max(abs(Vext.DV-DV)) < 1e-10);
    
    Vext.set('externalPotential',[],'atom',Va);
    Vext.initialize(DFT);   Vext.getPotentialDerivative(1);                 % Store result in Vext.DV

    DV                  =   DV + Va{1}.getPotentialDerivative(1,disc.x(:)) + ...
                                 Va{2}.getPotentialDerivative(1,disc.x(:)) + ...
                                 Va{3}.getPotentialDerivative(1,disc.x(:));
    R                   =   max(abs( Vext.DV-DV )) < 1e-10;
    obj.showResult('list of atomic centers' ,R);
    if ~R, fprintf('      (check QMol_Va_Gaussian and QMol_Va_softCoulomb components)\n');    end

end
function runUsePotential(obj) %=========================================
%runUsePotential unit tests using the external potential
    
    % Initialization
    x                   =   QMol_disc('x',-10:.1:20);
    DFT                 =   QMol_DFT_spinRes('discretization',x);
    x.initialize(DFT);
    
    V                   =   @(x) -3*x./sqrt((x-sqrt(2)).^2 + exp(2));
    Vext                =   QMol_DFT_Vext('externalPotential',V);
    Vext.initialize(DFT);

    % External potential operator (spin restricted) ~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('External potential operator (spin restricted)');

    Vks                 =   Vext.getPotential;
    R                   =   isempty(Vks.potentialUp)                        &&  ... spin restricted potential
                            isempty(Vks.potentialDown)                      &&  ...
                            all(size(Vks.potential) == size(x.x(:)))        &&  ... same size
                            all(Vks.potential == V(x.x(:)));                      % same values
    obj.showResult('getPotential (create new potential)',R);

    Vks                 =   Vext.getPotential(Vks,true);
    R                   =   isempty(Vks.potentialUp)                        &&  ... spin restricted potential
                            isempty(Vks.potentialDown)                      &&  ...
                            all(size(Vks.potential) == size(x.x(:)))        &&  ... same size
                            all(Vks.potential == V(x.x(:))+V(x.x(:)));            % same values
    obj.showResult('getPotential (add to existing potential)',R);

    Vks                 =   Vext.getPotential(Vks);
    R                   =   isempty(Vks.potentialUp)                        &&  ... spin restricted potential
                            isempty(Vks.potentialDown)                      &&  ...
                            all(size(Vks.potential) == size(x.x(:)))        &&  ... same size
                            all(Vks.potential == V(x.x(:)));                      % same values
    obj.showResult('getPotential (overwrite existing potential)',R);


    % External potential operator (spin polarized) ~~~~~~~~~~~~~~~~~~~~~~~~
    DFT                 =   QMol_DFT_spinPol('discretization',x);
    x.initialize(DFT);      Vext.reset;     Vext.initialize(DFT);

    obj.showSection('External potential operator (spin polarized)');

    Vks                 =   Vext.getPotential;
    R                   =   isempty(Vks.potential)                          &&  ... spin polarized potential
                            all(size(Vks.potentialUp) == size(x.x(:)))      &&  ... same size
                            all(size(Vks.potentialDown) == size(x.x(:)))    &&  ... 
                            all(Vks.potentialUp == V(x.x(:)))               &&  ... same values
                            all(Vks.potentialDown == V(x.x(:)));
    obj.showResult('getPotential (create new potential)',R);

    Vks                 =   Vext.getPotential(Vks,true);
    R                   =   isempty(Vks.potential)                          &&  ... spin polarized potential
                            all(size(Vks.potentialUp) == size(x.x(:)))      &&  ... same size
                            all(size(Vks.potentialDown) == size(x.x(:)))    &&  ... 
                            all(Vks.potentialUp == V(x.x(:))+V(x.x(:)))     &&  ... same values
                            all(Vks.potentialDown == V(x.x(:))+V(x.x(:)));
    obj.showResult('getPotential (add to existing potential)',R);

    Vks                 =   Vext.getPotential(Vks);
    R                   =   isempty(Vks.potential)                          &&  ... spin polarized potential
                            all(size(Vks.potentialUp) == size(x.x(:)))      &&  ... same size
                            all(size(Vks.potentialDown) == size(x.x(:)))    &&  ... 
                            all(Vks.potentialUp == V(x.x(:)))               &&  ... same values
                            all(Vks.potentialDown == V(x.x(:)));
    obj.showResult('getPotential (overwrite existing potential)',R);

    % External potential derivative (spin restricted) ~~~~~~~~~~~~~~~~~~~~~
    DFT                 =   QMol_DFT_spinRes('discretization',x);
    x.initialize(DFT);      Vext.reset;     Vext.initialize(DFT);

    obj.showSection('External potential derivative operator (spin restricted)');

    DV                  =   @(x) -3./sqrt((x-sqrt(2)).^2 + exp(2)) + 3*x.*(x-sqrt(2)).*((x-sqrt(2)).^2 + exp(2)).^-1.5;
    Vext.set('externalPotentialDerivative',DV);
    Vext.initialize(DFT);

    Vks                 =   Vext.getPotentialDerivative(1);
    R                   =   isempty(Vks.potentialGradientUp)                    &&  ... spin restricted potential
                            isempty(Vks.potentialGradientDown)                  &&  ...
                            all(size(Vks.potentialGradient) == size(x.x(:)))    &&  ... same size
                            all(Vks.potentialGradient == DV(x.x(:)));                 % same values
    obj.showResult('getPotentialDerivative (create new potential)',R);

    Vks                 =   Vext.getPotentialDerivative(1,Vks,true);
    R                   =   isempty(Vks.potentialGradientUp)                    &&  ... spin restricted potential
                            isempty(Vks.potentialGradientDown)                  &&  ...
                            all(size(Vks.potentialGradient) == size(x.x(:)))    &&  ... same size
                            all(Vks.potentialGradient == DV(x.x(:))+DV(x.x(:)));      % same values
    obj.showResult('getPotentialDerivative (add to existing potential)',R);

    Vks                 =   Vext.getPotentialDerivative(1,Vks);
    R                   =   isempty(Vks.potentialGradientUp)                    &&  ... spin restricted potential
                            isempty(Vks.potentialGradientDown)                  &&  ...
                            all(size(Vks.potentialGradient) == size(x.x(:)))    &&  ... same size
                            all(Vks.potentialGradient == DV(x.x(:)));                 % same values
    obj.showResult('getPotentialDerivative (overwrite existing potential)',R);

    % External potential derivative (spin polarized) ~~~~~~~~~~~~~~~~~~~~~
    DFT                 =   QMol_DFT_spinPol('discretization',x);
    x.initialize(DFT);      Vext.reset;     Vext.initialize(DFT);

    obj.showSection('External potential derivative operator (spin polarized)');

    Vks                 =   Vext.getPotentialDerivative(1);
    R                   =   isempty(Vks.potentialGradient)                      &&  ... spin polarized potential
                            all(size(Vks.potentialGradientUp) == size(x.x(:)))  &&  ... same size
                            all(size(Vks.potentialGradientDown) == size(x.x(:)))&&  ... 
                            all(Vks.potentialGradientUp == DV(x.x(:)))          &&  ...
                            all(Vks.potentialGradientDown == DV(x.x(:)));
    obj.showResult('getPotentialDerivative (create new potential)',R);

    Vks                 =   Vext.getPotentialDerivative(1,Vks,true);
    R                   =   isempty(Vks.potentialGradient)                      &&  ... spin polarized potential
                            all(size(Vks.potentialGradientUp) == size(x.x(:)))  &&  ... same size
                            all(size(Vks.potentialGradientDown) == size(x.x(:)))&&  ... 
                            all(Vks.potentialGradientUp == DV(x.x(:))+DV(x.x(:)))&& ... same values
                            all(Vks.potentialGradientDown == DV(x.x(:))+DV(x.x(:)));
    obj.showResult('getPotentialDerivative (add to existing potential)',R);

    Vks                 =   Vext.getPotentialDerivative(1,Vks);
    R                   =   isempty(Vks.potentialGradient)                      &&  ... spin polarized potential
                            all(size(Vks.potentialGradientUp) == size(x.x(:)))  &&  ... same size
                            all(size(Vks.potentialGradientDown) == size(x.x(:)))&&  ... 
                            all(Vks.potentialGradientUp == DV(x.x(:)))          &&  ...
                            all(Vks.potentialGradientDown == DV(x.x(:)));
    obj.showResult('getPotentialDerivative (overwrite existing potential)',R);

    % External energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DFT                 =   QMol_DFT_spinRes('discretization',x);
    x.initialize(DFT);      Vext.reset;     Vext.initialize(DFT);

    obj.showSection('External energy');
    x0                  =   5;
    s                   =   2;
    rho                 =   QMol_DFT_density('rho',exp(-(x.x(:)-x0).^2 * .5/s^2));
    E                   =   sum(V(x.x(:)).*rho.rho)*(x.x(2)-x.x(1));
    obj.showResult('getEnergy (spin restricted)',abs(E-Vext.getEnergy(rho)) < 1e-10);
    
    DFT                 =   QMol_DFT_spinPol('discretization',x);
    x.initialize(DFT);      Vext.reset;     Vext.initialize(DFT);
    
    rho.set('rho',[],'rhoUp',x.x(:).*exp(-(x.x(:)-x0).^2 * .5/s^2),...
                     'rhoDw',x.x(:).^2.*exp(-(x.x(:)-x0).^2 * .5/s^2))
    E                   =   [sum(V(x.x(:)).*rho.rhoUp),sum(V(x.x(:)).*rho.rhoDw)]*(x.x(2)-x.x(1));
    obj.showResult('getEnergy (spin polarized)',all(abs(E-Vext.getEnergy(rho)) < 1e-10));

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

