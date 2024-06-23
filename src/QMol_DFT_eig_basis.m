classdef QMol_DFT_eig_basis < QMol_suite
%QMol_DFT_eig_basis eigen-state solver for density-functional-theory (DFT)
%   level of theory discretized on a basis.
    
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT_eig_basis})
function showInfo
    fprintf('  * QMol_DFT_eig_basis:\n'); 
    fprintf('      > DFT eigen-state solver for basis-set models\n      > MATLAB eig function\n'); 
    QMol_DFT_eig_basis.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % General options
    fprintf('  * Eigen-state solver for DFT Hamiltonians            MATLAB eig function\n')
    fprintf('    using a direct diagonalization of the DFT Hamiltonian matrix.\n');
    
    % Version
    obj.version;
    ref                 =   {};
end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Linked objects
    DFT
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object
    
    % Run-time variables
    obj.DFT             =   [];
    
    % Initialization status
    obj.isInit          =   false;
end
end
methods (Access=public)
function initialize(obj,DFT)
%initialize initializes the eigen solver
    
    % Initialization
    obj.reset;

    % Set links
    obj.DFT             =   DFT;

    % House keeping
    obj.isInit       =   true;
end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_DFT_eig_basis';
    PropNames           =   {};
end
end
%% Eigenstate computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function E = computeEigenstates(obj,V)
%computeEigenstates computes the eigenstates, and corresponding energies,
%   associated with the member DFT and input potential objects. The result
%   eigenstates are stored in the member DFT's Kohn-Sham orbitals and only
%   the energies are returned. Omitting or passing an empty potential
%   object triggers the computation of the Kohn-Sham potential associated
%   with the current state of the DFT object.
    
    % Initialization
    if nargin < 2,  V   =   [];         end
    if isempty(V)
        obj.DFT.Vks     =   obj.DFT.getPotential([],obj.DFT.Vks);           % Use DFT-object memory slot
        V               =   obj.DFT.Vks;
    end
    KSO                 =   obj.DFT.KSO;
    
    % Compute eigen states
    if obj.DFT.isSpinPol
        % Spin up component
        if isempty(V.Vimp)
            mH          =   obj.DFT.disc.mT + V.mVup;
        else
            v           =   eye(obj.DFT.disc.nV);
            mH          =   obj.DFT.disc.mT;        for k = 1:obj.DFT.disc.nV
            mH(:,k)     =   mH(:,k) + V.applyPotential(v(:,k),true);        end
        end

        [VV,Eup]        =   eig(mH);

        [~,ind]         =   sort(real(diag(Eup)));    
        ind             =   ind(1:numel(obj.DFT.occ{1}));
        Eup             =   real(diag(Eup(ind,ind)));
        KSO.KSOup       =   obj.DFT.disc.DFT_normalizeOrbital(VV(:,ind));

        % Spin down component
        if isempty(V.Vimp)
            mH          =   obj.DFT.disc.mT + V.mVdw;
        else
            v           =   eye(obj.DFT.disc.nV);
            mH          =   obj.DFT.disc.mT;        for k = 1:obj.DFT.disc.nV
            mH(:,k)     =   mH(:,k) + V.applyPotential(v(:,k),false);       end
        end

        [VV,Edw]        =   eig(mH);

        [~,ind]         =   sort(real(diag(Edw)));    
        ind             =   ind(1:numel(obj.DFT.occ{2}));
        Edw             =   real(diag(Edw(ind,ind)));
        KSO.KSOdw       =   obj.DFT.disc.DFT_normalizeOrbital(VV(:,ind));


        % Output energies are combined in a cell
        E               =   {Eup,Edw};
    else
        % Spin restricted
        if isempty(V.Vimp)
            mH          =   obj.DFT.disc.mT + V.mV;
        else
            v           =   eye(obj.DFT.disc.nV);
            mH          =   obj.DFT.disc.mT;        for k = 1:obj.DFT.disc.nV
            mH(:,k)     =   mH(:,k) + V.applyPotential(v(:,k));             end
        end

        [VV,E]          =   eig(mH);

        [~,ind]         =   sort(real(diag(E)));
        ind             =   ind(1:numel(obj.DFT.occ));
        E               =   real(diag(E(ind,ind)));
        KSO.KSO         =   obj.DFT.disc.DFT_normalizeOrbital(VV(:,ind));

    end

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

