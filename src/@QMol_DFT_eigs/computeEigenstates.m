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
    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility

    opt                 =  {'IsFunctionSymmetric',  obj.issym,  ...
                            'Tolerance',            obj.tol,    ...
                            'MaxIterations',        obj.maxit,  ...
                            'SubspaceDimension',    min(obj.p,obj.DFT.disc.nV-1),...  % Cap Krylov space to "basis" size
                            'Display',              obj.disp    };

    % Hamiltonian operator
    disc                =   obj.DFT.disc;                                   % Copy handle to discretization object
    dim                 =   disc.DFT_sizeOrbital;                           % Size of the Kohn-Sham orbitals
    isRes               =   ~( numel(dim) == 2   &&   dim(2) == 1 );        % Do we need to resize?

    function Hp = Hop(p,S,isUp)
        % Reshape wave function
        if isRes,   p   =   reshape(p,dim);     end

        % Apply Hamiltonian operator                                          Eigen states are real valued
        if obj.DFT.isSpinPol
            Hp          =   real(disc.DFT_operatorHamiltonian(V,p,isUp,S));
        else
            Hp          =   real(disc.DFT_operatorHamiltonian(V,p,     S));
        end

        % Reshape back the wave function
        if isRes,   Hp  =   Hp(:);               end
    end

    % Compute eigenstates
    KSO                 =   obj.DFT.KSO;

    if obj.DFT.isSpinPol
        % SPIN UP COMPONENT
        N               =   0;
        KSO.KSOup       =   NaN(prod(dim),sum(obj.symState{1}(1,:)));       % Orbitals expanded as column vectors
        Eup             =   NaN(sum(obj.symState{1}(1,:)));                 % Will take the diagonal at the end
        
        % Parse symmetry groups
        for k = 1:size(obj.symState{1},2)
            % Miscellaneous
            ind         =   N + (1:obj.symState{1}(1,k));                   % Index of where to store results
            N           =   N + obj.symState{1}(1,k);                       % Update last filled index counter
            p0          =   disc.DFT_randomOrbital( ...                     % Random initial guess with proper symmetry
                                randStr,obj.symState{1}(2:end,k));

            % Compute eigenstates and energies
            [KSO.KSOup(:,ind),Eup(ind,ind)] ...
                        =   eigs(@(p) Hop(p,obj.symState{1}(2:end,k),true), ... Hamiltonian operator with proper symmetry
                                prod(dim),obj.symState{1}(1,k),'SA',      ... Parameters for eigs function
                                'StartVector',p0,opt{:});                   % Initial condition and other options
        end
        
        % Sort orbitals by increasing energy values
        if size(obj.symState{1},2) > 1
            [Eup,ind]   =   sort(diag(Eup));
            KSO.KSOup   =   KSO.KSOup(:,ind);
        else
            Eup         =   diag(Eup);
        end

        % SPIN DOWN COMPONENT
        N               =   0;
        KSO.KSOdw       =   NaN(prod(dim),sum(obj.symState{2}(1,:)));       % Orbitals expanded as column vectors
        Edw             =   NaN(sum(obj.symState{2}(1,:)));                 % Will take the diagonal at the end
        
        % Parse symmetry groups
        for k = 1:size(obj.symState{2},2)
            % Miscellaneous
            ind         =   N + (1:obj.symState{2}(1,k));                   % Index of where to store results
            N           =   N + obj.symState{2}(1,k);                       % Update last filled index counter
            p0          =   disc.DFT_randomOrbital( ...                     % Random initial guess with proper symmetry
                                randStr,obj.symState{2}(2:end,k));

            % Compute eigenstates and energies
            [KSO.KSOdw(:,ind),Edw(ind,ind)] ...
                        =   eigs(@(p) Hop(p,obj.symState{2}(2:end,k),false), ... Hamiltonian operator with proper symmetry
                                prod(dim),obj.symState{2}(1,k),'SA',      ... Parameters for eigs function
                                'StartVector',p0,opt{:});                   % Initial condition and other options
        end
        
        % Sort orbitals by increasing energy values
        if size(obj.symState{2},2) > 1
            [Edw,ind]   =   sort(diag(Edw));
            KSO.KSOdw   =   KSO.KSOdw(:,ind);
        else
            Edw         =   diag(Edw);
        end



        % Reshape and normalize orbitals
        if isRes
            KSO.KSOup   =   reshape(KSO.KSOup,[dim, sum(obj.symState{1}(1,:))]);
            KSO.KSOdw   =   reshape(KSO.KSOdw,[dim, sum(obj.symState{2}(1,:))]);
        end
        KSO.KSOup       =   disc.DFT_normalizeOrbital(KSO.KSOup);
        KSO.KSOdw       =   disc.DFT_normalizeOrbital(KSO.KSOdw);

        if sum(obj.symState{1}(1,:)) ~= numel(obj.DFT.occ{1})   ||          ...
           sum(obj.symState{2}(1,:)) ~= numel(obj.DFT.occ{2})
            % Mismatch in number of eigenstates computed and orbitals         This is VERY discouraged practice!
            disc.DFT_allocateOrbital([numel(obj.DFT.occ{1}) numel(obj.DFT.occ{2})],KSO); % KSO allocation can handle size mismatch

            if sum(obj.symState{1}(1,:)) > numel(obj.DFT.occ{1})
                Eup     =   Eup(1:numel(obj.DFT.occ{1}));                   % Truncate energies
            else
                Eup     =   [Eup; zeros(numel(obj.DFT.occ{1})-numel(Eup))]; % Pad energies with zeros
            end

            if sum(obj.symState{2}(1,:)) > numel(obj.DFT.occ{2})
                Edw     =   Edw(1:numel(obj.DFT.occ{2}));                   % Truncate energies
            else
                Edw     =   [Edw; zeros(numel(obj.DFT.occ{2})-numel(Edw))]; % Pad energies with zeros
            end
        end

        % Output energies are combined in a cell
        E               =   {Eup,Edw};

    else % SPIN RESTRICTED
        KSO             =   obj.DFT.KSO;                                    % Copy handle to KSO object

        N               =   0;
        KSO.KSO         =   NaN(prod(dim),sum(obj.symState(1,:)));          % Orbitals expanded as column vectors
        E               =   NaN(sum(obj.symState(1,:)));                    % Will take the diagonal at the end

        % Parse symmetry groups
        for k = 1:size(obj.symState,2)
            % Miscellaneous
            ind         =   N + (1:obj.symState(1,k));                      % Index of where to store results
            N           =   N + obj.symState(1,k);                          % Update last filled index counter
            p0          =   disc.DFT_randomOrbital( ...                     % Random initial guess with proper symmetry
                                randStr,obj.symState(2:end,k));

            % Compute eigenstates and energies
            [KSO.KSO(:,ind),E(ind,ind)] ...
                        =   eigs(@(p) Hop(p,obj.symState(2:end,k)),       ... Hamiltonian operator with proper symmetry
                                prod(dim),obj.symState(1,k),'SA',         ... Parameters for eigs function
                                'StartVector',p0,opt{:});                   % Initial condition and other options
        end

        % Sort orbitals by increasing energy values
        if size(obj.symState,2) > 1
            [E,ind]     =   sort(diag(E));
            KSO.KSO     =   KSO.KSO(:,ind);
        else
            E           =   diag(E);
        end

        % Reshape and normalize orbitals
        if isRes, KSO.KSO = reshape(KSO.KSO,[dim, sum(obj.symState(1,:))]); end
        KSO.KSO         =   disc.DFT_normalizeOrbital(KSO.KSO);

        if sum(obj.symState(1,:)) ~= numel(obj.DFT.occ)
            % Mismatch in number of eigenstates computed and orbitals         This is VERY discouraged practice!
            disc.DFT_allocateOrbital(numel(obj.DFT.occ),KSO);               % KSO allocation can handle size mismatch

            if sum(obj.symState(1,:)) > numel(obj.DFT.occ)
                E       =   E(1:numel(obj.DFT.occ));                        % Truncate energies
            else
                E       =   [E; zeros(numel(obj.DFT.occ)-numel(E))];        % Pad energies with zeros
            end
        end
        
    end

end