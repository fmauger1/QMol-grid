function computeEigenstates(obj)
%computeEigenstates computes the eigenstates, and corresponding energies,
%   associated with the member SE object. The result eigenstates are stored
%   in the member SE's wave function.
    
    % Initialization
    V                   =   obj.SE.V;                                       % Potential object
    Nb                  =   obj.SE.N;
    randStr             =   RandStream('dsfmt19937','Seed',0);              % For reproducibility

    opt                 =  {'IsFunctionSymmetric',  obj.issym,  ...
                            'Tolerance',            obj.tol,    ...
                            'MaxIterations',        obj.maxit,  ...
                            'SubspaceDimension',    min(obj.p,obj.SE.disc.nV-1),...  % Cap Krylov space to "basis" size
                            'Display',              obj.dispIter};

    % Hamiltonian operator
    disc                =   obj.SE .disc;                                   % Copy handle to discretization object
    dim                 =   disc.SE_sizeWaveFunction;                       % Size of the wave functions
    isRes               =   ~( numel(dim) == 2   &&   dim(2) == 1 );        % Do we need to resize?

    function Hp = Hop(p,S)
        % Reshape wave function
        if isRes,   p   =   reshape(p,dim);     end

        % Apply Hamiltonian operator
        Hp              =   real(disc.SE_operatorHamiltonian(V,p,S));

        % Reshape back the wave function
        if isRes,   Hp  =   Hp(:);               end
    end

    % Compute eigenstates
    wfcn                =   obj.SE.wfcn;

    N                   =   0;
    wfcn.wfcn           =   NaN(prod(dim),sum(obj.symState(1,:)));          % Orbitals expanded as column vectors
    E                   =   NaN(sum(obj.symState(1,:)));                    % Will take the diagonal at the end

        % Parse symmetry groups
        for k = 1:size(obj.symState,2)
            % Miscellaneous
            ind         =   N + (1:obj.symState(1,k));                      % Index of where to store results
            N           =   N + obj.symState(1,k);                          % Update last filled index counter
            p0          =   disc.SE_randomWaveFunction( ...                 % Random initial guess with proper symmetry
                                randStr,obj.symState(2:end,k));

            % Compute eigenstates and energies
            [wfcn.wfcn(:,ind),E(ind,ind)] ...
                        =   eigs(@(p) Hop(p,obj.symState(2:end,k)),       ... Hamiltonian operator with proper symmetry
                                prod(dim),obj.symState(1,k),'SA',         ... Parameters for eigs function
                                'StartVector',p0,opt{:});                   % Initial condition and other options
        end

    % Sort orbitals by increasing energy values
    if size(obj.symState,2) > 1
        [~,ind]     =   sort(diag(E));
        wfcn.wfcn   =   wfcn.wfcn(:,ind);
    end

    % Reshape and normalize orbitals
    if isRes, wfcn.wfcn = reshape(wfcn.wfcn,[dim, sum(obj.symState(1,:))]); end
    wfcn.wfcn           =   disc.SE_normalizeWaveFunction(wfcn.wfcn);

    if sum(obj.symState(1,:)) ~= Nb
        % Mismatch in number of eigenstates computed and orbitals             This is VERY discouraged practice!
        disc.SE_allocateWaveFunction(nb,wfcn);                              % KSO allocation can handle size mismatch

    end

end

