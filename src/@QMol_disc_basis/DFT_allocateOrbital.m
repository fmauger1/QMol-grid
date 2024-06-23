function KSO = DFT_allocateOrbital(obj,N,KSO,randStr,isR)
%DFT_allocateOrbital allocates a Kohn-Sham-orbital (KSO) object with the 
%   specified number of orbitals. Optionally, if a KSO object in given as
%   input, it is (re)initialized instead.
    
    % Initialization
    if nargin < 5,  isR =   false;
    if nargin == 2, KSO =   [];             end, end
    if iscell(N),   N   =   [N{1} N{2}];    end

    % Allocate KSO object (to zero/rand for compatibility with resizing)
    if obj.QM.isSpinPol     % Spin polarized
        % Allocate or resize orbital(s) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if nargin < 4, if isempty(KSO)                                      %#ok<ALIGN> 
            % Create (spin-polarized) object from scratch
            KSO         =   QMol_DFT_orbital_basis('isSpinPol',true,      ...
                                'KSOup',zeros(size(obj.v,2),N(1)),   ...
                                'KSOdw',zeros(size(obj.v,2),N(2)));
        elseif ~KSO.isSpinPol || any(size(obj.v,2) ~= [size(KSO.KSOup,1) size(KSO.KSOdw,1)])
            % Different spin polarization or incompatible sizes
            KSO.set('isSpinPol',true,'KSO',[],                      ...
                    'KSOup',zeros(size(obj.v,2),N(1)),'KSOdw',zeros(size(obj.v,2),N(2)));
        else
            % Initialization
            SU          =   size(KSO.KSOup);
            SD          =   size(KSO.KSOdw);

            % Resize (if needed)
            if SU(2) > N(1),    KSO.KSOup   =   KSO.KSOup(:,1:N(1));
            else,               KSO.KSOup   =  [KSO.KSOup zeros(size(obj.v,2),N(1)-SU(2))];     end

            if SD(2) > N(2),    KSO.KSOdw   =   KSO.KSOdw(:,1:N(2));
            else,               KSO.KSOdw   =  [KSO.KSOdw zeros(size(obj.v,2),N(2)-SD(2))];   end, end

        % Allocate random orbital ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else, if isempty(KSO),  KSO         =   QMol_DFT_orbital_basis();                               end
            if isR,             KSO.set('isSpinPol',true,'KSO',[],          ...
                                        'KSOup',rand(randStr,[size(obj.v,2),N(1)])-.5,                  ...
                                        'KSOdw',rand(randStr,[size(obj.v,2),N(2)])-.5);                 %#ok<ALIGN> 
            else,               KSO.set('isSpinPol',true,'KSO',[],          ...
                                        'KSOup',rand(randStr,[size(obj.v,2),N(1)],'like',1i)-.5-.5i,    ...
                                        'KSOdw',rand(randStr,[size(obj.v,2),N(2)],'like',1i)-.5-.5i);   end
            % Normalize orbitals
            KSO.KSOup   =   obj.DFT_normalizeOrbital(KSO.KSOup);
            KSO.KSOdw   =   obj.DFT_normalizeOrbital(KSO.KSOdw);
        end
    else                    % Spin restricted
        % Allocate or resize orbital(s) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if nargin < 4, if isempty(KSO)                                      %#ok<ALIGN> 
            % Create (spin-polarized) object from scratch
            KSO         =   QMol_DFT_orbital_basis('isSpinPol',false,'KSO',zeros(size(obj.v,2),N));

        elseif KSO.isSpinPol || size(obj.v,2) ~= size(KSO.KSO,1)
            % Different spin polarization or incompatible sizes
            KSO.set('isSpinPol',false,'KSO',zeros(size(obj.v,2),N),'KSOup',[],'KSOdw',[]);
        else
            % Initialization
            S           =   size(KSO.KSO);

            % Resize (if needed)
            if S(2) > N(1),     KSO.KSO     =   KSO.KSO(:,1:N);
            else,               KSO.KSO     =  [KSO.KSO zeros(size(obj.v,2),N(1)-S(2))]; end, end

        % Allocate random orbital ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else, if isempty(KSO),  KSO         =   QMol_DFT_orbital_basis();                           end
            if isR,             KSO.set('isSpinPol',false,'KSOup',[],'KSOdw',[],                    ...
                                        'KSO',rand(randStr,[size(obj.v,2),N])-.5);                  %#ok<ALIGN> 
            else,               KSO.set('isSpinPol',false,'KSOup',[],'KSOdw',[],                    ...
                                        'KSO',rand(randStr,[size(obj.v,2),N],'like',1i)-.5-.5i);    end
            % Normalize orbitals
            KSO.KSO     =   obj.DFT_normalizeOrbital(KSO.KSO);
        end
    end

    % Initialize orbital
    KSO.initialize(obj);

end