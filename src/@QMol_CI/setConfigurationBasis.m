function setConfigurationBasis(obj)
%setConfigurationBasis sets the configuration-state basis
%   Use setConfigurationBasis to set the configurationBasis property
%   associated with the type and active-space specifications in the CI 
%   object. This is typically done between the creation of the CI object
%   and the calculation of the CI matrix, using the computeCImatrix method.
%
%   CI.setConfigurationBasis
%
%   Supported type of configuration-state builders are
%     CIS       single excitations
%     CISD      single and double excitations
%     RAS       restricted active space
%
%     In multi-reference CIS(D) models, configurationBasis first contains
%     the refernce state(s), followed by the single excitations from those
%     references, followed by the double excitations, if any. Duplicate
%     states are then removed with cleanConfigurationBasis.
%
%   Options for building the configuration-state basis are
%     reference > Idexes for the spin orbitals making up the reference
%                 state. For instance, [-1 -2 1 2] corresponds to
%                 the reference state with two electron in the first two
%                 spatial orbitals. For CIS(D) calculations, the single
%                 (and double) excitations are defined with respect to that
%                 reference. For RAS calculations, the reference is used to
%                 determine the number of electrons in each of the up and
%                 down spin channels.
%               > Multiple references can be specified in successive rows.
%                 For instance, [-1 -2 1 2; -1 -3 1 3] corresponds to the
%                 first reference state with two electron in the first two
%                 spatial orbitals and the second reference with two
%                 electrons in the first and third spatial orbitals. For
%                 CIS(D) calculations, the configuration state basis will
%                 include all the unique single (and double) excitations
%                 from each of the references.
%               > When a reference is specified, numberElectron is updated
%                 to reflect the number of electrons associated with the
%                 configuration basis.
%               > If no reference is specified, the configuration assumes a
%                 single reference close shell molecule, with as many
%                 electrons in the up and down spin channels. The first 
%                 numberElectron/2 spatial orbitals are used for the 
%                 reference. For instance, if numberElectron is equal to 4,
%                 the reference is set to [-1 -2 1 2].
%     active    > Indexes for the active spin orbitals, which may contain 0
%                 or 1 electrons. If left empty, all the orbitals in the
%                 orbitalBasis are used for the active space. For CIS(D),
%                 spin orbitals of the reference state are implicitly
%                 included in the active-orbital list. For RAS, one must
%                 specify all the spin orbitals to include in the active
%                 space, but frozen ones can be omitted.
%               > Multi-reference CIS(D) calculations us the same active
%                 space for all references.
%     frozen    > Indexes of the frozen spin orbitals, which must contain 1
%                 electron. For CIS(D) calculations, only frozen orbitals
%                 that are in the reference are considered. For RAS, frozen
%                 orbitals are implicitly included in the active space. If
%                 left empty, no orbitals are frozen.
%     noDouble  > Indexes of the spatial orbitals that may not contain 2
%                 electrons, i.e., have an electron in both the up and down
%                 spin channels. If left empty, no restriction on double
%                 occupation is imposed. For CIS(D) calculations, the
%                 constraint is not imposed on the reference configuration.
%     noEmpty   > Indexes of the spatial orbitals that may not be left
%                 empty, i.e., must have at least 1 electron in the up or
%                 down spin channel. If left empty, no restriction on empty
%                 occupation is imposed. For CIS(D) calculations, the
%                 constraint is not imposed on the reference configuration.
%
%   Note: setConfigurationBasis resets the QMol_CI_conv object and clears
%   the CImatrix and DXmatrix properties before setting the configuration 
%   state basis.

% Algorithm mapping:
%     algo == 0 > RAS
%     algo == 1 > CIS
%     algo == 2 > CISD

    % Initialization ======================================================
        % (Re)initialize everythings
        obj.clear('CI','DX');   if obj.disp,    QMol_doc.showHeader;        end     % clear resets the object
        
        
        if obj.disp,    QMol_doc.showSection('Build the configuration-interaction (CI) model [Visentin 2025]'); %#ok<ALIGN> 
                        obj.initialize(2);
                        QMol_doc.showBibliography({'Visentin 2025'});
                        QMol_doc.showFunding;
                        QMol_doc.showFooter;
    
                        QMol_doc.showSection('Set the configuration-state basis');
                        fprintf('  * Initialization                                                   ');

        else,           obj.initialize(0);                                  end

        % Configuration basis parameters
        [algo,locN,locRef,locAct,locFrz] = obj.getConfigurationBasisParameters(true);
        obj.N           =   locN;

        if obj.disp,    fprintf(' DONE\n');                                 end

        function [iExP,nExP,iExN,nExN,iiHo,nHoP,nHoN] = getCISDparam(k)
            % Allowed spin-orbital excitations
            iExP        =   locAct{k}(locAct{k} > 0);       nExP    =   numel(iExP);
            iExN        =   locAct{k}(locAct{k} < 0);       nExN    =   numel(iExN);
            % Allowed spin-orbital holes
            [~,iiHo]    =   setdiff(locRef(k,:),locFrz);
            nHoP        =   sum(locRef(k,iiHo) > 0);        nHoN    =   sum(locRef(k,iiHo) < 0);
        end

        % CIS(D) parameters
        if algo > 0
            % Number of basis vectors
            nbRef       =   size(locRef,1);

            nCSB        =   nbRef;                                          % number of reference states
            for kk = 1:nbRef
                [~,nExP,~,nExN,~,nHoP,nHoN] = getCISDparam(kk);
                nCSB    =   nCSB + nHoP*nExP + nHoN*nExN;                   % number of single excitations for the reference
                if algo > 1                                                 % double excitations
                    nCSB=   nCSB + nHoP*(nHoP-1)/2 * nExP*(nExP-1)/2;       %  > both in up spin channel
                    nCSB=   nCSB + nHoN*(nHoN-1)/2 * nExN*(nExN-1)/2;       %  > both in down spin channel
                    nCSB=   nCSB + nHoP*nExP * nHoN*nExN;                   %  > in different spin channels
                end
            end
        % RAS parameters
        else
            % Active orbitals
            locAct      =   locAct{1};              nbRef   =   1;          % single reference/effective active space for RAS
            iActP       =   locAct(locAct > 0);     nActP   =   sum(locRef > 0) - sum(locFrz > 0);
            iActN       =   locAct(locAct < 0);     nActN   =   sum(locRef < 0) - sum(locFrz < 0);
            % Check compatibility of results
            if nActP < 0 || nActN < 0,  error('QMol:CI:RAS','Incompatible reference and frozen spin orbitals');     end
            if numel(iActP) < nActP,    error('QMol:CI:RAS','Incompatible reference and active space (spin up)');   end
            if numel(iActN) < nActN,    error('QMol:CI:RAS','Incompatible reference and active space (spin down)'); end
        end

    % Set configuration basis =============================================
    function n = setReference(n) %~~~~~~~~~~
        % Copy references
        obj.CSB(n:n+nbRef-1,:)  =   locRef;

        % Update counter
        n                       =   n + nbRef;
    end
    function n = setSingles(kk,n) %~~~~~~~~~
        % Initialization
        N                       =   nHoP*nExP + nHoN*nExN;
        obj.CSB(n:n+N-1,:)      =   repmat(locRef(kk,:),[N,1]);
                                                                                    for l = iiHo.'
        % CIS determinants                                     
        if locRef(kk,l) > 0, obj.CSB(n:n+nExP-1,l) = iExP;  n = n+nExP;                 % up spin excitation
        else,                obj.CSB(n:n+nExN-1,l) = iExN;  n = n+nExN;     end,    end % down spin excitation

    end
    function n = setDoubles(kk,n) %~~~~~~~~~
        % Initialization
        N                       =   nHoP*(nHoP-1)/2 * nExP*(nExP-1)/2 + ...    > both in up spin channel
                                    nHoN*(nHoN-1)/2 * nExN*(nExN-1)/2 + ...    > both in down spin channel
                                    nHoP*nExP * nHoN*nExN;                  %  > in different spin channels
        obj.CSB(n:n+N-1,:)      =   repmat(locRef(kk,:),[N,1]);

        for ik = 1:numel(iiHo), k = iiHo(ik);       for il = ik+1:numel(iiHo),  l = iiHo(il); %#ok<ALIGN>
            if locRef(kk,k) > 0     % 1st excitation in the up-spin channel
                if locRef(kk,l) > 0     % 2nd excitation in the same (up-spin) channel
                    for m = 1:nExP,     obj.CSB(n:n+nExP-m-1,k) =   iExP(m);            %#ok<ALIGN>
                                        obj.CSB(n:n+nExP-m-1,l) =   iExP(m+1:end);
                                        n       =   n + nExP-m;                     end
                else                    % 2nd excitation in the opposite (down-spin channel)
                    for m = 1:nExP,     obj.CSB(n:n+nExN-1,k)   =   iExP(m);            %#ok<ALIGN>
                                        obj.CSB(n:n+nExN-1,l)   =   iExN;
                                        n       =   n + nExN;                       end
                end
            else                    % 1st excitation in the down-spin channel
                if locRef(kk,l) > 0     % 2nd excitation in the opposite (up-spin) channel
                    for m = 1:nExN,     obj.CSB(n:n+nExP-1,k)   =   iExN(m);            %#ok<ALIGN>
                                        obj.CSB(n:n+nExP-1,l)   =   iExP;
                                        n       =   n + nExP;                       end
                else                    % 2nd excitation in the same (down-spin channel)
                    for m = 1:nExN,     obj.CSB(n:n+nExN-m-1,k) =   iExN(m);            %#ok<ALIGN>
                                        obj.CSB(n:n+nExN-m-1,l) =   iExN(m+1:end);
                                        n       =   n + nExN-m;                     end
                end
            end
        end, end
    end

    % CIS(D) type calculation
    if algo > 0
        % Initialization
        obj.CSB         =   NaN(nCSB,locN);

        % References
        n               =   setReference(1);

        % Single excitation
        for kk = 1:nbRef,   [iExP,nExP,iExN,nExN,iiHo,nHoP,nHoN] = getCISDparam(kk);
            n           =   setSingles(kk,n);
        end
                                                                            if algo > 1
        % Double excitations
        for kk = 1:nbRef,   [iExP,nExP,iExN,nExN,iiHo,nHoP,nHoN] = getCISDparam(kk);
            n           =   setDoubles(kk,n);
        end,                                                                end
    % RAS type calculation
    else
        % Initialization
        iRasP           =   nchoosek(iActP,nActP);      nRasP   =   size(iRasP,1);
        iRasN           =   nchoosek(iActN,nActN);      nRasN   =   size(iRasN,1);
        nCSB            =   nRasP*nRasN;

        obj.CSB         =   NaN(nCSB,locN);

        n               =   1;

        % Frozen spin down
        nFr                     =   sum(locFrz < 0);                        if nFr > 0
        fr                      =   locFrz(locFrz < 0);
        obj.CSB(:,n:n+nFr-1)    =   repmat(fr(:).',[nCSB,1]);
        n                       =   n + nFr;                                end

        % Down spin channel
        obj.CSB(:,n:n+nActN-1)  =   repmat(iRasN,[nRasP,1]);
        n                       =   n + nActN;

        % Frozen spin up
        nFr                     =   sum(locFrz > 0);                        if nFr > 0
        fr                      =   locFrz(locFrz > 0);
        obj.CSB(:,n:n+nFr-1)    =   repmat(fr(:).',[nCSB,1]);
        n                       =   n + nFr;                                end

        % Up spin channel
        for kk = 1:nRasN,    obj.CSB(kk:nRasN:end,n:n+nActP-1)    =   iRasP;  end

    end

    % Eliminate forbidden configurations ==================================
        % Initialization
        ind             =   true(nCSB,1);

        % Forbiddent double excitation
        for kk = 1:numel(obj.noDbl)
            ind         =   ind & ~(any(obj.CSB == obj.noDbl(kk),2) & ...
                                    any(obj.CSB ==-obj.noDbl(kk),2) );
        end

        % Forbidden emmpty orbitals
        for kk = 1:numel(obj.noEmpt)
            ind         =   ind & ~(all(obj.CSB ~= obj.noEmpt(kk),2) & ...
                                    all(obj.CSB ~=-obj.noEmpt(kk),2) );
        end

        % CIS(D) case
        if algo > 0,        ind(1:nbRef)    =   true;                       end % exclude reference from constraints

        % Remove forbidden configurations
        obj.CSB         =   obj.CSB(ind,:);

        % For MR-CIS(D) remove duplicate configurations
        if nbRef > 1,       obj.cleanConfigurationBasis;                    end

    % Finalization ========================================================
        % Update configuration-state basis tracker
        obj.isSetConfigBasis=   true;

        % Show configuration-state 
        if obj.disp
            obj.showConfigBasisSet;     obj.version;    fprintf('\n');
            QMol_doc.showFooter;                        fprintf('\n');
        end
end
%% Internal functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
