function [algo,locN,locRef,locAct,locFrz] = getConfigurationBasisParameters(obj,isBld)
%getConfigurationBasisParameters get the configuration basis parameters
%   Provides a unified (and consitent) interface for defining the effective
%   parameters for the configuration basis associated with the flavor of CI
%
%   The input parameters isBld specifies whether the configuration basis is
%   being built (for run-time documentation purposes)

    % Identify the flavor of CI
    switch lower(obj.type)
        case 'cis',     algo    =   1;
        case 'cisd',    algo    =   2;
        case 'ras',     algo    =   0;
        otherwise
            if isBld,   fprintf('*FAIL\n');                         end
            error('QMol:CI:type',['Unknown CI type ' obj.type '. Supported types are CIS.'])
    end

    % Default values
    if isempty(obj.ref),                                                            if mod(obj.N,2) ~= 0
        % Default reference assumes close-shell system
        error('QMol:CI:reference','Open-shell systems must define a reference.'),   end

        % Set default reference
        locN            =   obj.N;
        locRef          =   [-(1:0.5*obj.N),1:0.5*obj.N];
        nbRef           =   1;                                              % single reference
    else
        % Update the number of electrons
        [nbRef,locN]    =   size(obj.ref);                                  % number of references and number of electrons
        locRef          =   obj.ref;                                                % use provided reference
    end

    if isempty(obj.act)
        locAct          =   [-(1:size(obj.SOB,2)),1:size(obj.SOB,2)];       % by default all the spin orbitals are active
    else
        locAct          =   obj.act;                                        % use provided active space
    end
    locAct              =   repmat({locAct},[nbRef,1]);                     % each reference has its own effective active space

    % CIS(D) parameters
    if algo > 0,                                                for k = 1:nbRef
        % Allowed spin-orbital excitations
        locAct{k}       =   setdiff(locAct{k},locRef(k,:));     end

        % Frozen orbitals
        locFrz          =   intersect(locRef(:).',obj.frz);

    % RAS parameters
    else
        % Frozen orbitals
        locFrz      =   obj.frz;
        
        % Active orbitals
        locAct      =   {setdiff(locAct{1},locFrz)};
    end


end

