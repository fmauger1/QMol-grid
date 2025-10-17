function analyzeCIspectrum(obj,n)
%analyzeCIspectrum CI spectrum analysis
%   After calculating the CI matrix, use analyzeCIspectrum to print out the
%   energies and wave function composition of eigen states of the CI model.
%
%   obj.analyzeCIspectrum(n) displays the analysis for the n lowest-energy
%   eigen states.

    % Initialization ======================================================
    QMol_doc.showSection('CI matrix spectrum');

    fprintf('  * CI matrix\n')
    fprintf('    size = %i x %i with %i nonzero elements\n',...
                [size(obj.CI) sum(abs(obj.CI) > 0,'all')])

    function plotComposition(V)
        for k = 1:numel(V), a = V(k)^2; if a > 0.01 %#ok<ALIGN>
            fprintf('    > %#7.3f %% [ ',100*a)
            fprintf(' %i ',obj.CSB(k,:))
            fprintf(' ]\n')
        end, end
    end

    % Calculate spectrum ==================================================
    if size(obj.CI,1) < 500   ||   (n > 0.5*size(obj.CI,1))
        % Full diagonalization
        [V,E]               =   eig(obj.CI);
        [E,ind]             =   sort(diag(E));                              % sort energies
        V                   =   V(:,ind);
    else
        % Partial diagonalization
        [V,E]               =   eigs(obj.CI,n,'smallestreal');
        E                   =   diag(E);
    end

    % Ground state configuration ==========================================
    fprintf('  * Ground state\n')
    fprintf('    Total energy      = %#10.3f a.u. = %#10.3f eV\n',E(1),convertUnit.au2ev(E(1)))
    plotComposition(V(:,1));

    % Excited state configuration =========================================
    for l = 2:n
        fprintf('  * Excited state\n')
        fprintf('    Total energy      = %#10.3f a.u. = %#10.3f eV\n',E(l),convertUnit.au2ev(E(l)))
        fprintf('    Excitation energy = %#10.3f a.u. = %#10.3f eV\n',E(l)-E(1),convertUnit.au2ev(E(l)-E(1)))
        plotComposition(V(:,l));
    end

    % Finalization ========================================================
    QMol_doc.showSection();

end

