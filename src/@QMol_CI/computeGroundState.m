function computeGroundState(obj,n)
%computeGroundState compute the CI ground and excited states
%   Use computeGroundState to calculate the ground and first excited states
%   of the CImatrix. The result stored in the waveFunction property. The CI
%   object must be initialized and contain a CImatrix.
%   
%   obj.computeGroundState or obj.computeGroundState(1) calculates the
%     ground state only.
%
%   obj.computeGroundState(n), where n is a positive integer, calculates
%     the ground and n-1 lovest-energy excited states.

    % Initialization
    if nargin == 1,     n   =   [];     end
    if isempty(n),      n   =   1;      end

    % Calculate the states
    if ~issparse(obj.CI)  && (size(obj.CI,1) < 500   ||   (n > 0.5*size(obj.CI,1)))
        % Full diagonalization
        [obj.wfcn,E]        =   eig(obj.CI);
        [~,ind]             =   sort(diag(E));                              % sort energies
        obj.wfcn            =   obj.wfcn(:,ind(1:n));
    else
        % Partial diagonalization
        [obj.wfcn,~]        =   eigs(obj.CI,n,'smallestreal');
    end
end

