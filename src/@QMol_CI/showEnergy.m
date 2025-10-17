function showEnergy(obj)
%showEnergy display the member wave function energy and error
%
%   obj.showEnergy 

    % Initialization
    [E,err]             =   obj.getEnergy;

    fprintf('  Wave fcn      Energy (-eV)         Error(a.u.)\n')
    fprintf('  --------     ------------          -----------\n');

    % Display results
    fprintf('    %3i        %9.3f             %#10.3e\n', ...
        [1:numel(E); convertUnit.au2ev(-E.'); err.']);

    % Finalization
    fprintf('  ----------------------------------------------\n');

end

