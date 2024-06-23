function showSymmetry(obj)
%showSymmetry displays the symmetry configuration (for showDocumentation)
    
    if isempty(obj.SE) % No SE object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Nothing to show
    elseif obj.SE.dim == 1 % 1D SE system ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Initialization
        msg             =   {'Ax','','Sx'};
        if size(obj.symState,2)==1 && obj.symState(2,1)==0, return; end % No symmetry 
        
        fprintf('    State sym. = %i %s',obj.symState(1,1),msg{obj.symState(2,1)+2});
        
        for k = 2:size(obj.symState,2)
            fprintf(' + %i %s',obj.symState(1,k),msg{obj.symState(2,k)+2});
        end

        % Final line break
        fprintf('\n');

    else % Unexpected case, throuw warning ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        warning('QMolGrid:SE_eigs:SEdim', ...
            ['Display of state symmetry not (yet) supported for ' num2str(obj.DFT.dim) 'D SE systems.'])
    end
end