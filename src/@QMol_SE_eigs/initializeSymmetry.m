function initializeSymmetry(obj)
%initializeSymmetry parses the symmetry configuration and sets the symState
%   member property accordingly
    
    if obj.SE.dim == 1 % 1D SE system ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if isempty(obj.sym),    obj.symState    =   [obj.SE.N;0];
        else,                   obj.symState    =   obj.parseSymmetry_1D(obj.sym);      end

    else % Unexpected case, throuw warning ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        warning('QMolGrid:SE_eigs:SEdim', ...
            ['State symmetry not (yet) supported for ' num2str(obj.SE.dim) 'D SE systems.'])
    end
end