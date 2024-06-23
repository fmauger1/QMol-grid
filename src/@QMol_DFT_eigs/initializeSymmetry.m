function initializeSymmetry(obj)
%initializeSymmetry parses the symmetry configuration and sets the symState
%   member property accordingly
    
    if obj.DFT.dim == 1 % 1D DFT system ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if obj.DFT.isSpinPol
            if isempty(obj.sym),    obj.symState    =  {[numel(obj.DFT.occ{1});0], ...
                                                        [numel(obj.DFT.occ{2});0]};         else %#ok<ALIGN> 
            
            % Spin up
            if isempty(obj.sym{1}), Sup             =   [numel(obj.DFT.occ{1});0];
            else,                   Sup             =   obj.parseSymmetry_1D(obj.sym{1});   end

            % Spin down
            if isempty(obj.sym{2}), Sdw             =   [numel(obj.DFT.occ{2});0];
            else,                   Sdw             =   obj.parseSymmetry_1D(obj.sym{2});   end

            % Put it all together
                                    obj.symState    =  {Sup,Sdw};                           end
        else % Spin restricted
            if isempty(obj.sym),    obj.symState    =   [numel(obj.DFT.occ);0];
            else,                   obj.symState    =   obj.parseSymmetry_1D(obj.sym);      end
        end

    else % Unexpected case, throuw warning ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        warning('QMolGrid:DFT_eigs:DFTdim', ...
            ['State symmetry not (yet) supported for ' num2str(obj.DFT.dim) 'D DFT systems.'])
    end
end