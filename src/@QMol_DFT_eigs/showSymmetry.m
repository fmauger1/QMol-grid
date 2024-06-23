function showSymmetry(obj)
%showSymmetry displays the symmetry configuration (for showDocumentation)
    
    if isempty(obj.DFT) % No DFT object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Nothing to show
    elseif obj.DFT.dim == 1 % 1D DFT system ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Initialization
        msg             =   {'Ax','','Sx'};
        
        if obj.DFT.isSpinPol
            % Initialization
            Sup         =   size(obj.symState{1},2)==1 && obj.symState{1}(2,1)==0;
            Sdw         =   size(obj.symState{2},2)==1 && obj.symState{2}(2,1)==0;

            if Sup && Sdw,  return;     end     % No symmetry on either spin component
            
            fprintf('    State sym. = ');

            % Up-spin component
            fprintf('(up) ');

            if Sup
                fprintf('no symmetry')
            else
                fprintf('%i %s',obj.symState{1}(1,1),msg{obj.symState{1}(2,1)+2});
            end

            for k = 2:size(obj.symState{1},2)
                fprintf(' + %i %s',obj.symState{1}(1,k),msg{obj.symState{1}(2,k)+2});
            end

            % Down-spin component
            fprintf(' | (down) ');

            if Sdw
                fprintf('no symmetry')
            else
                fprintf('%i %s',obj.symState{2}(1,1),msg{obj.symState{2}(2,1)+2});
            end

            for k = 2:size(obj.symState{2},2)
                fprintf(' + %i %s',obj.symState{2}(1,k),msg{obj.symState{2}(2,k)+2});
            end

            % Final line break
            fprintf('\n');

        else
            % Initialization
            if size(obj.symState,2)==1 && obj.symState(2,1)==0, return; end % No symmetry 
            
            fprintf('    State sym. = %i %s',obj.symState(1,1),msg{obj.symState(2,1)+2});
            
            for k = 2:size(obj.symState,2)
                fprintf(' + %i %s',obj.symState(1,k),msg{obj.symState(2,k)+2});
            end

            % Final line break
            fprintf('\n');

        end

    else % Unexpected case, throuw warning ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        warning('QMolGrid:DFT_eigs:DFTdim', ...
            ['Display of state symmetry not (yet) supported for ' num2str(obj.DFT.dim) 'D DFT systems.'])
    end
end