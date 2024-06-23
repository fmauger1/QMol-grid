function setExpV(obj,c)
%setExpV

    % Update parameters: c_ = [T,V]
    obj.c_(2)           =   c;

    % Compute potential
    obj.DFT.Vks         =   obj.DFT.getPotential([],obj.DFT.Vks);

    % Add CAP
    if ~isempty(obj.ABC),   obj.ABC.getPotential(obj.DFT.Vks);                  end

    % Update expV
    if obj.FG == 1,                                                             if obj.DFT.isSpinPol
        % Electric field
        obj.expVup      =   exp(-1i*c*obj.dt_*(obj.DFT.Vks.Vup+obj.FE*obj.X));
        obj.expVdw      =   exp(-1i*c*obj.dt_*(obj.DFT.Vks.Vdw+obj.FE*obj.X));  else
        obj.expV        =   exp(-1i*c*obj.dt_*(obj.DFT.Vks.V  +obj.FE*obj.X));  end
    else,                                                                       if obj.DFT.isSpinPol
        % No electric field
        obj.expVup      =   exp(-1i*c*obj.dt_*(obj.DFT.Vks.Vup));
        obj.expVdw      =   exp(-1i*c*obj.dt_*(obj.DFT.Vks.Vdw));               else
        obj.expV        =   exp(-1i*c*obj.dt_*(obj.DFT.Vks.V  ));               end
    end

end