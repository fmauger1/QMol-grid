function setExpV(obj,c)
%setExpV

    % Update parameters: c_ = [T,V]
    obj.c_(2)           =   c;

    % Update expV
    if obj.FG == 1,                                                                 if ~isempty(obj.ABC)
        % Electric field
        obj.expV        =   exp(-1i*c*obj.dt_*(obj.SE.V.V+obj.FE*obj.X+obj.ABC.V)); else
        obj.expV        =   exp(-1i*c*obj.dt_*(obj.SE.V.V+obj.FE*obj.X));           end
    else,                                                                           if ~isempty(obj.ABC)
        % No electric field
        obj.expV        =   exp(-1i*c*obj.dt_*(obj.SE.V.V             +obj.ABC.V)); else
        obj.expV        =   exp(-1i*c*obj.dt_*(obj.SE.V.V));                        end
    end

end