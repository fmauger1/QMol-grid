function setExpT(obj,d)
%setExpT
    
    % Update parameters: c_ = [T,V]
    obj.c_(1)           =   d;
    
    % Update expT
    if obj.FG == 2,         obj.DFT.disc.setTv(obj.FA);                     % Update velovity-Gauge kinetic operator
        obj.expT        =   exp(-1i*d*obj.dt_*obj.DFT.disc.Tv);             else
        obj.expT        =   exp(-1i*d*obj.dt_*obj.DFT.disc.T );
    end
end