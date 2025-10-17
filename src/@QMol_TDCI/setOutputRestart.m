function setOutputRestart(obj,opt)
%setOutputRestart

switch lower(opt)
case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Clean up any old data
    obj.oRest           =   [];
    
    % Time sampling
    if obj.sRest,       obj.oRest.ind   =   [obj.getOutputIndex(obj.sRestT),NaN];   if obj.oRest.ind(1) == 1 %#ok<ALIGN> 
                        obj.oRest.ind   =   obj.oRest.ind(2:end);                   end     % No restart on initial time
    else,               obj.oRest.ind   =   NaN;                                    end
    obj.oRest.n     =   1;

case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Don't need oRest anymore
    obj.oRest           =   [];

    % Delete the restart file
    if exist(obj.sRestF,'file')
        N   =   [];     save(obj.sRestF,'N'),   delete(obj.sRestF);         % In case the restart file is sent to the recycle bin
    end

otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unexpected option
    error('QMol:TDCI:setOutputRestart',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end