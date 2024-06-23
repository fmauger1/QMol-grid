function setOutputSE(obj,opt)
%initializeOutputSE

% NOTE: does not call setOutputChildren

switch lower(opt)
case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Clean up any old data
    obj.oSE             =   [];
    
    % Time sampling
    if obj.sSE
        obj.oSE.ind     =   [obj.getOutputIndex(obj.sSET) NaN];
        obj.oSE.n       =   1;
    else
        obj.oSE.ind     =   NaN;
        obj.oSE.n       =   1;
    end

case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Don't need oDFT anymore
    obj.oSE             =   [];

otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unexpected option
    error('QMol:TDSE:setOutputSE',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end