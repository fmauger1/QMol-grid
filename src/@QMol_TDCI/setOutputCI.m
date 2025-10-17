function setOutputCI(obj,opt)
%setOutputCI

% NOTE: does not call setOutputChildren

switch lower(opt)
case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Clean up any old data
    obj.oCI             =   [];
    
    % Time sampling
    if obj.sCI
        obj.oCI.ind     =   [obj.getOutputIndex(obj.sCIT) NaN];
        obj.oCI.n       =   1;
    else
        obj.oCI.ind     =   NaN;
        obj.oCI.n       =   1;
    end

case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Don't need oDFT anymore
    obj.oCI             =   [];

otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unexpected option
    error('QMol:TDCI:setOutputCI',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end