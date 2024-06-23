function setOutputDFT(obj,opt)
%initializeOutputDFT

% NOTE: does not call setOutputChildren

switch lower(opt)
case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Clean up any old data
    obj.oDFT            =   [];
    
    % Time sampling
    if obj.sDFT
        obj.oDFT.ind    =   [obj.getOutputIndex(obj.sDFTT) NaN];
        obj.oDFT.n      =   1;
    else
        obj.oDFT.ind    =   NaN;
        obj.oDFT.n      =   1;
    end

case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Don't need oDFT anymore
    obj.oDFT            =   [];

otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unexpected option
    error('QMol:TDDFT:setOutputDFT',['Unknown option ' opt]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end