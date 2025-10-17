function setOutputExternalField(obj,sN,opt)
%setOutputExternalField

    switch lower(opt)
    case {'init','initialize','initialization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(obj.FA),    obj.(sN).potentialVector        =   NaN(numel(obj.FA), numel(obj.(sN).time));   end
        if ~isempty(obj.FE),    obj.(sN).electricField          =   NaN(numel(obj.FE), numel(obj.(sN).time));   end
        if ~isempty(obj.FDE),   obj.(sN).electricFieldDerivative=   NaN(numel(obj.FDE),numel(obj.(sN).time));   end
        
    case {'clean','finalize','finalization'} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % nothing to clean
     
    otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Unexpected option
        error('QMol:TDCI:setOutputExternalField',['Unknown option ' opt]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end