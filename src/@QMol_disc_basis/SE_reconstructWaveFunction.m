function wfcn = SE_reconstructWaveFunction(obj,proj,wfcn)
%SE_reconstructWaveFunction reconstructs the (discretization) of the 
%   wave functions associated with the proj projection coefficients
    
    % Initialization
    if nargin == 2,     wfcn    =   [];     end
    
    if isempty(wfcn)
        wfcn            =   QMol_SE_wfcn('wfcn',zeros(numel(obj.x),size(proj.wfcn,2)));
    else
        wfcn.wfcn       =   zeros(numel(obj.x),size(proj.wfcn,2));
    end

    % Reconstruct wave function
    for k = 1:obj.nV, for l = 1:size(proj.wfcn,2)                           %#ok<ALIGN> 
        wfcn.wfcn(:,l)  =   wfcn.wfcn(:,l) + proj.wfcn(k,l)*obj.v(:,k);
    end, end

    % Initialize (mirror projection input)
    if proj.isInit,     wfcn.initialize(obj);         end

end