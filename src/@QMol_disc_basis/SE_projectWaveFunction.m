function proj = SE_projectWaveFunction(obj,wfcn,proj)
%SE_projectWaveFunction computes the projection coefficients associated 
%   with the (discretized) wave functions

    
    % Initialization
    if nargin == 2,     proj    =   [];                     end
    if isempty(proj),   proj    =   QMol_SE_wfcn_basis;     end
    
    proj.wfcn       =   NaN(obj.nV,size(wfcn.wfcn,2));

    % Compute projections
    for k = 1:obj.nV, if isreal(obj.v)                                  %#ok<ALIGN> 
        for l = 1:size(wfcn.wfcn,2),    proj.wfcn(k,l)  =   sum(     obj.v(:,k) .*wfcn.wfcn(:,l)) * obj.dx; end, else
        for l = 1:size(wfcn.wfcn,2),    proj.wfcn(k,l)  =   sum(conj(obj.v(:,k)).*wfcn.wfcn(:,l)) * obj.dx; end, end
    end

    % Initialize (mirror wave function input)
    if wfcn.isInit,      proj.initialize(obj);        end

end