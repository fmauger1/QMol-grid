function wfcn = SE_allocateWaveFunction(obj,N,wfcn,randStr,isR)
%SE_allocateWaveFunction allocates a Schrodinger-equation wave function(s)
%   (wfcn) object with the specified number of wave functions. Optionally, 
%   if a wfcn object in given as input, it is (re)initialized instead.
    
    % Initialization
    if nargin < 5,  isR =   false;
    if nargin == 2, wfcn=   [];             end, end

    % Allocate wfcn object (to zero/rand for compatibility with resizing)
        % Allocate or resize wfcn(s) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if nargin < 4, if isempty(wfcn)                                     %#ok<ALIGN> 
            % Create object from scratch
            wfcn         =   QMol_SE_wfcn_basis('wfcn',zeros(size(obj.v,2),N));

        elseif size(obj.v,2) ~= size(wfcn.wfcn,1)
            % Incompatible sizes
            wfcn.set('wfcn',zeros(size(obj.v,2),N));
        else
            % Initialization
            S           =   size(wfcn.wfcn);

            % Resize (if needed)
            if S(2) > N(1),     wfcn.wfcn    =   wfcn.wfcn(:,1:N);
            else,               wfcn.wfcn    =  [wfcn.wfcn zeros(size(obj.v,2),N(1)-S(2))]; end,    end

        % Allocate random orbital ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else, if isempty(wfcn), wfcn         =   QMol_SE_wfcn_basis();                              end
            if isR,             wfcn.set('wfcn',rand(randStr,[size(obj.v,2),N])-.5);
            else,               wfcn.set('wfcn',rand(randStr,[size(obj.v,2),N],'like',1i)-.5-.5i);  end
            % Normalize orbitals
            wfcn.wfcn   =   obj.SE_normalizeWaveFunction(wfcn.wfcn);
        end

    % Initialize orbital
    wfcn.initialize(obj);

end