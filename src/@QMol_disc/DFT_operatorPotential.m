function Hp = DFT_operatorPotential(obj,V,varargin)
%DFT_operatorPotential potential operator V | psi >
    
    % Apply potential operator
    if V.isSpinPol      % DFT_operatorPotential(obj,V,p,isUp,S)
        if nargin < 5,  Hp  =   V.applyPotential(varargin{:});      S = 0;
        else,           Hp  =   V.applyPotential(varargin{1:2});    S = varargin{3};    end
    else                % DFT_operatorPotential(obj,V,p,S)
        if nargin < 4,  Hp  =   V.applyPotential(varargin{:});      S = 0;
        else,           Hp  =   V.applyPotential(varargin{1});      S = varargin{2};    end
    end

    % Symmetry (only is needed)
    if S ~= 0,              Hp  =   obj.DFT_applySymmetry(S,Hp);            end
end