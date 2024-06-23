function Hp = DFT_operatorHamiltonian(obj,V,varargin)
%DFT_operatorHamiltonian DFT-Hamiltonian operator H | psi >
    
    % Apply Hamiltonian operator
    if isempty(obj.Tv), if obj.QM.isSpinPol %#ok<ALIGN>     % DFT_operatorHamiltonian(obj,V,p,isUp,S)
        if nargin < 5,  Hp  =   ifft(obj.T .* fft(varargin{1})) + V.applyPotential(varargin{:});    S = 0;
        else,           Hp  =   ifft(obj.T .* fft(varargin{1})) + V.applyPotential(varargin{1:2});  S = varargin{3};    end
    else                                                    % DFT_operatorHamiltonian(obj,V,p,S)
        if nargin < 4,  Hp  =   ifft(obj.T .* fft(varargin{1})) + V.applyPotential(varargin{:});    S = 0;
        else,           Hp  =   ifft(obj.T .* fft(varargin{1})) + V.applyPotential(varargin{1});    S = varargin{2};    end
    end, else,          if obj.QM.isSpinPol %#ok<ALIGN>     % DFT_operatorHamiltonian(obj,V,p,isUp,S)
        if nargin < 5,  Hp  =   ifft(obj.Tv.* fft(varargin{1})) + V.applyPotential(varargin{:});    S = 0;
        else,           Hp  =   ifft(obj.Tv.* fft(varargin{1})) + V.applyPotential(varargin{1:2});  S = varargin{3};    end
    else                                                    % DFT_operatorHamiltonian(obj,V,p,S)
        if nargin < 4,  Hp  =   ifft(obj.Tv.* fft(varargin{1})) + V.applyPotential(varargin{:});    S = 0;
        else,           Hp  =   ifft(obj.Tv.* fft(varargin{1})) + V.applyPotential(varargin{1});    S = varargin{2};    end
    end, end

    % Symmetry (only is needed)
    if S ~= 0,              Hp  =   obj.DFT_applySymmetry(S,Hp);            end
end