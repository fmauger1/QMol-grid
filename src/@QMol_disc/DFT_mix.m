function X = DFT_mix(obj,X,varargin)
%DFT_mix linear mixing of input pairs of weights and DFT-data objects (all
%   one-body densities of all Kohn-Sham potentials)
    
    if isa(X,'QMol_DFT_density')
        % Mixing densities -----------------
        if obj.QM.isSpinPol
            X.rhoUp     =   varargin{1}(1)*varargin{2}.rhoUp;
            X.rhoDw     =   varargin{1}(2)*varargin{2}.rhoDw;
            for k = 5:2:nargin
                X.rhoUp =   X.rhoUp + varargin{k-2}(1)*varargin{k-1}.rhoUp;
                X.rhoDw =   X.rhoDw + varargin{k-2}(2)*varargin{k-1}.rhoDw;
            end
        else
            X.rho       =   varargin{1}   *varargin{2}.rho;
            for k = 5:2:nargin
                X.rho   =   X.rho   + varargin{k-2}   *varargin{k-1}.rho;
            end
        end
    elseif isa(X,'QMol_DFT_Vks')
        % Mixing potentials ----------------
        if obj.QM.isSpinPol
            X.Vup       =   varargin{1}(1)*varargin{2}.Vup;
            X.Vdw       =   varargin{1}(2)*varargin{2}.Vdw;
            for k = 5:2:nargin
                X.Vup   =   X.Vup + varargin{k-2}(1)*varargin{k-1}.Vup;
                X.Vdw   =   X.Vdw + varargin{k-2}(2)*varargin{k-1}.Vdw;
            end
        else
            X.V         =   varargin{1}   *varargin{2}.V;
            for k = 5:2:nargin
                X.V     =   X.V   + varargin{k-2}   *varargin{k-1}.V;
            end
        end
    else
        warning('QMol:disc:DFT_mix', ...
                ['Cannot mix ' class(X) '; no mixing performed.'])
    end

end