classdef QMol_DFT_Vc_LDA_soft < QMol_suite
%QMol_DFT_Vc_LDA_soft implementation of the local-density-approximation
%   (LDA) correlation functional for a soft-Coulomb electron-electron
%   interaction potential.
%
%   The LDA correlation functional is defined for the soft-Coulomb
%   potential Vee(x) = 1 / sqrt( x^2 + 1 ), and parameterized following
%       N. Helbig, J.I. Fuks, M. Casula, M.J. Verstraete, M.A.L. Marques, 
%       I.V. Tokatly, and A. Rubio, "Density functional theory beyond the 
%       linear regime: Validating an adiabatic local density
%       approximation," Physical Review A 83, 032503 (2011).

%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release
%   01.21.001   07/10/2024  F. Mauger
%       Correct typo in run-time documentation

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.001','07/10/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_DFT_Vc_LDA_soft})
function showInfo
    fprintf( '  * QMol_DFT_Vc_LDA_soft:\n');
    fprintf(['      > LDA correlation\n',...
             '      > Soft-Coulomb electron interaction\n']);
    QMol_DFT_Vc_LDA_soft.version;
end
end
methods (Access=public)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Header
    fprintf('  * Correlation functional               local-density approximation (LDA)\n');
    fprintf('    for a one-dimensional (1D) soft-Coulomb electron-electron interaction\n');
    fprintf('    potential of the form Vee(x) = 1 / sqrt( x^2 + 1 ) [Helbig 2011]');
    ref                 =   {'Helbig 2011'};

    % SIC (if any)
    if obj.SIC == 0
        fprintf('.\n');
    elseif obj.SIC == 1
        fprintf(', with\n    an average-density self-interaction correction (ADSIC) [Legrand 2002].\n');
        ref             =   [ref, {'Legrand 2002'}];
    end

    % Version
    obj.version;

end
function mem = getMemoryProfile(~,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Display components
    mem                 =   0;

    if opt
        QMol_DFT_profiler.showMemoryFootprint('Correlation functional (LDA soft Coulomb)',mem,1);
    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess={?QMol_DFT_Vc_LDA_soft,?QMol_DFT})
    % Linked objects
    DFT                             % DFT.disc always defines the domain
end
properties (Transient,Hidden,GetAccess=?QMol_suite,SetAccess=private)
    % Run-time variables
    SIC                             % 0 = none, 1 = ADSIC
end
properties(Constant,Access=public)
    type                =   'LDA_C'
end
properties (Constant,Access=private)
    % Unpolarized component (B0 = 0)
    A0                  =   18.40
    C0                  =   7.501
    D0                  =   0.10185
    E0                  =   0.012827
    alpha0              =   1.511
    beta0               =   0.258
    m0                  =   4.424
    % Polarized component (B1 = 0)
    A1                  =   5.24
    C1                  =   1.568
    D1                  =   0.1286
    E1                  =   0.00320
    alpha1              =   0.0538
    beta1               =   1.56e-5
    m1                  =   2.958
    % Threshold density for zero potential
    tol                 =   1e-10
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Run-time variables
    obj.DFT             =   [];
    obj.SIC             =   [];
    
    % Initialization status
    obj.isInit          =   false;
end
function initialize(obj,DFT,SIC)
%initialize initializes the object
    
    % Initialization needed?
    obj.reset;
    
    % Set links
    obj.DFT             =   DFT;

    % Identify the flavor of SIC
    if nargin == 2,     SIC = [];       end
    if isempty(SIC),    SIC = 'none';   end

    switch lower(SIC)
        case {'none','off'}
            obj.SIC     =   0;
        case {'adsic','average density','averagedensity','average_density'}
            obj.SIC     =   1;
        otherwise
            warning('QMol:QMol_DFT_Vx_LDA_soft:SIC', ...
                    ['Unknown flavor of self-interaction correction (SIC) ' SIC '\n; SIC ignored for the LDA-exchange potential.'])
            obj.SIC     =   0;
    end
    
    % Miscellaneous
    obj.isInit          =   true;

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    ClassName           =   'QMol_DFT_Vx_LDA_soft';
    PropNames           =   {};
end
end
%% Get LDA exchange potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function V = getPotential(obj,rho,V,isAdd)
%getPotential returns the discretization of the LDA-correlation potential
%   associated with the member properties and input one-body density (rho).
%   Optionally, provide the potential object(s), where the discretization
%   should be stored, and if the exchange potential should be added to it.
    
    % Initialization
    if nargin < 4,  isAdd   =   false;                                      %#ok<ALIGN> 
    if nargin < 3,  V       =   [];
    if nargin < 2,  rho     =   [];             end, end, end
    
    if isempty(V)   ||   ~isAdd  
        V               =   obj.DFT.disc.DFT_allocatePotential(V);          % If any potential, reset it
    end
    
    if isempty(rho)
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);                % Use DFT one-body density storage
        rho             =   obj.DFT.rho;
    end

    % Compute the LDA-correlation potential
    N                   =   obj.DFT.Ntot;
    if obj.SIC == 1 % ADSIC ----------------
        if rho.isSpinPol
            % Spin polarized
            r_s         =   .5 ./ max(rho.rhoUp+rho.rhoDw,obj.tol);
            xi          =   (max(rho.rhoUp,obj.tol)-max(rho.rhoDw,obj.tol)) ...
                          ./(max(rho.rhoUp,obj.tol)+max(rho.rhoDw,obj.tol));
            V.add(-r_s .*( (1-xi.^2).*obj.d_eps_c0(r_s) + xi.^2 .* obj.d_eps_c1(r_s) ));
            V.add((1 - xi.*(2-xi)) .* obj.eps_c0(r_s) + xi.*(2-xi) .* obj.eps_c1(r_s), ...
                  (1 + xi.*(2+xi)) .* obj.eps_c0(r_s) - xi.*(2+xi) .* obj.eps_c1(r_s));
            
            r_s         =   .5 ./ max(rho.rhoUp/N(1),obj.tol);
            xi          =   .5 ./ max(rho.rhoDw/N(2),obj.tol);               % Reuse xi to prevent new memory allocation
            V.add(-obj.eps_c1(r_s) + r_s.*obj.d_eps_c1(r_s), ...
                  -obj.eps_c1(xi ) + xi .*obj.d_eps_c1(xi ));
        else
            % Spin restricted
            rs          =   .5 ./ max(rho.rho,obj.tol);
            V.add(obj.eps_c0(rs  ) -   rs .* obj.d_eps_c0(rs  ) ...
                 -obj.eps_c1(rs*N) + N*rs .* obj.d_eps_c1(rs*N));
        end
    else % No SIC --------------------------
        if rho.isSpinPol
            % Spin polarized
            r_s         =   .5 ./ max(rho.rhoUp+rho.rhoDw,obj.tol);
            xi          =   (max(rho.rhoUp,obj.tol)-max(rho.rhoDw,obj.tol)) ...
                          ./(max(rho.rhoUp,obj.tol)+max(rho.rhoDw,obj.tol));
            
            V.add(-r_s .*( (1-xi.^2).*obj.d_eps_c0(r_s) + xi.^2 .* obj.d_eps_c1(r_s) ));
            V.add((1 - xi.*(2-xi)) .* obj.eps_c0(r_s) + xi.*(2-xi) .* obj.eps_c1(r_s), ...
                  (1 + xi.*(2+xi)) .* obj.eps_c0(r_s) - xi.*(2+xi) .* obj.eps_c1(r_s));
        else
            % Spin restricted
            rs          =   .5 ./ max(rho.rho,obj.tol);
            V.add(obj.eps_c0(rs) - rs .* obj.d_eps_c0(rs));
        end
    end

end
function DV = getPotentialDerivative(obj,dim,rho,DV,isAdd)                   % We don't need D_rho
%getPotentialDerivative returns the discretization of the derivative of the
%   LDA-correlation potential associated with the input one-body density
%   (rho) and its derivative (D_rho). Optionally, provide the potential
%   object(s), where the discretization should be stored, and if the
%   exchange potential should be added to it.
%   The potential derivative is computed using fast-Fourier transforms
%   (since the potential is local).

    % Initialization
    if nargin < 5,  isAdd   =   false;                                      %#ok<ALIGN> 
    if nargin < 4,  DV      =   [];
    if nargin < 3,  rho     =   [];     end, end, end
    
    if isempty(DV)   ||   ~isAdd  
        DV              =   obj.DFT.disc.DFT_allocatePotentialGradient(DV); % If any potential, reset it
    end
    
    if isempty(rho)
        obj.DFT.rho     =   obj.DFT.getDensity(obj.DFT.rho);                % Use DFT one-body density storage
        rho             =   obj.DFT.rho;
    end

    % Compute the LDA-correlation potential derivative
    N                   =   obj.DFT.Ntot;
    if obj.SIC == 1 % ADSIC ----------------
        if rho.isSpinPol
            % Spin polarized
            r_s         =   .5 ./ max(rho.rhoUp+rho.rhoDw,obj.tol);
            xi          =   (max(rho.rhoUp,obj.tol)-max(rho.rhoDw,obj.tol)) ...
                          ./(max(rho.rhoUp,obj.tol)+max(rho.rhoDw,obj.tol));
            if dim == 1
                DV.add(1,real(ifft(obj.DFT.disc.D .* fft(-r_s .*( (1-xi.^2).*obj.d_eps_c0(r_s) + xi.^2 .* obj.d_eps_c1(r_s) ) + (1 - xi.*(2-xi)) .* obj.eps_c0(r_s) + xi.*(2-xi) .* obj.eps_c1(r_s)) )), ...
                         real(ifft(obj.DFT.disc.D .* fft(-r_s .*( (1-xi.^2).*obj.d_eps_c0(r_s) + xi.^2 .* obj.d_eps_c1(r_s) ) + (1 + xi.*(2+xi)) .* obj.eps_c0(r_s) - xi.*(2+xi) .* obj.eps_c1(r_s)) )));
            else
                warning('QMol:QMol_DFT_Vx_LDA_soft:getPotentialDerivative',['Unexpected dimension (' num2str(dim) ') for potential-derivative computation.'])
            end
            
            r_s         =   .5 ./ max(rho.rhoUp/N(1),obj.tol);
            xi          =   .5 ./ max(rho.rhoDw/N(2),obj.tol);               % Reuse xi to prevent new memory allocation
            if dim == 1
                DV.add(1,real(ifft(obj.DFT.disc.D .* fft(-obj.eps_c1(r_s) + r_s.*obj.d_eps_c1(r_s)) )), ...
                         real(ifft(obj.DFT.disc.D .* fft(-obj.eps_c1(xi ) + xi .*obj.d_eps_c1(xi )) )));
            else
                warning('QMol:QMol_DFT_Vx_LDA_soft:getPotentialDerivative',['Unexpected dimension (' num2str(dim) ') for potential-derivative computation.'])
            end
        else
            % Spin restricted
            rs          =   .5 ./ max(rho.rho,obj.tol);
            if dim == 1
                DV.add(1,real(ifft(obj.DFT.disc.D .* fft(obj.eps_c0(rs  ) -   rs .* obj.d_eps_c0(rs  ) -obj.eps_c1(rs*N) + N*rs .* obj.d_eps_c1(rs*N)) )));
            else
                warning('QMol:QMol_DFT_Vx_LDA_soft:getPotentialDerivative',['Unexpected dimension (' num2str(dim) ') for potential-derivative computation.'])
            end
        end
    else % No SIC --------------------------
        if rho.isSpinPol
            % Spin polarized
            r_s         =   .5 ./ max(rho.rhoUp+rho.rhoDw,obj.tol);
            xi          =   (max(rho.rhoUp,obj.tol)-max(rho.rhoDw,obj.tol)) ...
                          ./(max(rho.rhoUp,obj.tol)+max(rho.rhoDw,obj.tol));
            if dim == 1
                DV.add(1,real(ifft(obj.DFT.disc.D .* fft(-r_s .*( (1-xi.^2).*obj.d_eps_c0(r_s) + xi.^2 .* obj.d_eps_c1(r_s) ) + (1 - xi.*(2-xi)) .* obj.eps_c0(r_s) + xi.*(2-xi) .* obj.eps_c1(r_s)) )), ...
                         real(ifft(obj.DFT.disc.D .* fft(-r_s .*( (1-xi.^2).*obj.d_eps_c0(r_s) + xi.^2 .* obj.d_eps_c1(r_s) ) + (1 + xi.*(2+xi)) .* obj.eps_c0(r_s) - xi.*(2+xi) .* obj.eps_c1(r_s)) )));
            else
                warning('QMol:QMol_DFT_Vx_LDA_soft:getPotentialDerivative',['Unexpected dimension (' num2str(dim) ') for potential-derivative computation.'])
            end
            
        else
            % Spin restricted
            rs          =   .5 ./ max(rho.rho,obj.tol);
            if dim == 1
                DV.add(1,real(ifft(obj.DFT.disc.D .* fft(obj.eps_c0(rs) - rs .* obj.d_eps_c0(rs)) )));
            else
                warning('QMol:QMol_DFT_Vx_LDA_soft:getPotentialDerivative',['Unexpected dimension (' num2str(dim) ') for potential-derivative computation.'])
            end
        end
    end

end
function setPotentialKernel(~)
%setPotentialKernel no kernel to set
    
end
function E = getEnergy(obj,rho)
%getEnergy returns the exchange energy associated with the member 
%   properties and input one-body density (rho). Empty or missing density 
%   uses the one-body-density of the member DFT object (discouraged).
    
    % Initialization
    if nargin < 2,      rho     =   [];                                     end
    if isempty(rho),    rho     =   obj.DFT.getDensity(obj.DFT.rho);        end

    % Compute exchange energy
    N                   =   obj.DFT.Ntot;
    if obj.SIC == 1 % ADSIC ----------------
        if rho.isSpinPol
            % Spin polarized
            r_s         =   .5 ./ max(rho.rhoUp+rho.rhoDw,obj.tol);
            xi          =   (max(rho.rhoUp,obj.tol)-max(rho.rhoDw,obj.tol)) ...
                          ./(max(rho.rhoUp,obj.tol)+max(rho.rhoDw,obj.tol));
                      
            E           =   sum( (rho.rhoUp+rho.rhoDw) .* ...
                                ((1-xi.^2) .* obj.eps_c0(r_s) + xi.^2 .* obj.eps_c1(r_s)) );
            
            r_s         =   .5 ./ max(rho.rhoUp/N(1),obj.tol);
            E           =   E - sum( rho.rhoUp .* obj.eps_c1(r_s) );
            
            r_s         =   .5 ./ max(rho.rhoDw/N(2),obj.tol);
            E           =   E - sum( rho.rhoDw .* obj.eps_c1(r_s) );
        else
            % Spin restricted
            r_s         =   .5 ./ max(rho.rho,obj.tol);
            E           =   sum(rho.rho .* (obj.eps_c0(r_s)-obj.eps_c1(r_s*N)));
        end
    else % No SIC --------------------------
        if rho.isSpinPol
            % Spin polarized
            r_s         =   .5 ./ max(rho.rhoUp+rho.rhoDw,obj.tol);
            xi          =   (max(rho.rhoUp,obj.tol)-max(rho.rhoDw,obj.tol)) ...
                          ./(max(rho.rhoUp,obj.tol)+max(rho.rhoDw,obj.tol));
                      
            E           =   sum( (rho.rhoUp+rho.rhoDw) .* ...
                                ((1-xi.^2) .* obj.eps_c0(r_s) + xi.^2 .* obj.eps_c1(r_s)) );
        else
            % Spin restricted
            r_s         =   .5 ./ max(rho.rho,obj.tol);
            E           =   sum(rho.rho .* obj.eps_c0(r_s));
        end
    end

    E                   =   E * obj.DFT.disc.dx;
end
end
%% LDA exchange functional components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=private)
    % Energy per particle and derivative ---
    function E = eps_c0(obj,r),     E   =   -.5 * obj.eps_N0(r) ./ obj.eps_D0(r) .* obj.eps_L0(r);      end
    function DE = d_eps_c0(obj,r),  DE  =   -.5 * obj.d_eps_N0(r) ./ obj.eps_D0(r).^2 .* obj.eps_L0(r)  ...
                                            -.5 * obj.eps_N0(r) ./ obj.eps_D0(r) .* obj.d_eps_L0(r);    end
    
    function E = eps_c1(obj,r),     E   =   -.5 * obj.eps_N1(r) ./ obj.eps_D1(r) .* obj.eps_L1(r);      end
    function DE = d_eps_c1(obj,r),  DE  =   -.5 * obj.d_eps_N1(r) ./ obj.eps_D1(r).^2 .* obj.eps_L1(r)  ...
                                            -.5 * obj.eps_N1(r) ./ obj.eps_D1(r) .* obj.d_eps_L1(r);    end

    % Components ---------------------------
    function N = eps_N0(obj,r),     N   =   r + obj.E0*r.^2;                                end
    function D = eps_D0(obj,r),     D   =   obj.A0 + obj.C0*r.^2 + obj.D0*r.^3;             end
    function L = eps_L0(obj,r),     L   =   log( 1 + obj.alpha0*r + obj.beta0*r.^obj.m0 );  end
    
    function N = d_eps_N0(obj,r),   N   =   obj.A0 + 2*obj.A0*obj.E0*r - obj.C0*r.^2 - 2*obj.D0*r.^3 - obj.D0*obj.E0*r.^4;                  end
    function L = d_eps_L0(obj,r),   L   =   ( obj.alpha0 + obj.beta0*obj.m0*r.^(obj.m0-1) )./( 1 + obj.alpha0*r + obj.beta0*r.^obj.m0 );    end
    
    function N = eps_N1(obj,r),     N   =   r + obj.E1*r.^2;                                end
    function D = eps_D1(obj,r),     D   =   obj.A1 + obj.C1*r.^2 + obj.D1*r.^3;             end
    function L = eps_L1(obj,r),     L   =   log( 1 + obj.alpha1*r + obj.beta1*r.^obj.m1 );  end
    
    function N = d_eps_N1(obj,r),   N   =   obj.A1 + 2*obj.A1*obj.E1*r - obj.C1*r.^2 - 2*obj.D1*r.^3 - obj.D1*obj.E1*r.^4;                  end
    function L = d_eps_L1(obj,r),   L   =   ( obj.alpha1 + obj.beta1*obj.m1*r.^(obj.m1-1) )./( 1 + obj.alpha1*r + obj.beta1*r.^obj.m1 );    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

