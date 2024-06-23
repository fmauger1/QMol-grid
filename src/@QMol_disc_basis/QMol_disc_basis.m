classdef QMol_disc_basis < QMol_disc
%QMol_disc Cartesian-grid domain discretization. It defines
%   > The domain discretization grid
%   > Common differential operators on that grid
%
% PROPERTIES
%     x (xspan) | QM | dx | D | T | Tv | nV (basisSize)
%     v (basis) | mD | mT | mTv
%     isBasis | isBasisPol
   
%   Version     Date        Author
%   01.21.000   06/17/2024  F. Mauger
%       Prepare 01.21 release

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.21.000','06/17/2024','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_disc})
function showInfo
    fprintf('  * QMol_disc_basis: Domain discretization\n      > Basis-set discretization\n'); 
    QMol_disc_basis.version;
end
end
methods (Access=?QMol_suite)
function ref = showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Initialization
    fprintf('  * Domain discretization                    basis set on a Cartesian grid\n');

    % Domain setting
    if ~obj.isInit
        fprintf('    Domain discretization has not be initialized\n')
    else
        N               =   factor(numel(obj.x));

        fprintf('    Grid  = %s:%s:%s\n',num2str(obj.x(1),3),num2str(obj.dx,3),num2str(obj.x(end),3));
        fprintf('    Size  = %u ',numel(obj.x));
        if isscalar(N)
            fprintf('(prime) points\n');
        else
            fprintf('(%u',N(1)); fprintf(' x %u',N(2:end)); fprintf(') points\n')
        end
        fprintf('    Basis = %u vectors\n',size(obj.v,2));
    end

    % Version
    obj.version;

    % References
    ref                 =   {};
end
end 
methods (Access=public)
function mem = getMemoryProfile(obj,opt)
%getMemoryProfile computes and returns an estimate of the total memory
%   footprint (in bytes) of the DFT object with all its components
%   initialized and used.
    
    % Initialization
    if nargin < 2,  opt =   false;  end

    % Parent memory profiling
    mem                 =   getMemoryProfile@QMol_disc(obj,opt);

    % Basis
    if isreal(obj.v),   m   =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.v),'real');
    else,               m   =   QMol_DFT_profiler.getMemoryFootprint(numel(obj.v),'imag');  end

    QMol_DFT_profiler.showMemoryFootprint('basis'       ,m,2);      mem     =   mem + m;
    
    % Operator matrices
    if isreal(obj.v),   m   =   QMol_DFT_profiler.getMemoryFootprint(size(obj.v,2)^2,'real');
    else,               m   =   QMol_DFT_profiler.getMemoryFootprint(size(obj.v,2)^2,'imag');  end

    QMol_DFT_profiler.showMemoryFootprint('gradient matrix'     ,m,2);
    QMol_DFT_profiler.showMemoryFootprint('Laplacian matrix'    ,m,2);
    mem                 =   mem + 2*m;

    m                   =   QMol_DFT_profiler.getMemoryFootprint(size(obj.v,2)^2,'imag');
    QMol_DFT_profiler.showMemoryFootprint('Laplacian matrix (velocity gauge)'    ,m,2);
    mem                 =   mem + m;

    % Orbital reconstruction
    if ~isempty(obj.QM)
        % DFT model
        if isa(obj.QM,'QMol_DFT')
            % Reconstructed orbitals
            obj.KSO     =   QMol_DFT_orbital('isSpinPol',obj.QM.isSpinPol);
            obj.KSO.disc=   obj;

            m           =   obj.KSO.getMemoryProfile(false);

            QMol_DFT_profiler.showMemoryFootprint('reconstructed orbitals'  ,m,2);
            mem         =   mem + m;

        % Nothing to do (other cases)
        end

    end

end
end
%% Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,GetAccess=public,SetAccess=?QMol_suite)
    v
end
properties (Dependent,GetAccess=public,SetAccess=?QMol_suite)
    basis                   % v
end
properties (Transient,GetAccess=public,SetAccess=private)
    % Matrix representation for the differential operators
    mD
    mT
    mTv
end
methods (Static,Access=public)
    function b = isBasis,       b   =   true;   end
    function b = isBasisPol,    b   =   false;  end
end
properties (Transient,Hidden,Access=?QMol_suite)
    KSO                     % reconstruction of the KSO on domain grid
end
properties (Constant,Access=private)
    eqTol               =   1e-10   % tolerance for equality of data
end
%% Alias handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    % x ~~~~~~~~~~~~~~~~
    function set.basis(obj,val),            obj.v       =   val;            end
    function val = get.basis(obj),          val         =   obj.v;          end
end
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function reset(obj)
%reset resets all temporary (transient) properties of the object. Should be
%   overloaded by each subclasses to perform proper reset actions.
    
    % Parent-class reset (also updates the isInit)
    reset@QMol_disc(obj);
    
    % Run-time variables
    obj.mD              =   [];
    obj.mT              =   [];
    obj.mTv             =   [];
    if ~isempty(obj.KSO)
        % Force delete previous reconstructed KSO on domain grid
        % (those are NOT the KSOs from the DFT object)
        delete(obj.KSO);    obj.KSO         =   [];
    end
end
function orthonormalizeBasis(obj,opt,ind)
%orthonormalizeBasis orthonormalize the member projection basis using the
%   specified scheme (default is Gram–Schmidt):
%   orthonormalizeBasis('overlap') uses a balanced orthonormalization
%     scheme, using the overlap matrix between the vectors.
%   orthonormalizeBasis('Gram–Schmidt') uses a Gram–Schmidt
%     orthonormalization algorithm. Optionally, include a second index
%     vector to specify the index ordering of the basis vectors in which
%     Gram–Schmidt procedure should be applied.
    
    % Initialization
    if nargin < 3,  ind =   [];
    if nargin < 2,  opt =   [];                 end, end

    n                   =   size(obj.v,2);
    if isempty(opt),opt =   'Gram–Schmidt';     end
    if isempty(ind),ind =   1:n;                end

    if ~obj.isInit, obj.dx= obj.x(2)-obj.x(1);  end     % object need not be initialized for basis orthonormalization
    
    % Orthonormalize basis
    switch lower(opt)
        case {'overlap','balanced'}
            % Compute overlap matrix
            M           =   NaN(n);

            for k = 1:n, for l = 1:k                                        %#ok<ALIGN> 
                M(k,l)  =   sum(conj(obj.v(:,k)).*obj.v(:,l)) * obj.dx;
                M(l,k)  =   conj(M(k,l));
            end, end
            
            % Orthonormalize basis
            [B,~]       =   eig(M);
            obj.v       =   obj.v * B;

            for k = 1: n   
                obj.v(:,k)= obj.v(:,k) / sqrt( sum(abs(obj.v(:,k)).^2) * obj.dx);
            end

        case {'gram–schmidt','gram schmidt','gramschmidt','gs'}
            % Gram–Schmidt orthonormalization
            for k = 1:n, for l = 1:k-1
                obj.v(:,ind(k)) =   obj.v(:,ind(k)) - obj.v(:,ind(l)) * sum(conj(obj.v(:,ind(l))).*obj.v(:,ind(k)))*obj.dx;     end
                obj.v(:,ind(k)) =   obj.v(:,ind(k)) / sqrt(sum(abs(obj.v(:,ind(k))).^2)*obj.dx);
            end

        otherwise
            warning('QMol:disc_basis:orthonormalizeBasis', ...
                ['Unknown orthonormalization method ' opt '; no orthonormalization performed.'])
    end

    % Reinitialize object (as needed)
    if obj.isInit,  obj.isInit  =   false;  obj.initialize(obj.QM);     end
end
function initialize(obj,QM)
%initialize initializes the object
    
    % Update links
    if nargin < 2,  QM  =   [];                 end
    
    % Initialization needed?
    if obj.isInit,   return; end

    % Parent-class initialization (also updates isInit)
    initialize@QMol_disc(obj,QM);

    % Matrix representation for the differential operators
    obj.nV              =   size(obj.v,2);
    obj.mD              =   NaN(obj.nV,obj.nV);
    obj.mT              =   NaN(obj.nV,obj.nV);

    for k = 1:obj.nV
        % Apply differential operator on basis
        Dv              =   ifft(obj.D.*fft(obj.v(:,k)));
        Tv              =   ifft(obj.T.*fft(obj.v(:,k)));
        if isreal(obj.v),   Tv  =   real(Tv);   Dv  =   real(Dv);   end
        
        % Compute matrix elements
        for l = 1:k
            obj.mD(l,k) =   sum(conj(obj.v(:,l)).*Dv)*obj.dx; obj.mD(k,l) =  -conj(obj.mD(l,k));
            obj.mT(l,k) =   sum(conj(obj.v(:,l)).*Tv)*obj.dx; obj.mT(k,l) =   conj(obj.mT(l,k));
        end
    end
end
function setTv(obj,A)
%setTv sets the Laplacian (kinetic) operator Tv, in the velocity gauge, for
%   the input potential vector A. Empty A clears Tv and mTv (reverting back
%   to the default T operator).
    
    % Initialization
    if nargin < 2,  A       =   [];                         end
    setTv@QMol_disc(obj,A);

    % Set mTv
    if isempty(obj.Tv),     obj.mTv     =   [];
    else,                   obj.mTv     =   NaN(obj.nV,obj.nV); 
        for k = 1:obj.nV,   Tv          =   ifft(obj.Tv.*fft(obj.v(:,k)));  %#ok<ALIGN> 
        for l = 1:k,        obj.mTv(l,k)=   sum(conj(obj.v(:,l)).*Tv)*obj.dx;
                            obj.mTv(k,l)=   conj(obj.mTv(l,k));             end, end
    end

end
end
methods (Static=true,Access=?QMol_suite)
function [ClassName,PropNames] = propertyNames()
%propertyNames returns the names of member properties that can be set
%   through name-value assignment
    
    % Parent-class components
    [~,PropNames]       =   QMol_disc.propertyNames;

    ClassName           =   'QMol_disc_basis';
    PropNames           =   [PropNames,{'v','basis'}];
end
end
%% Operator overload %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public)
function t = eq(disc1,disc2)
%eq overloads the == operator
    t                   =   numel(disc1.x) == numel(disc2.x)                        && ...  same number of elements
                            all(abs(disc1.x(:)-disc2.x(:)) <= disc1.eqTol,'all')    && ...  all elements are the same (within tolerance)
                            disc1.isBasis   &&   disc2.isBasis                      && ...  both are basis discretizations
                            all(size(disc1.v) == size(disc2.v),'all')               && ...  same number of basis elements
                            all(abs(disc1.v-disc2.v) <= disc1.eqTol,'all');               % same basis vectors (in the same order)
end
end
%% DFT specific methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public,Static)
    function s = DFT_classOrbital,              s   =   'QMol_DFT_orbital_basis';   end
    function s = DFT_classPotential,            s   =   'QMol_DFT_Vks_basis';       end
    function s = DFT_classPotentialGradient,    s   =   'N/A';                      end     % Not sure yet what's needed
end
methods (Access=public)
function DVks = DFT_allocatePotentialGradient(obj,DVks)                     %#ok<INUSD>       Not sure yet what's needed
%DFT_allocatePotentialGradient
    
    error('QMol:disc_basis:potentialGradient',...
        'Kohn-Sham potential gradient is not yet defined for basis-set discretizations');
end
end
methods (Access=public)
    % Kohn-Sham orbitals
    proj    =   DFT_projectOrbital(obj,KSO,proj)
    KSO     =   DFT_reconstructOrbital(obj,proj,KSO)
end
%% Schrodinger-equation specific methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access=public,Static)
    function s = SE_classWaveFunction,          s   =   'QMol_SE_wfcn_basis';       end
end
methods (Access=public)
    % Wave functions
    proj    =   SE_projectWaveFunction(obj,wfcn,proj)
    wfcn    =   SE_reconstructWaveFunction(obj,proj,wfcn)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

