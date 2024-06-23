classdef cubeObject
%cubeObject Volumetric data object, typically converted from OCTOPUS or
%   NWChem softwares.

% 2-Clause BSD License
%
% Copyright (c) 2022, Francois Mauger, all right reserved
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% F. Mauger
%   Version 01.01.000 - 07/23/2022 - Development version (from 01.00.002)
%   Version 01.02.000 - 09/30/2022 - Add exportCubeFile method
%                                    Correct importCubeFile method

%% Version info (for dependent use) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Constant,Access=public)
    version             =   '01.02.000'
    lastMod             =   '10/01/2022'
    author              =   'F. Mauger'
end

%% Constructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Access=public)
    % Volume data
    XGrid               =   []
    YGrid               =   []
    ZGrid               =   []
    
    XVector             =   []
    YVector             =   []
    ZVector             =   []
    
    Values              =   []
    
    % Time
    Time                =   []
    
    % Atomic centers
    AtomPosition        =   []
    AtomNumber          =   []
    AtomCharge          =   []
end
%% Methods (initialization) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Access = public)
function obj = cubeObject(varargin)
%cubeObject Constructor

    % If any property input, pass them to the set member method
    if nargin > 1,  obj = cubeObject.set(varargin{:});   end
end
end
methods (Static=true,Access=private)
function obj = set(varargin)
%setProperties sets named member properties to defined values

    % Create output default object
    obj                 =   cubeObject();
    
    % Test number of inputs
    if mod(nargin,2) == 0
        warning('cubeObject:setProperty:missingPropertyValue', ...
                'Inputs should be of maned property-value pair type. Lase entry ignored.')
    end
    
    % Set properties
    PN                  =  {'XGrid','YGrid','ZGrid','XVector','YVector','ZVector','Values', ...
                            'Time', ...
                            'AtomPosition','AtomNumber','AtomCharge'};
    
    for k = 1:floor(.5*nargin)
        % Identify property
        l               =   find(strcmpi(varargin{2*k-1},PN));
        
        % Set property to proper value
        if ~isempty(l)
            obj.(PN{l}) =   varargin{2*k};
        else
            warning('cubeObject:setProperty:propertyName', ...
                ['Unknown property ' varargin{2*k-1} ' -- entry ignored'])
        end
    end
    
end
end
%% Methods (import data from cube file) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true)
function obj = importCubeFile(fileName)
    % Create empty output object
    obj                 =   cubeObject();
    
    % Open input files
    FID                 =   fopen(fileName,'r');
    
    if FID < 0
        FID             =   fopen([fileName '.cube'],'r');
        
        if FID < 0
            try fclose(FID); end                                            %#ok<TRYNC>
            
            error('cubeObject:importCubeFile:inputFile', ...
                ['Unable to open ' fileName ' or ' fileName '.cube files for importing cube-file data.'])
        end
    end
    
    % First 2 (comment) lines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    com_1               =   fgetl(FID);
    com_2               =   fgetl(FID);
    
    if contains(lower(com_1),'octopus')                                     % Octopus generated
        sftwCase        =   1;
    elseif contains(lower(com_1),'nwchem')                                  % NWChem generated
        sftwCase        =   2;
    elseif contains(lower(com_1),'cubeObject')                              % cubeObject generated
        sftwCase        =   3;
    else
        sftwCase        =  -1;
    end
    
    if (sftwCase == 2 || sftwCase == 3) && strcmpi(com_2(1:4),' t =')       % Try getting time
        com_2           =   strsplit(com_2(5:end));
        
        if isempty(com_2{1})
            com_2       =   com_2(2:end);
        end
        
        obj.Time        =   str2double(com_2{1});
        
        switch lower(com_2{2})
            case {'au','a.u.'}
                            % Already in atomic units, nothing to do
            case 'fs'
                obj.Time=   convertUnit.second2au(obj.Time * 1e-15);
            otherwise
                warning('cubeObject:setProperty:timeUnits', ...
                    ['Unknown time units ' com_2{2} '; no conversion performed.'])
                
        end
        
    end
    
    % Number of atomic centers and cube properties ~~~~~~~~~~~~~~~~~~~~~~~~
    loc                 =   textscan(FID,'%d %f %f %f',4);
    
    X0                  =   [loc{2}(1) loc{3}(1) loc{4}(1)];
    
    if sftwCase == 1        % Octopus has the sign convention backward
        loc{1}(2:end)   =   -sign(loc{1}(2))*abs(loc{1}(2:end));
    end
    
    if loc{1}(2) > 0        % Already in atomic units
        obj.XVector   	=   [loc{2}(2) loc{3}(2) loc{4}(2)];
    else
        obj.XVector   	=   convertUnit.meter2au([loc{2}(2) loc{3}(2) loc{4}(2)] * 1e-10);
        X0(1)           =   convertUnit.meter2au(X0(1) * 1e-10);
    end
    obj.XGrid           =   X0(1) + double(0:abs(loc{1}(2))-1) * sqrt(sum(obj.XVector.^2));
    
    if loc{1}(3) > 0        % Already in atomic units
        obj.YVector   	=   [loc{2}(3) loc{3}(3) loc{4}(3)];
    else
        obj.YVector   	=   convertUnit.meter2au([loc{2}(3) loc{3}(3) loc{4}(3)] * 1e-10);
        X0(2)           =   convertUnit.meter2au(X0(2) * 1e-10);
    end
    obj.YGrid           =   X0(2) + double(0:abs(loc{1}(3))-1) * sqrt(sum(obj.YVector.^2));
    
    if loc{1}(4) > 0        % Already in atomic units
        obj.ZVector   	=   [loc{2}(4) loc{3}(4) loc{4}(4)];
    else
        obj.ZVector   	=   convertUnit.meter2au([loc{2}(4) loc{3}(4) loc{4}(4)] * 1e-10);
        X0(3)           =   convertUnit.meter2au(X0(3) * 1e-10);
    end
    obj.ZGrid           =   X0(3) + double(0:abs(loc{1}(4))-1) * sqrt(sum(obj.ZVector.^2));
    
    volUnits            =   sign(loc{1}(2:4)).';
    
    % Atomic centers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    loc                 =   textscan(FID,'%d %f %f %f %f',loc{1}(1));
    
    obj.AtomNumber      =   loc{1}.';
    if sftwCase == 1        % OCTOPUS doesn't seem to properly define the atomic/nuclear charge
        obj.AtomCharge  =   double(obj.AtomNumber);
    else
        obj.AtomCharge  =   loc{2}.';
    end
    
    if volUnits(1) < 0
        loc{3}          =   convertUnit.meter2au(loc{3} * 1e-10);
    end
    if volUnits(2) < 0
        loc{4}          =   convertUnit.meter2au(loc{4} * 1e-10);
    end
    if volUnits(3) < 0
        loc{5}          =   convertUnit.meter2au(loc{5} * 1e-10);
    end
    
    obj.AtomPosition    =   [loc{3} loc{4} loc{5}];
    
    % Volumetric data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nbX                 =   numel(obj.XGrid);
    nbY                 =   numel(obj.YGrid);
    nbZ                 =   numel(obj.ZGrid);
    
    obj.Values          =   NaN(nbY,nbX,nbZ);
    
    for kx=1:nbX
        for ky=1:nbY
            obj.Values(ky,kx,:) = cell2mat(textscan(FID,'%f',nbZ));
        end
    end

    % Convert density to a.u.
    if any(volUnits < 0)
        vol             =   (convertUnit.au2meter(1)*1e10)^sum(volUnits < 0);
        obj.Values      =   obj.Values * vol;
    end
    
    % Close file
    fclose(FID);
end
end
%% Methods (export data to cube file) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
function exportCubeFile(obj,fileName,varargin)
%exportCubeFile writes content of cubeObject into a plain-text cube file
%
%   optional name-value pairs:
%       units for time
%       units for space
    
    % Parse initial conditions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    uT              =   1;              % Time units (1 = a.u. and fs, 2 = a.u, 3 = fs)
    uL              =   1;              % Length units (1 = a.u., 2 = A)

    if mod(nargin,2) == 1
        warning('cubeObject:exportCubeFile:outputOptions',...
            'Expecting name-value pairs for output options; last option ignored.')
    end

    for k = 1:2:nargin-3,   switch lower(varargin{k})                       %#ok<ALIGN> 
        % Time stamp units
        case {'timeunits','time units','time_units','timeunit','time unit','time_unit'}
            switch lower(varargin{k+1})
                case {'a.u. and fs','a.u. & fs','a.u._and_fs','a.u._&_fs'}
                    uT  =   1;
                case {'a.u.'}
                    uT  =   2;
                case {'fs'}
                    uT  =   3;
                otherwise
                    warning('cubeObject:exportCubeFile:outputTimeUnits',...
                        ['Unknown units for time stamp ' varargin{k+1} '; default a.u. and fs printed instead.'])
                    uT  =   1;
            end
        
        % Length units
        case {'lengthunits','length units','length_units','lengthunit','length unit','length_unit'}
            switch lower(varargin{k+1})
                case {'a.u.'}
                    uL  =   1;
                case {'a','angstroms','angstrom'}
                    uL  =   2;
                otherwise
                    warning('cubeObject:exportCubeFile:outputLengthUnits',...
                        ['Unknown length units ' varargin{k+1} '; default a.u. used instead.'])
                    uL  =   1;
            end

        % Unknown option
        otherwise
            warning('cubeObject:exportCubeFile:outputOptions',...
                ['Unknown option ' varargin{k} '; option ignored.'])
    end, end
    
    if uL == 2, vol     =   convertUnit.meter2au(1e-10)^3;
    else,       vol     =   1;                                              end

    % Open output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if numel(fileName) > 5   &&   strcmpi(fileName(end-4:end),'.cube')
        FID             =   fopen(fileName,'w');
    else
        FID             =   fopen([fileName '.cube'],'w');
    end
    
    if FID < 0,     error('cubeObject:exportCubeFile:outputFile', ...
                        ['Unable to open ' fileName ' or ' fileName '.cube files for exporting cube-file data.'])
    end

    % Comment lines (header info) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf(FID,'Cube file generated with MATLAB cubeObject\n');
    if ~isempty(obj.Time)
        if     uT == 2,     fprintf(FID,'t = %14.7E au\n', obj.Time);
        elseif uT == 3,     fprintf(FID,'t = %14.7E fs\n', convertUnit.au2second(obj.Time)*1e15);
        else,               fprintf(FID,'t = %14.7E au = %14.7E fs\n', ...
                                obj.Time,convertUnit.au2second(obj.Time)*1e15);
        end
    else
        fprintf(FID,'(no time stamp)\n');
    end

    % Atomic centers and cube propreties ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    X0                  =   [obj.XGrid(1) obj.YGrid(1) obj.ZGrid(1)];
    if uL == 2;     X0  =   convertUnit.au2meter(X0)*1e10;                  end

    fprintf(FID,'%5u %+11s %+11s %+11s\n',numel(obj.AtomCharge),...
        num2str(X0(1),'%.6f'),num2str(X0(2),'%.6f'),num2str(X0(3),'%.6f'));
    

    N                   =   [numel(obj.XGrid) numel(obj.YGrid) numel(obj.ZGrid)];
    V                   =   [obj.XVector(:)'; obj.YVector(:)'; obj.ZVector(:)'];
    if uL == 2,     N   =   -N;     V   =   convertUnit.au2meter(V)*1e10;   end

    fprintf(FID,'%5i %+11s %+11s %+11s\n',N(1),num2str(V(1,1),'%.6f'),num2str(V(1,2),'%.6f'),num2str(V(1,3),'%.6f'));
    fprintf(FID,'%5i %+11s %+11s %+11s\n',N(2),num2str(V(2,1),'%.6f'),num2str(V(2,2),'%.6f'),num2str(V(2,3),'%.6f'));
    fprintf(FID,'%5i %+11s %+11s %+11s\n',N(3),num2str(V(3,1),'%.6f'),num2str(V(3,2),'%.6f'),num2str(V(3,3),'%.6f'));


    A                   =   obj.AtomPosition;
    if uL == 2;     A   =   convertUnit.au2meter(A)*1e10;                   end
    for k = 1:numel(obj.AtomCharge)
        fprintf(FID,'%5u %+11s %+11s %+11s %+11s\n',obj.AtomNumber(k),num2str(obj.AtomCharge(k),'%.6f'),...
            num2str(A(k,1),'%.6f'),num2str(A(k,2),'%.6f'),num2str(A(k,3),'%.6f'));
    end

    % Volumetric data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    nbX                 =   numel(obj.XGrid);
    nbY                 =   numel(obj.YGrid);
    nbZ                 =   numel(obj.ZGrid);

    for kx = 1:nbX, for ky = 1:nbY, for kz = 1:nbZ                          %#ok<ALIGN> 
        fprintf(FID,' %13.6E',vol*obj.Values(ky,kx,kz));
        if mod(kz,6) == 0,  fprintf(FID,'\n');              end
    end,    fprintf(FID,'\n');  end, end

    % Close file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fclose(FID);
end
end
end

