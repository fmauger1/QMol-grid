classdef (Abstract) convertUnit
%convertUnit miscellaneous unit conversion

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
%   Version 01.01.000 - 07/23/2022 - Development version (from 01.00.003)
%   Version 01.02.000 - 09/30/2022

%% Version info (for dependent use) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Constant,Access=public)
    version             =   '01.02.000'
    lastMod             =   '09/30/2022'
    author              =   'F. Mauger'
end
%% Conversion from atomic units %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true)
function cval = au2ev(val)
    % Energy conversion from atomic units to electron volts
    %   1 [a.u.] = 27.211386245988 [eV]
    cval                =   val*27.211386245988;
end
function cval = au2kg(val)
    % Energy conversion from atomic units to kilograms
    %   1 [a.u.] = 9.1093837015e-31   [kg]
    cval                =   val*9.1093837015e-31;
end
function cval = au2meter(val)
    % Distance conversion from atomic units to meters
    %   1 [a.u.] = 5.29177210903e-11 [m]
    cval                =   val*5.29177210903e-11;
end
function cval = au2second(val)
    % Time conversion from atomic units to seconds
    %   1 [a.u] = 2.4188843265857e-17 [s]
    cval                =   val*2.4188843265857e-17;
end
end

%% Conversion from eV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true)
function cval = ev2au(val)
    % Energy conversion from electron volts to atomic units
    %   1 [eV] = 1/27.211386245988 [a.u.]
    cval                =   val/27.211386245988;
end
end

%% Conversion from kg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true)
function cval = kg2au(val)
    % Energy conversion from kilograms to atomic units
    %   1 [eV] = 1/9.1093837015e-31 [a.u.]
    cval                =   val/9.1093837015e-31;
end
end

%% From meter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true) 
function cval = meter2au(val)
    % Distance conversion from meters to atomic units
    %   1 [m] = 1/5.29177210903e-11 [a.u.]
    cval                =   val/5.29177210903e-11;
end
end

%% From second %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true)
function cval = second2au(val)
    % Time conversion from seconds to atomic units
    %   1 [s] = 1/2.4188843265857e-17 [a.u]
    cval                =   val/2.4188843265857e-17;
end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

