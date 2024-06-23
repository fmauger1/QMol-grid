classdef (Abstract) convertLaser
%convertLaser miscellaneous conversion for laser fields

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
%   Version 01.02.000 - 09/30/2022

%% Version info (for dependent use) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Constant,Access=public)
    version             =   '01.02.000'
    lastMod             =   '09/30/2022'
    author              =   'F. Mauger'
end

%% Electric field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true)
function cval = amplitude2intensity(val,ept)
    % Conversion of the electric field amplitude E0 with ellipticity ept,
    % expressed in atomic units, to the corresponding intensity I0,
    % expressed in watt per squared centimeters
    %   I0 [W/cm^2] = (1+ept^2)*(E0 [a.u.]/5.33803e-9).^2
    % If unspecified, the default ellipticity corresponds to linear
    % polarization (ept = 0).
    if nargin < 2, ept = 0; end
    cval                =   (1+ept.^2).*(val/.00000000533803).^2;
end
function cval = intensity2amplitude(val,ept)
    % Conversion of the electric field intensity I0 with ellipticity ept,
    % expressed in watt per squared centimeters, to the corresponding
    % amplitude E0, expressed in atomic units
    %   E0 [a.u.] = 5.33803e-9*sqrt(I0 [W/cm^2]/(1+ept^2)).
    % If unspecified, the default ellipticity corresponds to linear
    % polarization (ept = 0).
    if nargin < 2, ept = 0; end
    cval                =   .00000000533803*sqrt(val./(1+ept.^2));
end
end

%% Wavelength/frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static=true)
function cval = frequency2wavelength(val)
    % Conversion of the electric field frequency omega, expressed in atomic
    % units, to the corresponding wavelength lambda, expressed in meters.
    %   lambda [m] = 2*pi*299792458/4.1341e16/omega [a.u.]
    cval                =   2*pi*299792458/41341000000000000./val;
end
function cval = wavelength2frequency(val)
    % Conversion of the electric field wavelength lambda, expressed in
    % meters, to the corresponding frequency omega, expressed in atomic
    % units.
    %   omega [a.u.] = 2*pi*299792458/4.1341e16/ lambda [m]
    cval                =   2*pi*299792458/41341000000000000./val;
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

