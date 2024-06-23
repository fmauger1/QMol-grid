classdef (Abstract) fourierTool
%fourierTool miscellaneous tools for Fourier transforms

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
%   Version 01.02.001 - 07/18/2023 - Add the option to specify the number
%                                    of points in the fftGrid
%                 002              - Add the scaled fast Fourier transform

%% Version info (for dependent use) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Constant,Access=public)
    version             =   '01.02.002'
    lastMod             =   '07/18/2023'
    author              =   'F. Mauger'
end

%% Frequency grid for fft computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)
function freq = fftGrid(tspan,n)
    % Computation of the grid of frequencies associated with the Cartesian
    % (all increasing, equally-spaces) discretization of a signal with the
    % fast Fourier transform function fft. The output is a row vector with
    % the same number of elements as the Cartesian discretixation.
    % Optionally specify the number of frequency points to include in the
    % grid as a second argument.
    dt                  =   tspan(2)-tspan(1);
    
    if nargin < 2,  n   =   [];                 end
    if isempty(n),  n   =   length(tspan);      end
    if mod(n,2) == 1
        freq            =   ifftshift(2*pi/dt/n*( -(n-1)/2 : (n-1)/2 ));
    else
        freq            =   ifftshift(2*pi/dt/n*( - n   /2 : n/2-1 ));
    end
end
end
%% Scaled fast Fourier transforms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)
function Y = fft(X,dt,n,dim)
    % Computation of the scaled fast Fourier transform. The output is
    % provided on the same frequency grid as for native MATLAB FFT, but
    % with an amplitude scaled to match the definition
    %   F[X](w) = 1/sqrt(2*pi) * int(X(t)*exp(-i*w*t),-inf,inf)
    % WARNING: the scaling implicitly assumes that the input signal starts
    % at time t=0. If not, an additional phase should be added to account
    % for the time shift.

    % Initialization
    if nargin < 4,      dim =   [];
    if nargin < 3,      n   =   [];                         end, end
    
    if isempty(dim),    dim     =   find(size(X) > 1,1);    end             % 1st dimension greater than 1
    if isempty(n),      n       =   size(X,dim);            end
    
    % Sacled Fourier transform
    Y                   =   dt/sqrt(2*pi) * fft(X,n,dim);
    
end
function Y = ifft(X,df,n,dim)
    % Computation of the scaled inverse fast Fourier transform. The output
    % is provided on the same grid as for native MATLAB IFFT, but
    % with an amplitude scaled to match the definition
    %   F^-1[X](t) = 1/sqrt(2*pi) * int(X(w)*exp(i*w*t),-inf,inf)

    % Initialization
    if nargin < 4,      dim =   [];
    if nargin < 3,      n   =   [];                         end, end
    
    if isempty(dim),    dim     =   find(size(X) > 1,1);    end             % 1st dimension greater than 1
    if isempty(n),      n       =   size(X,dim);            end
    
    % Sacled Fourier transform
    Y                   =   n*df/sqrt(2*pi) * ifft(X,n,dim);
    
end
end
end

