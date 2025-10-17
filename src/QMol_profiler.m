classdef (Abstract) QMol_profiler < QMol_suite
%QMol_profiler parent memory and/or execution time profiling 
%   This class provides common formatting for the various profilers
%   associated with different simulation frameworks.
    
%   Version     Date        Author
%   01.23.000   05/25/2025  F. Mauger
%       Creation (from QMol_DFT_profiler version 01.21.000)

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.000','05/25/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_profiler})
function showInfo
    fprintf('  * QMol_profiler:\n');
    fprintf('      > Memory and execution time profiling interface\n'); 
    QMol_profiler.version;
end
end
%% Memory profiling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=?QMol_suite)
function m = getMemoryFootprint(N,opt)
%getMemoryFootprint obtain the memory footprint of an array with N elements
%   of the type opt. The memory footpring m is expressed in bites

    % Compute size
    switch lower(opt)
        case {'real','double'}
            m           =   8 * N;
        case {'imag','complex'}
            m           =   16 * N;
        otherwise
            warning('QMol:profiler:arrayType',['Unknown array type: ' opt '. No memory evaluation performed']);
            m           =   0;
    end

end
function showMemoryFootprint(msg,mem,lvl)
%showMemoryFootprint displays the memory foortrint 
    
    % Initialization
    if nargin < 3,  lvl = 1;        end

    % Display memory footprint
    switch lvl
        case 1
            fprintf('  * %-60s ',msg);      QMol_profiler.showMemoryFormatted(mem);
        case 2
            fprintf('    > %-58s ',msg);    QMol_profiler.showMemoryFormatted(mem);
        otherwise
            fprintf('    > %-59s ',msg);    QMol_profiler.showMemoryFormatted(mem);
            warning('QMol:profiler:showLevel','Maximum show level is 2; Entry downgraded to level 2');
    end
    fprintf('\n');
end
% end
% methods (Static,Access=private)
function showMemoryFormatted(mem)
%showMemoryFormatted display the formatted memory footprint
    
    if mem == 0   ||   isnan(mem)   || isinf(mem)                           % Do not display
    elseif mem <= 1024,       fprintf('%6u  B',mem);
    elseif mem <= 1024^2,   fprintf('%6.1f KB',mem/1024);
    elseif mem <= 1024^3,   fprintf('%6.1f MB',mem/1024^2);
    elseif mem <= 1024^4,   fprintf('%6.1f GB',mem/1024^3);
    else,                   fprintf('%6.1f TB',mem/1024^4);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

