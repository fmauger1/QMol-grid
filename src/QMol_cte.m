classdef (Abstract) QMol_cte < QMol_suite
%QMol_cte common constants used by several components of QMol-grid
%   This class provides common constants shared by various components of
%   the QMol-grid package. End users are not expected to interact directly
%   with this class and should ignore it.
%
%   WARNING: editing the values for the constants listed in this class will
%     lead to hard-to-track side effects throughout the QMol-grid package.
%     Only do this if you are absolutely certain that it is required.
    
%   Version     Date        Author
%   01.23.000   05/25/2025  F. Mauger
%       Creation 

%% Documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static,Access=private)
function version
    QMol_doc.showVersion('01.23.000','05/25/2025','F. Mauger')
end
end
methods (Static,Access={?QMol_doc,?QMol_cte})
function showInfo
    fprintf('  * QMol_cte:\n');
    fprintf('      > Common constants, including\n'); 
    fprintf('        + Symplectic split-operator propagators\n');
    QMol_cte.version;
end
end
%% Symplectic propagators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties (Hidden,Constant)
    % Strang splitting (a.k.a Verlet) scheme, 2-component split
    symp_2S2_sC1      =   0.5
    symp_2S2_sC2      =   0.5
    symp_2S2_sD1      =   0.5
    % Strang splitting (a.k.a Verlet) scheme, n-component split
    symp_2Sn_a1        =   0.5
    symp_2Sn_a2        =   0.5
    % Forest-Ruth scheme, 2-component split
    symp_4FR2_sC1       =   .5               / (2-2^(1/3))
    symp_4FR2_sC2       =   .5 * (1-2^(1/3)) / (2-2^(1/3))
    symp_4FR2_sC3       =   .5 * (1-2^(1/3)) / (2-2^(1/3))
    symp_4FR2_sC4       =   .5               / (2-2^(1/3))
    symp_4FR2_sD1       =   1                / (2-2^(1/3))
    symp_4FR2_sD2       =          -2^(1/3)  / (2-2^(1/3))
    symp_4FR2_sD3       =   1                / (2-2^(1/3))
    % Forest-Ruth scheme, n-component split
    symp_4FRn_a1        =   .5 / (2 - 2^(1/3))
    symp_4FRn_a2        =   .5 / (2 - 2^(1/3))
    symp_4FRn_a3        =   .5 - 1 / (2 - 2^(1/3))
    symp_4FRn_a4        =   .5 - 1 / (2 - 2^(1/3))
    symp_4FRn_a5        =   .5 / (2 - 2^(1/3))
    symp_4FRn_a6        =   .5 / (2 - 2^(1/3))
    % 4th order Blanes and Moan, 2-component split
    symp_4BM2_sC1       =   0.0792036964311957
    symp_4BM2_sD1       =   0.209515106613362
    symp_4BM2_sC2       =   0.353172906049774
    symp_4BM2_sD2       =  -0.143851773179818
    symp_4BM2_sC3       =  -0.0420650803577195
    symp_4BM2_sD3       =   0.434336666566456
    symp_4BM2_sC4       =   0.2193769557534996
    symp_4BM2_sD4       =   0.434336666566456
    symp_4BM2_sC5       =  -0.0420650803577195
    symp_4BM2_sD5       =  -0.143851773179818
    symp_4BM2_sC6       =   0.353172906049774
    symp_4BM2_sD6       =   0.209515106613362
    symp_4BM2_sC7       =   0.0792036964311957
    % 4th order Blanes and Moan, n-component split
    symp_4BMn_a1        =   0.0792036964311957
    symp_4BMn_a2        =   0.1303114101821663
    symp_4BMn_a3        =   0.2228614958676077
    symp_4BMn_a4        =  -0.3667132690474257
    symp_4BMn_a5        =   0.3246481886897062
    symp_4BMn_a6        =   0.1096884778767498
    symp_4BMn_a7        =   0.1096884778767498
    symp_4BMn_a8        =   0.3246481886897062
    symp_4BMn_a9        =  -0.3667132690474257
    symp_4BMn_a10       =   0.2228614958676077
    symp_4BMn_a11       =   0.1303114101821663
    symp_4BMn_a12       =   0.0792036964311957
    % 6th order Blanes and Moan, 2-component split
    symp_6BM2_sC1       =   0.0502627644003922
    symp_6BM2_sD1       =   0.148816447901042
    symp_6BM2_sC2       =   0.413514300428344
    symp_6BM2_sD2       =  -0.132385865767784
    symp_6BM2_sC3       =   0.0450798897943977
    symp_6BM2_sD3       =   0.067307604692185
    symp_6BM2_sC4       =  -0.188054853819569
    symp_6BM2_sD4       =   0.432666402578175
    symp_6BM2_sC5       =   0.54196067845078
    symp_6BM2_sD5       =  -0.016404589403618
    symp_6BM2_sC6       =  -0.7255255585086897
    symp_6BM2_sD6       =  -0.016404589403618
    symp_6BM2_sC7       =   0.54196067845078
    symp_6BM2_sD7       =   0.432666402578175
    symp_6BM2_sC8       =  -0.188054853819569
    symp_6BM2_sD8       =   0.067307604692185
    symp_6BM2_sC9       =   0.0450798897943977
    symp_6BM2_sD9       =  -0.132385865767784
    symp_6BM2_sC10      =   0.413514300428344
    symp_6BM2_sD10      =   0.148816447901042
    symp_6BM2_sC11      =   0.0502627644003922
    % 6th order Blanes and Moan, n-component split
    symp_6BMn_a1        =   0.050262764400392
    symp_6BMn_a2        =   0.098553683500650
    symp_6BMn_a3        =   0.314960616927694
    symp_6BMn_a4        =  -0.447346482695478
    symp_6BMn_a5        =   0.492426372489876
    symp_6BMn_a6        =  -0.425118767797691
    symp_6BMn_a7        =   0.237063913978122
    symp_6BMn_a8        =   0.195602488600053
    symp_6BMn_a9        =   0.346358189850727
    symp_6BMn_a10       =  -0.362762779254345
    symp_6BMn_a11       =  -0.362762779254345
    symp_6BMn_a12       =   0.346358189850727
    symp_6BMn_a13       =   0.195602488600053
    symp_6BMn_a14       =   0.237063913978122
    symp_6BMn_a15       =  -0.425118767797691
    symp_6BMn_a16       =   0.492426372489876
    symp_6BMn_a17       =  -0.447346482695478
    symp_6BMn_a18       =   0.314960616927694
    symp_6BMn_a19       =   0.098553683500650
    symp_6BMn_a20       =   0.050262764400392
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

