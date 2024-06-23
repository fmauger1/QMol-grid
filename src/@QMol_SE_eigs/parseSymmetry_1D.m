function S = parseSymmetry_1D(~,symList)
%parseSymmetry_1D parse the symmetry-configuration list for 1D DFT systems
    
    % Initialization
    symList             =   strip(strsplit(symList,'+'));

    % Parse symmetry groups
    S                   =   zeros(2,numel(symList));

    for k = 1:numel(symList)
        % Detect symmetry condition
        if contains(symList{k},'sx','IgnoreCase',true)
            S(2,k)      =   1;
            symList{k}  = strip(regexprep(symList{k},'sx','','ignorecase'));
        elseif contains(symList{k},'ax','IgnoreCase',true)
            S(2,k)      =  -1;
            symList{k}  = strip(regexprep(symList{k},'ax','','ignorecase'));
        else
            warning('QMolGrid:SE_eigs:symmetry', ...
                ['Unable to identify the symmetry condition for ' symList{k} '. Symmetry group ignored'])
        end
        
        % Detect how may states
        if isempty(symList{k})
            S(1,k)      =   1;
        else
            S(1,k)      =   str2double(symList{k});
        end
    end
end