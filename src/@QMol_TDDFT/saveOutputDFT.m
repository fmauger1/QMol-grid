function saveOutputDFT(obj,t)
%saveOutputDFT
    
    % Extract relevant data
    DFT                 =   obj.DFT;
    externalField       =   obj.EF;

    % Save data
    save([obj.sDFTF '_' num2str(obj.oDFT.n,'%06u') '.mat'],'DFT','t','externalField');

    % Update counter
    obj.oDFT.n          =   obj.oDFT.n + 1;

end