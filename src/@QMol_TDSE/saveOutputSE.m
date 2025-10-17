function saveOutputSE(obj,t)
%saveOutputSE
    
    % Extract relevalt data
    SE                  =   obj.SE;
    externalField       =   obj.EF;

    % Save data
    save([obj.sSEF '_' num2str(obj.oSE.n,'%06u') '.mat'],'SE','t','externalField');

    % Update counter
    obj.oSE.n           =   obj.oSE.n + 1;

end