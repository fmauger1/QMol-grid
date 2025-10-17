function saveOutputCI(obj,t)
%saveOutputCI
    
    % Extract relevalt data
    CI                  =   obj.CI;
    externalField       =   obj.EF;

    % Save data
    save([obj.sCIF '_' num2str(obj.oCI.n,'%06u') '.mat'],'CI','t','externalField');

    % Update counter
    obj.oCI.n           =   obj.oCI.n + 1;

end