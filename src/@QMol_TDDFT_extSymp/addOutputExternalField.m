function addOutputExternalField(obj,sN,~,~)
%addOutputExternalField

    % Add external-field information to output structure
    if ~isempty(obj.FA),    obj.(sN).potentialVector(:,obj.(sN).n)          =   obj.FA;     end
    if ~isempty(obj.FE),    obj.(sN).electricField(:,obj.(sN).n)            =   obj.FE;     end
    if ~isempty(obj.FDE),   obj.(sN).electricFieldDerivative(:,obj.(sN).n)  =   obj.FDE;    end
end