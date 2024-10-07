function values = getParamValues(params,paramLabels,variableParam,iter)
    values = zeros(1,length(params.values));
    for i = 1:length(params.values)
        if paramLabels(i) == variableParam
            paramValueCell = params(variableParam);
            values(i) = paramValueCell{1}(iter);
        else
            valueToAdd = params(paramLabels(i));
            values(i) = valueToAdd{1};
        end
        
    end
end