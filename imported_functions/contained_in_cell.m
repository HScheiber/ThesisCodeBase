function [output,i] = contained_in_cell(str,pattern)
    output=false;
    for i=1:length(pattern)
        if strcmp(str,pattern{i})
            output=true;
            return
        end
    end
    i = [];
    return
end