% Required for parallelization
function Output = Crystal_Info(Cry,Structure,Item)
    if isempty(Item)
        Output = Cry.(Structure);
    else
        Output = Cry.(Structure).(Item);
    end
end