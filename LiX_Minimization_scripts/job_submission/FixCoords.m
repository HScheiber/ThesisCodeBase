function Output = FixCoords(Input)
    Output = Input;
    while Output < 0
        Output = Output + 1;
    end
    while Output >= 1.0
        Output = Output - 1;
    end
    return
end