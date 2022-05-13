function output = Remove_Infinities(input_vector)
    W = fields(input_vector);
    for i=1:length(W)
        input_vector.(W{i})(isnan(input_vector.(W{i}))) = 0;
        input_vector.(W{i})(isinf(input_vector.(W{i}))) = 0;
    end

    output = input_vector;
end