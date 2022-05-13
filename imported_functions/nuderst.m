function out_vec = nuderst(args)

    if ~isnumeric(args)
        error('Input must be numeric.')
    end
    
    out_vec = nan(size(args));
    for ii = 1:size(args,1)
        for jj = 1:size(args,2)
            out_vec(ii,jj) = max(abs(args(ii,jj))*1e-4,1e-7);
        end
    end
end