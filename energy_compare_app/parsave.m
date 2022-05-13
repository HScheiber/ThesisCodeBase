function parsave(fname,Data,DataName)
    if strcmp(DataName,'Data')
        save(fname,'Data')
    elseif strcmp(DataName,'Metal_BSSE_Data')
        Metal_BSSE_Data = Data;
        save(fname,'Metal_BSSE_Data')
    elseif strcmp(DataName,'Halide_BSSE_Data')
        Halide_BSSE_Data = Data;
        save(fname,'Halide_BSSE_Data')
    elseif strcmp(DataName,'BSSE_Corrected_Data')
        BSSE_Corrected_Data = Data;
        save(fname,'BSSE_Corrected_Data')
    end
end