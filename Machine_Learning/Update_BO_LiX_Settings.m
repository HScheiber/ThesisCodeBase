function Settings = Update_BO_LiX_Settings(Settings)
    if ~isfield(Settings,'GPActiveSetSize')
        Settings.GPActiveSetSize = 1000;
    end
    if ~isfield(Settings,'FinalGPFitActiveSetSize')
        Settings.FinalGPFitActiveSetSize = 6000;
    end
    if ~isfield(Settings,'OuterDir')
        Settings.OuterDir = pwd;
    end
end