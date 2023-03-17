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
    if ~isfield(Settings,'BestPointCriterion')
        Settings.BestPointCriterion = 'min-observed';
    end
    if ~isfield(Settings,'Initialize_From_Model_ExcludeError')
        Settings.InitializeExcludeError = false;
    end
    if ~isfield(Settings,'InitializeRealizeError')
        Settings.InitializeRealizeError = false;
    end
end