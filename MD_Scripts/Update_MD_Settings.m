function Settings = Update_MD_Settings(Settings)
    if ~isfield(Settings,'Metal')
        [Settings.Metal,Settings.Halide] = Separate_Metal_Halide(Settings.Salt);
    end
    if ~isfield(Settings,'Longest_Cutoff')
        Settings.Longest_Cutoff = max([Settings.MDP.RList_Cutoff Settings.MDP.RCoulomb_Cutoff Settings.MDP.RVDW_Cutoff]);
    end
    if ~isfield(Settings,'gmx')
        [~,Settings] = MD_Batch_Template(Settings);
    end
    if strcmp(Settings.MDP.CoulombType,'PME') && Settings.GaussianCharge
        Settings.MDP.CoulombType = 'PME-User';
    end
    if ~isfield(Settings,'Metal')
        [Settings.Metal,Settings.Halide] = Separate_Metal_Halide(Settings.Salt);
    end
    if ~isfield(Settings,'Comb_rule')
        Settings.Comb_rule = 'Lorentz-Berthelot';
    end
end