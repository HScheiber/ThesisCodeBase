
[Settings.home,Settings.project,Settings.computer,Settings.slurm,Settings.BO_Models,...
    Settings.qsub,Settings.passlog,Settings.pipe,Settings.wsl,~] = find_home;
[~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts,~] = MD_Batch_Template(Settings.JobSettings);
Settings.WorkDir = pwd;