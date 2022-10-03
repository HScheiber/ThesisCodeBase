% Script to test the reproducibility of Calc_Solid_Properties and
% Calc_Liquid_Properties
Settings = load('Calc_Settings.mat').Settings;

N_Reps = 50;
Settings.Verbose = true;
Settings.JobSettings.Cores = 8;
Settings.JobSettings.MPI_Ranks = 8;
Settings.JobSettings.OMP_Threads = 1;
Settings.JobSettings.dd  = [];
Settings.JobSettings.npme = [];
Settings.Diary_Loc='';

[Settings.home,Settings.project,Settings.computer,Settings.slurm,Settings.BO_Models,...
    Settings.qsub,Settings.passlog,Settings.pipe,Settings.wsl,~] = find_home;
[~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts,Settings.MLModelDir] = MD_Batch_Template(Settings.JobSettings);

workdir = pwd;
Output_Sol = cell(1,N_Reps);
Output_Liq = cell(1,N_Reps);
for idx = 1:N_Reps
    Settings.WorkDir = tempname(workdir);
    Settings.OuterDir = tempname(workdir);
    Settings.scratch_dir = tempname(workdir);
    Output_Liq{idx} = Calc_Liquid_Properties_at_MP(Settings);
    Settings.WorkDir = tempname(workdir);
    Settings.OuterDir = tempname(workdir);
    Settings.scratch_dir = tempname(workdir);
    Output_Sol{idx} = Calc_Solid_Properties_at_MP(Settings);
end

save(fullfile(workdir,'Reproducibility.mat'),'Output_Sol','Output_Liq');