% Script to test the reproducibility of Calc_Solid_Properties and
% Calc_Liquid_Properties
Settings = load('Calc_Settings.mat').Settings;

N_Reps = 1000;
Settings.Verbose = true;
Settings.Cores = 8;
Settings.MPI_Ranks = 8;
Settings.OMP_Threads = 1;
Settings.dd  = [];
Settings.npme = [];
Settings.Diary_Loc='';

[Settings.home,Settings.project,Settings.computer,Settings.slurm,Settings.BO_Models,...
    Settings.qsub,Settings.passlog,Settings.pipe,Settings.wsl,~] = find_home;
[~,Settings] = MD_Batch_Template(Settings);

workdir = pwd;
load(fullfile(workdir,'Reproducibility.mat'),'Output_Sol','Output_Liq');
for idx = 206:N_Reps
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