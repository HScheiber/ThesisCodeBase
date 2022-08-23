
Settings.JobSettings.Cores = 8;
Settings.JobSettings.MPI_Ranks = 2;
Settings.JobSettings.OMP_Threads = 4;
Settings.JobSettings.dd  = [];
Settings.JobSettings.npme = [];

[Settings.home,Settings.project,Settings.computer,Settings.slurm,Settings.BO_Models,...
    Settings.qsub,Settings.passlog,Settings.pipe,Settings.wsl,~] = find_home;
[~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts,Settings.MLModelDir] = MD_Batch_Template(Settings.JobSettings);
Settings.WorkDir = pwd;
Settings.OuterDir = pwd;
% Settings.Output_Energies = 1;
% Settings.Calc_Energies = 1;
% Settings.Output_Coords = 1;
Settings.QECompressibility = 1e-7;
Settings.ScaleInitialLiqDensity = 0.8;
Settings.MDP.dt = 0.001;
Settings.Output_Coords = 1000;
Settings.Delete_Equil = false;
Settings.Verbose = true;
% Calc_Liquid_Properties_at_MP(Settings)
% Calc_Solid_Properties_at_MP(Settings)
Settings.Equilibrate_Liquid = 20;
Settings.CheckAmorphousLiquid = true;
Settings.AmorphousDiffThreshold = 1e-6;
Settings.Liquid_Test_Time = 100; % ps
Settings.Liquid_Equilibrate_Time = 25; % ps
Settings.Finite_T_Data = Initialize_Finite_T_Data(Settings);
%[Tm_estimate,WorkDir,Aborted,T_dat] = Find_Melting_Point(Settings);

Output = Calc_Liquid_Properties_at_MP(Settings);

[U_MX, U_MM, U_XX] = JC_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);

[U_MX, U_MM, U_XX] = TF_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);

[U_MX, U_MM, U_XX] = BH_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);