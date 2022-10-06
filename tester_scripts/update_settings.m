
Settings.MDP.Initial_T = 1326.4;
Settings.T0 = 1326.4;
Settings.Target_T = 1326.4;

Settings.JobSettings.Cores = 8;
Settings.JobSettings.MPI_Ranks = 8;
Settings.JobSettings.OMP_Threads = 1;
Settings.JobSettings.dd  = [];
Settings.JobSettings.npme = [];

[Settings.home,Settings.project,Settings.computer,Settings.slurm,Settings.BO_Models,...
    Settings.qsub,Settings.passlog,Settings.pipe,Settings.wsl,~] = find_home;
Settings.scratch_dir = pwd;
[~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts,Settings.MLModelDir] = MD_Batch_Template(Settings.JobSettings);
%[WorkDir,Settings.JobName,Settings.Full_Model_Name] = GetMDWorkdir(Settings);
Settings.WorkDir = pwd;
Settings.OuterDir = pwd;
% Settings.Output_Energies = 1;
% Settings.Calc_Energies = 1;
% Settings.Output_Coords = 1;
%Settings.QECompressibility = 1e-7;
%Settings.ScaleInitialLiqDensity = 0.8;
%Settings.MDP.dt = 0.001;
%Settings.Output_Coords = 100;
%Settings.Delete_Equil = false;
Settings.Verbose = true;
Settings.Diary_Loc = '';
%Settings.CheckAmorphousHalide = true;
% Calc_Liquid_Properties_at_MP(Settings)
% Calc_Solid_Properties_at_MP(Settings)
% Settings.Equilibrate_Liquid = 20;
% Settings.CheckAmorphousLiquid = true;
% Settings.AmorphousDiffThreshold = 1e-6;
% Settings.Liquid_Test_Time = 100; % ps
% Settings.Liquid_Equilibrate_Time = 25; % ps
Settings.Equilibrate_Liquid = 20; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
Settings.Finite_T_Data = Initialize_Finite_T_Data(Settings);
Settings.Longest_Cutoff = max([Settings.MDP.RList_Cutoff Settings.MDP.RCoulomb_Cutoff Settings.MDP.RVDW_Cutoff]);
[Tm_estimate,WorkDir,Aborted,T_dat] = Find_Melting_Point(Settings);

OutputLiq = Calc_Liquid_Properties_at_MP(Settings);
OutputSol = Calc_Solid_Properties_at_MP(Settings);


Settings.MinMDP.Verbose = true;
Output = Structure_Minimization(Settings);


[U_MX, U_MM, U_XX] = JC_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);

[U_MX, U_MM, U_XX] = TF_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);

[U_MX, U_MM, U_XX] = BH_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);

[U_MX, U_MM, U_XX] = BD_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);

[U_MX, U_MM, U_XX] = BE_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);

[U_MX, U_MM, U_XX] = Mie_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);