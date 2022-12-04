
Settings.MDP.Initial_T = 1183.9;
Settings.T0 = 1183.9;
Settings.Target_T = 1183.9;


% Settings.dd  = [];
% Settings.npme = [];
% Settings.Hours= 3;
% Settings.N_Calc= 8;
% Settings.Mins= 0;
% Settings.Nodes=0;
% Settings.Cores=12;
% Settings.Mempernode='0';
% Settings.SinglePrecision= 0;
% Settings.BigNode= 0;
% Settings.MPI_Ranks= 12;
% Settings.OMP_Threads= 1;
% Settings.npme= [];
% Settings.dd= [];
% Settings.dds= 0.8000;
% Settings.DLB= 1;
% Settings.TunePME= 1;

Settings.Cores = 12;
Settings.MPI_Ranks = 12;
Settings.OMP_Threads = 1;

% Settings.MP_Liquid_Test_Time = 100; % ps
% Settings.MP_Equilibrate_Solid = Settings.Equilibrate_Solid; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
% Settings.MP_Equilibrate_Liquid = 100; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
[Settings.home,Settings.project,Settings.computer,Settings.slurm,Settings.BO_Models,...
    Settings.qsub,Settings.passlog,Settings.pipe,Settings.wsl,~] = find_home;
Settings.scratch_dir = pwd;
[~,Settings] = MD_Batch_Template(Settings);
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
Settings.MP_Liquid_Test_Time = 2; %100 ps. Time used for calculation of liquid MSD in melting point calculations.
Settings.MP_Equilibrate_Solid = 2; %15 number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
Settings.MP_Equilibrate_Liquid = 2; %20 number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
Settings.Finite_T_Data = Initialize_Finite_T_Data(Settings);
Settings.Longest_Cutoff = max([Settings.MDP.RList_Cutoff Settings.MDP.RCoulomb_Cutoff Settings.MDP.RVDW_Cutoff]);
[Tm_estimate,WorkDir,Aborted,T_dat] = Find_Melting_Point(Settings);

OutputLiq = Calc_Liquid_Properties_at_MP(Settings);
OutputSol = Calc_Solid_Properties_at_MP(Settings);


Settings.MinMDP.Verbose = true;
Output = Structure_Minimization(Settings);


