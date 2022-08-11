
[Settings.home,Settings.project,Settings.computer,Settings.slurm,Settings.BO_Models,...
    Settings.qsub,Settings.passlog,Settings.pipe,Settings.wsl,~] = find_home;
[~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts,~] = MD_Batch_Template(Settings.JobSettings);
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
Calc_Liquid_Properties_at_MP(Settings,'Verbose',true)
Calc_Solid_Properties_at_MP(Settings,'Verbose',true)



[U_MX, U_MM, U_XX] = JC_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);

[U_MX, U_MM, U_XX] = TF_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);

[U_MX, U_MM, U_XX] = BH_Potential_Generator(Settings,'Plotswitch',true,'PlotType','full',...
    'Startpoint',0.001);