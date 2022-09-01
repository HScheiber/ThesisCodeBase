% JC 2 x 2 matrix:   sigma_M      sigma_X                (units: nm)
%                    epsilon_M    epsilon_X              (units: kJ/mol)

% TF 4 x 3 array:   alpha_MM    alpha_XX    alpha_MX   (units: nm^-1)
%                   B_MM        B_XX        B_MX       (units: kJ/mol)
%                   C_MM        C_XX        C_MX       (units: (kJ nm^6)/mol)
%                   D_MM        D_XX        D_MX       (units: (kJ nm^8)/mol)
clear;
Alph_conv = 10; % Ang^-1 -> nm^-1
B_conv = 1e4; % 10^4 kJ -> kJ
C_conv = 0.1^6; % A^6 -> nm^6
D_conv = 0.1^8; % A^8 -> nm^8

% JC_Parameters = [1.409/10 4.022/10 (1.409+4.022)/20; ...
%               1.4088967 0.0309637 sqrt(1.4088967*0.0309637)];

% JC_Parameters = [0.11312014, 0.47257793 (0.11312014+0.47257793)/20; ...
%               1.76859616 0.04949224 sqrt(1.76859616*0.04949224)];

% TF_Parameters = [3.3445*Alph_conv  3.3445*Alph_conv 3.3445*Alph_conv; ...
%                 0.95535*B_conv    4.0616*B_conv    2.2115*B_conv; ...
%                 4.3962*C_conv     873.21*C_conv    48.177*C_conv; ...
%                 1.8066*D_conv     1023.8*D_conv    36.133*D_conv];
Settings = Initialize_MD_Settings;
Startpoint = 0;
Settings.Table_Length = 4.01;
Settings.Table_StepSize = 0.0005;
Settings.Salt = 'LiI';
Settings.Structure = 'Rocksalt';
Settings.Theory = 'TF';
Settings.JobName = 'Test';
Settings.WorkDir = pwd;
Settings.MDP.vdw_modifier = 'potential-shift';
Settings.MDP.RVDW_Cutoff = 1.9;
Settings.TF_Paramset = 0;
Settings.Parallel_Min = false;
Settings.MDP.Maintain_Symmetry = true;


% Allowed plot types: 'full', 'full-derivative','lj', 'lj-derivative',
% 'dispersion', 'dispersion-derivative', 'repulsive',
% 'repulsive-derivative'
PlotType = 'full';
plot_PES = true;
Calc_energy = false;
save_table = false;

Settings.Model = 'ED5';% 'CR3';
Settings = Load_Model_Params(Settings);
%Settings.CR_Damp = Init_CRDamping_Object;
% CRDamping.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
% CRDamping.MX.b = 100; % sigmoid "steepness" for damping
% CRDamping.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
% CRDamping.MM.b  = 75; % 75
% CRDamping.XX.r_d = 0.20; 
% CRDamping.XX.b  = 100;

% CRDamping.MM.r_d = 0.1;
% CRDamping.MM.b = 300;
% CRDamping.XX.r_d = 0.2;
% CRDamping.XX.b = 300;
% CRDamping.MX.r_d = 0.13;
% CRDamping.MX.b = 300;


if plot_PES
    if strcmp(Settings.Theory,'JC')
        
        switch Settings.Theory
            case 'JC'
                Settings.WaterModel = 'SPC/E';
            case 'JC3P'
                Settings.WaterModel = 'TIP3P';
            case 'JC4P'
                Settings.WaterModel = 'TIP4PEW';
            case 'JCSD'
                Settings.WaterModel = 'SD';
        end
        
        JC_Potential_Parameters(Settings,'PlotType',PlotType,'Plotswitch',true,...
            'Startpoint',Startpoint);
        
        [U_MX, U_MM, U_XX] = JC_Potential_Generator(Settings,'PlotType',PlotType,'Plotswitch',true,...
            'Startpoint',Startpoint);
    elseif strcmp(Settings.Theory,'TF')

        TF_Potential_Parameters(Settings,'PlotType',PlotType,'Plotswitch',true,...
            'Startpoint',Startpoint);
        
        [U_MX, U_MM, U_XX] = TF_Potential_Generator(Settings,'Plotswitch',true,'PlotType',PlotType,...
            'Startpoint',Startpoint);
    elseif strcmp(Settings.Theory,'HS')
        [U_MX, U_MM, U_XX] = HS_Potential_Generator(Settings,'Plotswitch',true,'PlotType',PlotType,...
            'Startpoint',Startpoint);
        
    elseif strcmp(Settings.Theory,'BH')
        BH_Potential_Parameters(Settings,'PlotType',PlotType,'Plotswitch',true,...
            'Startpoint',Startpoint);
        
        [U_MX, U_MM, U_XX] = BH_Potential_Generator(Settings,'Plotswitch',true,'PlotType',PlotType,...
            'Startpoint',Startpoint);
    end
    
    if save_table
        [Settings.Metal,Settings.Halide] = Separate_Metal_Halide(Settings.Salt);
        
        TableName = [Settings.JobName '_Table'];
        TableFile_MX = fullfile(Settings.WorkDir,[TableName '.xvg']);
        TableFile_MM = fullfile(Settings.WorkDir,[TableName '_' Settings.Metal '_' Settings.Metal '.xvg']);
        TableFile_XX = fullfile(Settings.WorkDir,[TableName '_' Settings.Halide '_' Settings.Halide '.xvg']);

        % Save tables into current directory
        fidMX = fopen(TableFile_MX,'wt');
        fwrite(fidMX,regexprep(U_MX,'\r',''));
        fclose(fidMX);

        fidMM = fopen(TableFile_MM,'wt');
        fwrite(fidMM,regexprep(U_MM,'\r',''));
        fclose(fidMM);

        fidXX = fopen(TableFile_XX,'wt');
        fwrite(fidXX,regexprep(U_XX,'\r',''));
        fclose(fidXX);
    end
end

if Calc_energy
    output = Structure_Minimization(Settings);
end