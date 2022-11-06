Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
T0s = [1121.3 883.2 825.2 742.2]; % Initial temperature (experimental MP)

%% BH models from GG and GI set
Theory = 'BH';
Models = {'GG1' 'GG2' 'GI2' 'GI3' };
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    T0 = T0s(sidx);
    
    if strcmp(Salt,'LiF')
        Shared_Settings.Liquid_Fraction = 0.53;
    else
        Shared_Settings.Liquid_Fraction = 0.51;
    end
    
    for midx = 1:length(Models)
        Model = Models{midx};
        
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = '10K-CA2'; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).N_atoms = 10000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
        Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = T0; % K, Initial temperature
        Settings_array(idx).c_over_a = 2;
    end
end

%% JC models from CM/CR set
% For LiF : test CM3, CM5 (fixed charge), CR2, CR3 (no fixed charge)
% For LiCl : test CM3, CM4
% For LiBr : test CM8, CM10
% For LiI : test CM3, CM4
Models = {{'CM3' 'CM5' 'CR2' 'CR3'} {'CM3' 'CM4'} {'CM8' 'CM10'} {'CM3' 'CM4'}};
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    Model_set = Models{sidx};
    T0 = T0s(sidx);
    
    if strcmp(Salt,'LiF')
        Shared_Settings.Liquid_Fraction = 0.53;
    else
        Shared_Settings.Liquid_Fraction = 0.51;
    end
    
    for midx = 1:length(Model_set)
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).Theory = 'JC'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model_set{midx}; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = '10K-CA2'; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).N_atoms = 10000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
        Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = T0; % K, Initial temperature
        Settings_array(idx).c_over_a = 2;
    end
end

%% Default JC/TF Models
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
T0s = [1121.3 883.2 825.2 742.2 1073.2]; % Initial temperature (experimental MP)
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    T0 = T0s(sidx);
    if strcmp(Salt,'LiF') || strcmp(Salt,'NaCl') 
        Shared_Settings.Liquid_Fraction = 0.53;
    else
        Shared_Settings.Liquid_Fraction = 0.51;
    end
        
    %% Default JC box 2:1
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = 'JC'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = '10K-CA2'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 10000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    Settings_array(idx).c_over_a = 2;
    
    %% JC with spherical interface (finite sized crystals)
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = 'JC'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = '10K-Bubble'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 50000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    Settings_array(idx).GenCluster = true;
    Settings_array(idx).Cluster_Frac = 0.1; % Cluster fraction, must be between (0,1)
    Settings_array(idx).Liquid_Interface = true;
    
    %% TF with spherical interface (finite sized crystals)
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = 'TF'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = '10K-Bubble'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 50000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    Settings_array(idx).GenCluster = true;
    Settings_array(idx).Cluster_Frac = 0.1; % Cluster fraction, must be between (0,1)
    Settings_array(idx).Liquid_Interface = true;
    
    %% Default JC box 2:1 Longer MaxCheckTime
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = 'JC'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = '10K-CA2-Long'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 10000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).MaxCheckTime = 2000; % ps

    %% Default JC def box 3:1
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = 'JC'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = '10K-CA3'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 10000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    Settings_array(idx).c_over_a = 3;
    
    %% Large JC def box 2:1
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = 'JC'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = '20K-CA2'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 20000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    Settings_array(idx).c_over_a = 2;
    
    %% Default TF box 2:1 BetaBeO
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = 'TF'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'BetaBeO'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = '10K-CA2'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 10000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    Settings_array(idx).c_over_a = 2;
    
    %% Default TF box 2:1 Rocksalt
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = 'TF'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = '10K-CA2'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 10000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    Settings_array(idx).c_over_a = 2;
    
end

%% Box Scaling Experiments
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
T0s = [1121.3 883.2 825.2 742.2 1073.2]; % Initial temperature (experimental MP)
Models = {{'JC' ''} {'TF' ''}};
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    T0 = T0s(sidx);
    if strcmp(Salt,'LiF') || strcmp(Salt,'NaCl') 
        Shared_Settings.Liquid_Fraction = 0.53;
    else
        Shared_Settings.Liquid_Fraction = 0.51;
    end
    
    for jdx = 1:length(Models)
        Theory = Models{jdx}{1};
        Model = Models{jdx}{2};
        
        %% 30 K atoms
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = '30K-CA2'; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).N_atoms = 30000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
        Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = T0; % K, Initial temperature
        Settings_array(idx).c_over_a = 2;
        
        %% 50 K atoms
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = '50K-CA2'; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).N_atoms = 50000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
        Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = T0; % K, Initial temperature
        Settings_array(idx).c_over_a = 2;

        %% Spherical interface 20 K atoms
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = '20K-Bubble'; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).N_atoms = 100000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
        Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = T0; % K, Initial temperature
        Settings_array(idx).GenCluster = true;
        Settings_array(idx).Cluster_Frac = 0.1; % Cluster fraction, must be between (0,1)
        Settings_array(idx).Liquid_Interface = true;
        
        %% Spherical interface 30 K atoms
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = '30K-Bubble'; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).N_atoms = 150000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
        Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = T0; % K, Initial temperature
        Settings_array(idx).GenCluster = true;
        Settings_array(idx).Cluster_Frac = 0.1; % Cluster fraction, must be between (0,1)
        Settings_array(idx).Liquid_Interface = true;
        
        %% Spherical interface 50 K atoms
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = '50K-Bubble'; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).N_atoms = 250000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
        Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = T0; % K, Initial temperature
        Settings_array(idx).GenCluster = true;
        Settings_array(idx).Cluster_Frac = 0.1; % Cluster fraction, must be between (0,1)
        Settings_array(idx).Liquid_Interface = true;
    end
end

%% JC EQ models
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theory = 'JC';
Models = {'EQ2' 'EQ2' 'EQ1' 'EQ5' };
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    T0 = T0s(sidx);
    Model = Models{sidx};
    
    if strcmp(Salt,'LiF')
        Shared_Settings.Liquid_Fraction = 0.53;
    else
        Shared_Settings.Liquid_Fraction = 0.51;
    end

    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = '20K-CA2'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 20000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    Settings_array(idx).c_over_a = 2;
end

%% BH FE models
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theory = 'BH';
Models = {'FE5' 'FE6' 'FE7' 'FE5' };
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    T0 = T0s(sidx);
    Model = Models{sidx};
    
    if strcmp(Salt,'LiF')
        Shared_Settings.Liquid_Fraction = 0.53;
    else
        Shared_Settings.Liquid_Fraction = 0.51;
    end
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = '20K-CA2'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 20000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    Settings_array(idx).c_over_a = 2;
end

%% Set1: vary the XY dimension
Salt = 'NaCl';
Models = {{'JC' ''} {'TF' ''}};
XY_Sizes = 4:14; % nm
Z_Size = 10; % nm
for kdx = 1:length(XY_Sizes)
    XY_Size = XY_Sizes(kdx);
    for jdx = 1:length(Models)
        Theory = Models{jdx}{1};
        Model = Models{jdx}{2};
        
        switch Theory
            case 'JC'
                T0 = 1073.2 + 210;
            case 'TF'
                T0 = 1073.2;
        end

        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).BracketThreshold = 5;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = ['Set1_XY_' num2str(XY_Size)  '_Z_' num2str(Z_Size)]; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = T0; % K, Initial temperature
        Settings_array(idx).Isotropy = 'semiisotropic';
        Settings_array(idx).Target_P = [1 1]; % Bar
        Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
        Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
        Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
        Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    end
end

%% Set2: vary the Z dimension
Salt = 'NaCl';
Models = {{'JC' ''} {'TF' ''}};
XY_Size = 6; % nm
Z_Sizes = [4.5 4.6 4.7 4.8 4.9 5 6 7 8 10 11 12 14]; % nm
for kdx = 1:length(Z_Sizes)
    Z_Size = Z_Sizes(kdx);
    for jdx = 1:length(Models)
        Theory = Models{jdx}{1};
        Model = Models{jdx}{2};
        
        switch Theory
            case 'JC'
                T0 = 1073.2 + 210;
            case 'TF'
                T0 = 1073.2;
        end
        
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).BracketThreshold = 5;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = ['Set2_XY_' num2str(XY_Size)  '_Z_' num2str(Z_Size)]; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = T0; % K, Initial temperature
        Settings_array(idx).Isotropy = 'semiisotropic';
        Settings_array(idx).Target_P = [1 1]; % Bar
        Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
        Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
        Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
        Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
        
        if Z_Size > 6
            Settings_array(idx).MeltFreezeThreshold = 1000;
        end
    end
end

%% Set3: Cluster, vary the radius R
Salt = 'NaCl';
Models = {{'JC' ''} {'TF' ''}};
R_Sizes = [1.9 2 2.1 2.2 2.3 2.4 2.5 3 4 5 6 7]; % nm
for kdx = 1:length(R_Sizes)
    R_Size = R_Sizes(kdx);
    for jdx = 1:length(Models)
        Theory = Models{jdx}{1};
        Model = Models{jdx}{2};
        
        switch Theory
            case 'JC'
                T0 = 1073.2 + 200;
            case 'TF'
                T0 = 1073.2;
        end
        
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        
        
        Settings_array(idx).N_Calc = 2; % Number of chained calculations
        Settings_array(idx).Hours = 12; % Max time for each job (hours)
        Settings_array(idx).Mins = 0; % Max time for job (minutes)
        Settings_array(idx).Nodes = 1; % Number of nodes to request for calculation.
        Settings_array(idx).Cores = 32; % Number of cores to request for calculation. Set to -1 for entire node
        Settings_array(idx).MPI_Ranks = 32; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
        Settings_array(idx).OMP_Threads = 1; % Set the number of OMP threads per MPI rank
        Settings_array(idx).npme = []; % Number of rank assigned to PME
        Settings_array(idx).dd = []; %[5 3 2] Domain decomposition
        
        Settings_array(idx).BracketThreshold = 5;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = ['Set3_R_' num2str(R_Size)]; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = T0; % K, Initial temperature
        Settings_array(idx).Isotropy = 'isotropic';
        Settings_array(idx).Target_P = 1; % Bar
        Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
        Settings_array(idx).GenCluster = true; % Switch from FLAT INTERFACE to SPHERICAL INTERFACE
        Settings_array(idx).Manual_Radius = R_Size; % Box radius in nm
        Settings_array(idx).Liquid_Fraction = 0.70;
        if R_Size < 3
            Settings_array(idx).MeltFreezeThreshold = 0.2;
        end
    end
end

%% Set4: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 1 ns with MeltFreezeThreshold = 0.15
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 210;
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:5;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    T0_i = T0 + (rand-1/2)*8; % add a random number between -4 and +4 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set4_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1;
    Settings_array(idx).MinStepSize = 0.25;
end

%% Set5: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 2 ps with MeltFreezeThreshold = 0.15
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 210;
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:5;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    T0_i = T0 + (rand-1/2)*8; % add a random number between -4 and +4 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).BracketThreshold = 5;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set5_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1;
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 2000; % ps. Max time for melting/freezing runs
end

%% Set6: Repeated Measurements [8x8x10 nm] to check for consistency at tmax = 2 ps with MeltFreezeThreshold = 0.15
Salt = 'NaCl'; 
Theory = 'JC';
T0 = 1073.2 + 210;
XY_Size = 8; % nm
Z_Size = 10; % nm
Reps = 1:5;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx))
    T0_i = T0 + (rand-1/2)*8; % add a random number between -4 and +4 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).BracketThreshold = 5;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set6_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1;
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 2000; % ps. Max time for melting/freezing runs
end

%% Set7: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 2 ps with MeltFreezeThreshold = 0.10
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:5;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).BracketThreshold = 5;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set7_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1;
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 2000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.10;
end

%% Set8: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 2 ps with MeltFreezeThreshold = 0.25
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:5;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set8_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1;
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 2000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Set9: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set9_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Set10: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 10 ns with MeltFreezeThreshold = 0.25
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set10_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 10000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Set11: Repeated Measurements [8x8x10 nm] to check for consistency at tmax = 10 ns with MeltFreezeThreshold = 0.25
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 8; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set11_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 10000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Set12: Repeated Measurements [6x6x7 nm] to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 6; % nm
Z_Size = 7; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set12_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Set13: Repeated Measurements [6x6x4.5 nm] to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 6; % nm
Z_Size = 4.5; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set13_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Set14: Repeated Measurements with 4000 atoms to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
% Reproducing results of: https://pubs.rsc.org/en/content/articlehtml/2019/cc/c9cc06177k
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set14_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.41; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.41; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.41; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 4;
    Settings_array(idx).N_atoms = 4000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Set15: Repeated Measurements with 2000 atoms to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
% Reproducing results of: https://aip.scitation.org/doi/10.1063/1.4745205
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+15)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set15_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Set16: Repeated Measurements with 1024 atoms to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
% Reproducing results of: https://aip.scitation.org/doi/10.1063/1.4745205
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+15)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set16_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.2; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.2; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.2; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 1024;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Set17: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
% Dispersion correction added
% Cutoff 1.9 nm
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+17)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set17_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
end

%% Set18: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
% Dispersion correction added
% Cutoff 1.4 nm
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+18)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set18_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
end

%% Set19: Repeated Measurements [6x6x14 nm] to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
% Dispersion correction NOT added
% Cutoff 1.9 nm
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 6; % nm
Z_Size = 14; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+19)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set19_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Set20: Reproducing lit value with 2000 atoms, tmax = 5 ns with MeltFreezeThreshold = 0.25
% Dispersion correction added
% Cutoff 1.4 nm
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+15)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    
    % MPI not usable here, change up the job settings
    Settings_array(idx).N_Calc = 2; % Number of chained calculations
    Settings_array(idx).Hours = 6; % Max time for each job (hours)
    Settings_array(idx).Mins = 0; % Max time for job (minutes)
    Settings_array(idx).Nodes = 0; % Number of nodes to request for calculation.
    Settings_array(idx).Cores = 8; % Number of cores to request for calculation. Set to -1 for entire node
    Settings_array(idx).MPI_Ranks = 1; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
    Settings_array(idx).OMP_Threads = 8; % Set the number of OMP threads per MPI rank
    Settings_array(idx).npme = []; % Number of rank assigned to PME
    Settings_array(idx).dd = []; %[5 3 2] Domain decomposition

    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set20_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 2.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = 50; % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 

    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = 50; %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'Ewald'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 2.5050e-05; %Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
    Settings_array(idx).MDP.Fourier_Spacing = 0.5; % [nm]
end

%% Set21: Reproducing lit value with 2000 atoms (except using PME), tmax = 5 ns with MeltFreezeThreshold = 0.25
% Dispersion correction added
% Cutoff 1.4 nm
Set = 21;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 2.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = 50; % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = 50; %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set22: Reproducing lit value with 2000 atoms (except using PME + Berendsen Barostat), tmax = 5 ns with MeltFreezeThreshold = 0.25
% Dispersion correction added
% Cutoff 1.4 nm
Set = 22;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Berendsen'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = 50; % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = 50; %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set23: Reproducing lit value with 2000 atoms (except using PME + V-Rescale Thermostat), tmax = 5 ns with MeltFreezeThreshold = 0.25
% Dispersion correction added
% Cutoff 1.4 nm
Set = 23;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 2.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = 50; % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = 10; %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set24: Reproducing lit value with 2000 atoms (except using PME + Reduced Fourierspacing), tmax = 5 ns with MeltFreezeThreshold = 0.25
% Dispersion correction added
% Cutoff 1.4 nm
Set = 24;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 2.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = 50; % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = 50; %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
    Settings_array(idx).MDP.Fourier_Spacing = 0.1; % used 0.1 in minimization. Default 0.12 nm. Grid dimensions in PME are controlled with fourierspacing
end

%% Set25: Testing Pre-Equilibration with my "default" parameters
% Dispersion correction not added
% Cutoff 1.9 nm
Set = 25;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290; % K
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+Set)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    Settings_array(idx).MDP.Disp_Correction = false; % Adds in long-range dispersion correction
    
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).Barostat = 'Parrinello-Rahman';
    
end

%% Set26: Reproducing lit value with 2000 atoms (except using PME), tmax = 5 ns with MeltFreezeThreshold = 0.25
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
Set = 26;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 2.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = 50; % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 5e-2;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = 50; %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set19: Repeated Measurements [6x6x14 nm] to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
% Dispersion correction NOT added
% Cutoff 1.9 nm
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 6; % nm
Z_Size = 14; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+19)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).N_Calc = 2; % Number of chained calculations
    Settings_array(idx).Hours = 6; % Max time for each job (hours)
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set19_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Set25: Testing Pre-Equilibration with my "default" parameters
% Dispersion correction not added
% Cutoff 1.9 nm
Set = 25;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290; % K
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+Set)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).N_Calc = 2; % Number of chained calculations
    Settings_array(idx).Hours = 6; % Max time for each job (hours)
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    Settings_array(idx).MDP.Disp_Correction = false; % Adds in long-range dispersion correction
    
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).Barostat = 'Parrinello-Rahman';
    
end

%% Set27: 2000 atoms, Prequilibration, v-rescale/PR, tau_t = 0.2 ps, Testing tau_p = 5 ps
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 27;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1285;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 5.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set28: 2000 atoms, Prequilibration, v-rescale/PR, tau_t = 0.2 ps, Testing tau_p = 2 ps
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 28;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 2.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set29: 2000 atoms, Prequilibration, v-rescale/PR, tau_t = 0.2 ps, Testing tau_p = 1 ps
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 29;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction during minimization
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set30: 2000 atoms, Prequilibration, v-rescale/PR, tau_t = 0.2 ps, Testing tau_p = 0.5 ps
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 30;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 0.5; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set31: 2000 atoms, Prequilibration, v-rescale/PR, tau_t = 0.2 ps, Testing tau_p = 0.1 ps
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 31;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 0.1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set32: 2000 atoms, Prequilibration, v-rescale/PR, tau_p = 1 ps, Testing tau_t = 10 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 32;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 10; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set33: 2000 atoms, Prequilibration, v-rescale/PR, tau_p = 1 ps, Testing tau_t = 5 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 33;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 5; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set34: 2000 atoms, Prequilibration, v-rescale/PR, tau_p = 1 ps, Testing tau_t = 2 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 34;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set35: 2000 atoms, Prequilibration, v-rescale/PR, tau_p = 1 ps, Testing tau_t = 1 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 35;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 1; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set36: 2000 atoms, Prequilibration, v-rescale/PR, tau_p = 1 ps, Testing tau_t = 0.1 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 36;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.1; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set37: 2000 atoms, Prequilibration, v-rescale/PR, tau_p = 1 ps, Testing tau_t = 0.01 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 37;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.01; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set38: 2000 atoms, Prequilibration, v-rescale/PR, tau_p = 1 ps, Testing tau_t = 0.001 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 38;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.001; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set39: 2000 atoms, Prequilibration, nose-hoover/PR, tau_p = 1 ps, Testing tau_t = 10 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 39;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 10; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set40: 2000 atoms, Prequilibration, nose-hoover/PR, tau_p = 1 ps, Testing tau_t = 5 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 40;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 5; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set41: 2000 atoms, Prequilibration, nose-hoover/PR, tau_p = 1 ps, Testing tau_t = 2 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 41;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set42: 2000 atoms, Prequilibration, nose-hoover/PR, tau_p = 1 ps, Testing tau_t = 1 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 42;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 1; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set43: 2000 atoms, Prequilibration, nose-hoover/PR, tau_p = 1 ps, Testing tau_t = 0.1 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 43;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.1; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set44: 2000 atoms, Prequilibration, nose-hoover/PR, tau_p = 1 ps, Testing tau_t = 0.01 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 44;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.01; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set45: 2000 atoms, Prequilibration, nose-hoover/PR, tau_p = 1 ps, Testing tau_t = 0.001 ps,
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 45;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.001; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set46: 2000 atoms, Prequilibration, nose-hoover/PR, tau_t = 1.0, Testing tau_p = 5 ps
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 46;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1285;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 5.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 1.0; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set47: 2000 atoms, Prequilibration, nose-hoover/PR, tau_t = 1.0, Testing tau_p = 2 ps
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 47;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 2.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 1.0; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set48: 2000 atoms, Prequilibration, nose-hoover/PR, tau_t = 1.0, Testing tau_p = 1 ps
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 48;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 1.0; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set49: 2000 atoms, Prequilibration, nose-hoover/PR, tau_t = 1.0, Testing tau_p = 0.5 ps
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 49;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 0.5; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 1.0; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set50: 2000 atoms, Prequilibration, nose-hoover/PR, tau_t = 1.0, Testing tau_p = 0.1 ps
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 50;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 25; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 25; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 0.1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 1.0; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.

end

%% Set51: Testing Pre-Equilibration + Initial equilibration with my "default" parameters
% Dispersion correction not added
% Cutoff 1.9 nm
% Now the compressibility is the same in all spatial directions
% Additional 1ps relaxation added to start of simulation
Set = 51;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290; % K
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+Set)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).N_Calc = 2; % Number of chained calculations
    Settings_array(idx).Hours = 6; % Max time for each job (hours)
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    Settings_array(idx).MDP.Disp_Correction = false; % Adds in long-range dispersion correction
    
    Settings_array(idx).MP_Equilibrate_Solid = 10; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 1; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    Settings_array(idx).Barostat = 'Parrinello-Rahman';
    
end

%% Set52: 2000 atoms default settings
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 52;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 10; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 1; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
    Settings_array(idx).MDP.Fourier_Spacing = 0.12;
end

%% Set53: 2000 atoms default settings, except removal of dispersion correction
% Dispersion correction removed
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 53;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = false; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 10; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 1; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
    Settings_array(idx).MDP.Fourier_Spacing = 0.12;
end

%% Set54: 2000 atoms default settings, except decreasing Ewald_rtol and Fourier_Spacing
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 54;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.0;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 10; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 1; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-7; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
    Settings_array(idx).MDP.Fourier_Spacing = 0.10;
end

%% Set55: Testing Pre-Equilibration + Initial equilibration with my "default" parameters + dispersion correct added
% Dispersion correction not added
% Cutoff 1.9 nm
% Now the compressibility is the same in all spatial directions
% Additional 1ps relaxation added to start of simulation
Set = 55;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290; % K
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+Set)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).N_Calc = 2; % Number of chained calculations
    Settings_array(idx).Hours = 6; % Max time for each job (hours)
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    
    Settings_array(idx).MP_Equilibrate_Solid = 10; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 1; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    Settings_array(idx).Barostat = 'Parrinello-Rahman';
    
end

%% Set56: Testing Pre-Equilibration + Initial equilibration with my "default" parameters + Rc = 1.4 + dispersion correction off
% Dispersion correction not added
% Cutoff 1.9 nm
% Now the compressibility is the same in all spatial directions
% Additional 1ps relaxation added to start of simulation
Set = 56;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290; % K
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+Set)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).N_Calc = 2; % Number of chained calculations
    Settings_array(idx).Hours = 6; % Max time for each job (hours)
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    Settings_array(idx).MDP.Disp_Correction = false; % Adds in long-range dispersion correction
    
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    
    Settings_array(idx).MP_Equilibrate_Solid = 10; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 1; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    Settings_array(idx).Barostat = 'Parrinello-Rahman';
    
end

%% Set57: Testing Pre-Equilibration + Initial equilibration with my "default" parameters + Rc = 1.4 + dispersion correction ON
% Dispersion correction not added
% Cutoff 1.9 nm
% Now the compressibility is the same in all spatial directions
% Additional 1ps relaxation added to start of simulation
Set = 57;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290; % K
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+Set)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).N_Calc = 2; % Number of chained calculations
    Settings_array(idx).Hours = 6; % Max time for each job (hours)
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    
    Settings_array(idx).MP_Equilibrate_Solid = 10; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 1; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    Settings_array(idx).Barostat = 'Parrinello-Rahman';
    
end

%% Set58: Testing Pre-Equilibration + Initial equilibration with my "default" parameters + decreasing Ewald_rtol and Fourier_Spacing
% Dispersion correction not added
% Cutoff 1.9 nm
% Now the compressibility is the same in all spatial directions
% Additional 1ps relaxation added to start of simulation
Set = 58;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290; % K
XY_Size = 6; % nm
Z_Size = 10; % nm
Reps = 1:10;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+Set)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).N_Calc = 3; % Number of chained calculations
    Settings_array(idx).Hours = 6; % Max time for each job (hours)
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    Settings_array(idx).MDP.Disp_Correction = false; % Adds in long-range dispersion correction
    
    Settings_array(idx).MP_Equilibrate_Solid = 10; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 1; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    Settings_array(idx).Barostat = 'Parrinello-Rahman';
    
    Settings_array(idx).MDP.Ewald_rtol = 1e-7; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
    Settings_array(idx).MDP.Fourier_Spacing = 0.10;
end

%% Set59: Testing Pre-Equilibration + Initial equilibration with my "default" parameters + Rc = 1.4 + dispersion correction ON + 20k atoms
% Dispersion correction not added
% Cutoff 1.9 nm
% Now the compressibility is the same in all spatial directions
% Additional 1ps relaxation added to start of simulation
Set = 59;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290; % K
Reps = 1:20;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx)+Set)
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).N_Calc = 2; % Number of chained calculations
    Settings_array(idx).Hours = 6; % Max time for each job (hours)
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 20000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    
    Settings_array(idx).MP_Equilibrate_Solid = 10; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 1; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
end

%% Set60: 2000 atoms default settings (full debugged 2022-04-29) - reproduciblility test x 100
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 60;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:100;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.00;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 2000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 5; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    Settings_array(idx).InitialMeshSize = 10;
    Settings_array(idx).MeshSizeMultiplier = 1;
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
    Settings_array(idx).MDP.Fourier_Spacing = 0.12;
end

%% Set61: 11000 atoms default settings (full debugged 2022-04-29) - reproduciblility test x 100
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 61;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:100;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Hours = 6; % Max time for each job (hours)
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.05;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 11000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 5; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    Settings_array(idx).InitialMeshSize = 10;
    Settings_array(idx).MeshSizeMultiplier = 1;
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
    Settings_array(idx).MDP.Fourier_Spacing = 0.12;
end

%% Set63: 11000 atoms default settings (full debugged 2022-04-29) - reproduciblility test with MeltFreezeThreshold = 15%
% Dispersion correction added
% Cutoff 1.4 nm
% Adding new Pre-equilibration strategy
% Compressibility in all dimensions is now fixed at the solid experimental isothermal compressibility
Set = 63;
Salt = 'NaCl';
Theory = 'JC';
T0 = 1290;
Reps = 1:50;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Set + Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Hours = 6; % Max time for each job (hours)
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set' num2str(Set) '_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
    Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
    Settings_array(idx).Cutoff_Buffer = 1.05;
    Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 11000;
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.15;
    
    Settings_array(idx).MaxWarn = 2;
    Settings_array(idx).MP_Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).MP_Equilibrate_Liquid = 5; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    Settings_array(idx).InitialMeshSize = 10;
    Settings_array(idx).MeshSizeMultiplier = 1;
    
    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;
    
    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    
    Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
    Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
    Settings_array(idx).MDP.Fourier_Spacing = 0.12;
end


%% Prod1: Melting points of all alkali halides with "default" parameters + Rc = 1.4 + dispersion correction ON
% Structure: rocksalt
% N = 20000 atoms
% Cutoff 1.4 nm with dispersion correction added
% The compressibility is the same in all spatial directions
% Additional 1ps relaxation added to start of simulation
Experiment = Load_Experimental_Data;
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
Structure = 'Rocksalt';
Theories = {'TF' 'JC' 'JC3P' 'JC4P' 'JCSD'};
Set = 1;

for jdx = 1:length(Salts)
    Salt = Salts{jdx};
    
    switch Salt
        case {'LiF' 'LiCl' 'LiBr' 'LiI'}
            Structures = {'Rocksalt' 'Wurtzite'};
        case {'CsCl' 'CsBr' 'CsI'}
            Structures = {'Rocksalt' 'CsCl'};
        otherwise
            Structures = {'Rocksalt'};
    end
    
    for kdx = 1:length(Theories)
        Theory = Theories{kdx};
        
        if strcmpi(Theory,'JCSD') && ~strcmpi(Salt,'NaCl')
            continue
        end
        
        switch Theory
            case 'TF'
                Theory_Settings = Shared_Settings;
                Theory_Settings.MPI_Ranks = 32; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
                Theory_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
                Theory_Settings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
                Theory_Settings.SinglePrecision = false; % choose true for single precision mode, false for double
                Theory_Settings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
                Theory_Settings.npme = 2; % Number of rank assigned to PME
                Theory_Settings.dd = [2 3 5]; % Domain decomposition
                
                Theory_Settings.MDP.RVDW_Cutoff = 1.9; % nm
                Theory_Settings.MDP.RCoulomb_Cutoff = 1.9; % nm
                Theory_Settings.MDP.RList_Cutoff = 1.9; % nm
                Theory_Settings.Cutoff_Buffer = 1.05; % Not using verlet-buffer-tolerance
                
            otherwise
                Theory_Settings = Shared_Settings;
                Theory_Settings.MPI_Ranks = 4; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
                Theory_Settings.OMP_Threads = 8; % Set the number of OMP threads per MPI rank
                Theory_Settings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
                Theory_Settings.SinglePrecision = false; % choose true for single precision mode, false for double
                Theory_Settings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
                Theory_Settings.npme = 2; % Number of rank assigned to PME
                Theory_Settings.dd = [1 1 2]; %[5 3 2] Domain decomposition
                
                Theory_Settings.MDP.Disp_Correction = true; % Adds in long-range dispersion correction
                Theory_Settings.MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction
                
                Theory_Settings.MDP.RVDW_Cutoff = 1.4; % nm
                Theory_Settings.MDP.RCoulomb_Cutoff = 1.4; % nm
                Theory_Settings.MDP.RList_Cutoff = 1.4; % nm
        end
        
        for mdx = 1:length(Structures)
            Structure = Structures{mdx};
            
            idx = idx+1;
            Settings_array(idx) = Theory_Settings;
            Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            [T0_i,Matchfound] = BestKnownMP(Settings_array(idx));
            if Matchfound
                Settings_array(idx).InitialMeshSize = 5; % K
            else
                Settings_array(idx).InitialMeshSize = 100; % K
            end
            Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
            Settings_array(idx).JobID = ['Prod' num2str(Set)]; % An ID that is tacked onto the folder name of all current jobs
            Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
            Settings_array(idx).T0 = T0_i; % K, Initial temperature
            if contained_in_cell(Structure,{'Rocksalt' 'Sphalerite' 'CsCl'})
                Settings_array(idx).Isotropy = 'semiisotropic';
                Settings_array(idx).Target_P = [1 1]; % Bar
            else
                Settings_array(idx).Isotropy = 'anisotropic';
                Settings_array(idx).Target_P = [1 1 1 1 1 1]; % Bar
            end
            Settings_array(idx).c_over_a = 2;
            Settings_array(idx).N_atoms = 20000;

            Settings_array(idx).BracketThreshold = 1; % K
            Settings_array(idx).MinStepSize = 0.25;
            Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs

            Settings_array(idx).MP_Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
            Settings_array(idx).MP_Equilibrate_Liquid = 5; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
            Settings_array(idx).PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.

            % Barostat Options
            Settings_array(idx).Isotropy = 'semiisotropic';
            Settings_array(idx).Target_P = [1 1]; % Bar
            Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
            Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
            Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
            Settings_array(idx).ScaleCompressibility = 1;

            % Thermostat Options
            Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
            Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
            Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
            Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
            Settings_array(idx).T0 = T0_i; % K, Initial temperature
        end
    end
end

%% Prod2: Melting points of all alkali halides with "default" parameters + Rc = 1.4 + dispersion correction ON
% Structure: rocksalt
% N = 11000 atoms
% Cutoff 1.4 nm with dispersion correction added
% The compressibility is the same in all spatial directions
% Additional 1ps relaxation added to start of simulation
Experiment = Load_Experimental_Data;
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
Theories = {'JC' 'JC3P' 'JC4P' 'JCSD'}; % 'JC' 'JC3P' 'JC4P' 'JCSD'
Set = 2;

for jdx = 1:length(Salts)
    Salt = Salts{jdx};
    
    switch Salt
        case {'LiF' 'LiCl' 'LiBr' 'LiI'}
            Structures = {'Rocksalt' 'Wurtzite'};
        case {'CsCl' 'CsBr' 'CsI'}
            Structures = {'Rocksalt' 'CsCl'};
        otherwise
            Structures = {'Rocksalt'};
    end
    
    for kdx = 1:length(Theories)
        Theory = Theories{kdx};
        
        if strcmpi(Theory,'JCSD') && ~strcmpi(Salt,'NaCl')
            continue
        end
        
        switch Theory
            case 'TF'
                Theory_Settings = Shared_Settings;
                Theory_Settings.MPI_Ranks = 32; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
                Theory_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
                Theory_Settings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
                Theory_Settings.SinglePrecision = false; % choose true for single precision mode, false for double
                Theory_Settings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
                Theory_Settings.npme = 2; % Number of rank assigned to PME
                Theory_Settings.dd = [2 3 5]; % Domain decomposition
                
                Theory_Settings.MDP.RVDW_Cutoff = 1.9; % nm
                Theory_Settings.MDP.RCoulomb_Cutoff = 1.9; % nm
                Theory_Settings.MDP.RList_Cutoff = 1.9; % nm
                Theory_Settings.Cutoff_Buffer = 1.05; % Not using verlet-buffer-tolerance
                
            otherwise
                Theory_Settings = Shared_Settings;
                Theory_Settings.MPI_Ranks = 4; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
                Theory_Settings.OMP_Threads = 8; % Set the number of OMP threads per MPI rank
                Theory_Settings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
                Theory_Settings.SinglePrecision = false; % choose true for single precision mode, false for double
                Theory_Settings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
                Theory_Settings.npme = 2; % Number of rank assigned to PME
                Theory_Settings.dd = [1 1 2]; %[5 3 2] Domain decomposition
                
                Theory_Settings.MDP.Disp_Correction = true; % Adds in long-range dispersion correction
                Theory_Settings.MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction
                
                Theory_Settings.MDP.RVDW_Cutoff = 1.4; % nm
                Theory_Settings.MDP.RCoulomb_Cutoff = 1.4; % nm
                Theory_Settings.MDP.RList_Cutoff = 1.4; % nm
        end
        
        for mdx = 1:length(Structures)
            Structure = Structures{mdx};
            
            idx = idx+1;
            Settings_array(idx) = Theory_Settings;
            Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            [T0_i,Matchfound] = BestKnownMP(Settings_array(idx));
            if Matchfound
                Settings_array(idx).InitialMeshSize = 5; % K
            else
                Settings_array(idx).InitialMeshSize = 100; % K
            end
            Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
            Settings_array(idx).JobID = ['Prod' num2str(Set)]; % An ID that is tacked onto the folder name of all current jobs
            Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
            Settings_array(idx).T0 = T0_i; % K, Initial temperature
            if contained_in_cell(Structure,{'Rocksalt' 'Sphalerite' 'CsCl'})
                Settings_array(idx).Isotropy = 'semiisotropic';
                Settings_array(idx).Target_P = [1 1]; % Bar
            else
                Settings_array(idx).Isotropy = 'anisotropic';
                Settings_array(idx).Target_P = [1 1 1 1 1 1]; % Bar
            end
            Settings_array(idx).c_over_a = 2;
            Settings_array(idx).N_atoms = 11000;

            Settings_array(idx).BracketThreshold = 1; % K
            Settings_array(idx).MinStepSize = 0.25;
            Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs

            Settings_array(idx).MP_Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
            Settings_array(idx).MP_Equilibrate_Liquid = 5; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
            Settings_array(idx).PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.

            % Barostat Options
            Settings_array(idx).Isotropy = 'semiisotropic';
            Settings_array(idx).Target_P = [1 1]; % Bar
            Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
            Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
            Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
            Settings_array(idx).ScaleCompressibility = 1;

            % Thermostat Options
            Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
            Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
            Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
            Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
            Settings_array(idx).T0 = T0_i; % K, Initial temperature
        end
    end
end

