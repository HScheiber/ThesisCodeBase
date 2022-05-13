%% Thermal Expansion Reference for Solids
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
Structures = {'NiAs' 'AntiNiAs' 'CsCl'}; %   'Rocksalt' 'Wurtzite' 'BetaBeO' 'NiAs'  'Sphalerite' 'FiveFive' 'AntiNiAs' 'CsCl' 
Models = {{'JC' ''} {'TF' ''} {'JC' 'EQ'} {'BH' 'FE'}};
Nums_JC = {'2' '2' '1' '5' ''};
Nums_BH = {'5' '6' '7' '5' ''};
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    Num_JC = Nums_JC{sidx};
    Num_BH = Nums_BH{sidx};
    
    for jdx = 1:length(Structures)
        Structure = Structures{jdx};
        
        for kdx = 1:length(Models)
            
            Theory = Models{kdx}{1};
            Model = Models{kdx}{2};

            if ~isempty(Model) && strcmp(Theory,'JC') && isempty(Num_JC)
                continue
            elseif ~isempty(Model) && strcmp(Theory,'BH') && isempty(Num_BH)
                continue
            elseif ~isempty(Model) && strcmp(Theory,'JC')
                Model = [Model Num_JC];
            elseif ~isempty(Model) && strcmp(Theory,'BH')
                Model = [Model Num_BH];
            end
        
            %% Thermal expansion for solid starting at 0 K and heating to 2200 K at 100 K / ns
            idx = idx+1;
            Settings_array(idx) = Shared_Settings;
            Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
            Settings_array(idx).JobID = 'ExpRef'; % An ID that is tacked onto the folder name of all current jobs
            Settings_array(idx).N_atoms = 5000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
            Settings_array(idx).Target_T = 0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            Settings_array(idx).MDP.Initial_T = 0; % Initial termpature at which to generate velocities
            Settings_array(idx).T0 = 0; % K, Initial temperature
            Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
            Settings_array(idx).Annealing_Times = [0  22000]; % [ps] A list with the number of annealing reference/control points used
            Settings_array(idx).Annealing_Temps = [0  2200]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
            Settings_array(idx).MDP.Trajectory_Time = 22; % ns
            if strcmp(Salt,'LiF')
                Settings_array(idx).JobLinks = Settings_array(idx).JobLinks + 4;
            end
        end
    end
end

%% Thermal Expansion Reference CsCl 
Salts = {'CsCl'};
Structures = {'CsCl'}; %   'Rocksalt' 'Wurtzite' 'BetaBeO' 'NiAs'  'Sphalerite' 'FiveFive' 'AntiNiAs' 'CsCl' 
Models = {'JC' 'TF'};
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for jdx = 1:length(Structures)
        Structure = Structures{jdx};
        
        for kdx = 1:length(Models)
            
            Theory = Models{kdx};
        
            %% Thermal expansion for solid starting at 0 K and heating to 2200 K at 100 K / ns
            idx = idx+1;
            Settings_array(idx) = Shared_Settings;
            Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
            Settings_array(idx).JobID = 'ExpRef'; % An ID that is tacked onto the folder name of all current jobs
            Settings_array(idx).N_atoms = 5000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
            Settings_array(idx).Target_T = 0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            Settings_array(idx).MDP.Initial_T = 0; % Initial termpature at which to generate velocities
            Settings_array(idx).T0 = 0; % K, Initial temperature
            Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
            Settings_array(idx).Annealing_Times = [0  22000]; % [ps] A list with the number of annealing reference/control points used
            Settings_array(idx).Annealing_Temps = [0  2200]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
            Settings_array(idx).MDP.Trajectory_Time = 22; % ns
            if strcmp(Salt,'LiF')
                Settings_array(idx).JobLinks = Settings_array(idx).JobLinks + 4;
            end
        end
    end
end

%% Default JC Models at high temperature: Aiming for vapour
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
T0 = 3000; % Initial temperature
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
        
    %% Generating gas (JC)
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = 'JC'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Liquid'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = 'GenGas'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 5000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
    Settings_array(idx).Annealing_Times = [0    500   20500]; % [ps] A list with the number of annealing reference/control points used
    Settings_array(idx).Annealing_Temps = [3000 3000  1000]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
    Settings_array(idx).MDP.Trajectory_Time = 20.5; % ns
    if strcmp(Salt,'LiF')
        Settings_array(idx).JobLinks = Settings_array(idx).JobLinks + 4;
    end
    
    %% Generating amorphous solid/glass with JC model
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = 'JC'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Liquid'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = 'GenGlass'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 5000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = 2000; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = 2000; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = 2000; % K, Initial temperature
    Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
    Settings_array(idx).Annealing_Times = [0    999.999 1000  31000]; % [ps] A list with the number of annealing reference/control points used
    Settings_array(idx).Annealing_Temps = [2000 2000    0     2000]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
    Settings_array(idx).MDP.Trajectory_Time = 31; % ns
    Settings_array(idx).Cutoff_Buffer = 1.25;
    if strcmp(Salt,'LiF')
        Settings_array(idx).JobLinks = Settings_array(idx).JobLinks + 4;
    end
    
    %% Generating amorphous solid/glass with TF model
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = 'TF'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Liquid'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = 'GenGlass'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).N_atoms = 5000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings_array(idx).Target_T = 2000; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = 2000; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = 2000; % K, Initial temperature
    Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
    Settings_array(idx).Annealing_Times = [0    999.999 1000  31000]; % [ps] A list with the number of annealing reference/control points used
    Settings_array(idx).Annealing_Temps = [2000 2000    0     2000]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
    Settings_array(idx).MDP.Trajectory_Time = 31; % ns
    Settings_array(idx).Cutoff_Buffer = 1.25;
    if strcmp(Salt,'LiF')
        Settings_array(idx).JobLinks = Settings_array(idx).JobLinks + 4;
    end
end

%% Thermal Expansion Reference for liquids
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
Models = {{'JC' ''} {'TF' ''} {'JC' 'EQ'} {'BH' 'FE'}};
Nums_JC = {'2' '2' '1' '5' ''};
Nums_BH = {'5' '6' '7' '5' ''};

T0 = 2200; % Initial temperature
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    Num_JC = Nums_JC{sidx};
    Num_BH = Nums_BH{sidx};
    
    for jdx = 1:length(Models)
        Theory = Models{jdx}{1};
        Model = Models{jdx}{2};
        
        if ~isempty(Model) && strcmp(Theory,'JC') && isempty(Num_JC)
            continue
        elseif ~isempty(Model) && strcmp(Theory,'BH') && isempty(Num_BH)
            continue
        elseif ~isempty(Model) && strcmp(Theory,'JC')
            Model = [Model Num_JC];
        elseif ~isempty(Model) && strcmp(Theory,'BH')
            Model = [Model Num_BH];
        end
    
        %% Thermal expansion for molten salt starting at 2200 K and cooling to 300 K at 100 K / ns
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Liquid'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = 'ExpRef'; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).N_atoms = 5000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
        Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = T0; % K, Initial temperature
        Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
        Settings_array(idx).Annealing_Times = [0    500   19500]; % [ps] A list with the number of annealing reference/control points used
        Settings_array(idx).Annealing_Temps = [T0   T0  300]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
        Settings_array(idx).MDP.Trajectory_Time = 19.5; % ns
        Settings_array(idx).Cutoff_Buffer = 1.25;
        if strcmp(Salt,'LiF')
            Settings_array(idx).JobLinks = Settings_array(idx).JobLinks + 4;
        end
        if strcmp(Salt,'LiF') && strcmp(Theory,'BH')
            Settings_array(idx).Cutoff_Buffer = 1.5;
        end
    end
end

%% Thermal Expansion Reference for Solids
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
Structures = {'Rocksalt' 'Wurtzite' 'BetaBeO' 'Sphalerite' 'FiveFive'}; %   'Rocksalt' 'Wurtzite' 'BetaBeO' 'NiAs'  'Sphalerite' 'FiveFive' 'AntiNiAs' 'CsCl' 
Models = {{'JC' ''} {'TF' ''} {'JC' 'EQ'} {'BH' 'FE'}};
Nums_JC = {'2' '2' '1' '5' ''};
Nums_BH = {'5' '6' '7' '5' ''};
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    Num_JC = Nums_JC{sidx};
    Num_BH = Nums_BH{sidx};
    
    for jdx = 1:length(Structures)
        Structure = Structures{jdx};
        
        for kdx = 1:length(Models)
            
            Theory = Models{kdx}{1};
            Model = Models{kdx}{2};

            if ~isempty(Model) && strcmp(Theory,'JC') && isempty(Num_JC)
                continue
            elseif ~isempty(Model) && strcmp(Theory,'BH') && isempty(Num_BH)
                continue
            elseif ~isempty(Model) && strcmp(Theory,'JC')
                Model = [Model Num_JC];
            elseif ~isempty(Model) && strcmp(Theory,'BH')
                Model = [Model Num_BH];
            end
        
            %% Thermal expansion for solid starting at 0 K and heating to 2200 K at 100 K / ns
            idx = idx+1;
            Settings_array(idx) = Shared_Settings;
            Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
            Settings_array(idx).JobID = 'ExpRef'; % An ID that is tacked onto the folder name of all current jobs
            Settings_array(idx).N_atoms = 5000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
            Settings_array(idx).Target_T = 0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            Settings_array(idx).MDP.Initial_T = 0; % Initial termpature at which to generate velocities
            Settings_array(idx).T0 = 0; % K, Initial temperature
            Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
            Settings_array(idx).Annealing_Times = [0  22000]; % [ps] A list with the number of annealing reference/control points used
            Settings_array(idx).Annealing_Temps = [0  2200]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
            Settings_array(idx).MDP.Trajectory_Time = 22; % ns
            if strcmp(Salt,'LiF')
                Settings_array(idx).JobLinks = Settings_array(idx).JobLinks + 4;
            end
        end
    end
end
