%%%%% Melting_Point_Batch_Submit_Cedar %%%%%%

%% Global Calculation settings
clear;

% Set up calculations below
skip_calculations = [];
check_complete = true; % Checks if job is already completed, skips completed jobs
check_running = true; % Checks if a job is already running, skips running jobs

% Load shared resource and mdrun settings
Shared_Settings = Initialize_MD_Settings;
Shared_Settings.Project_Directory_Name = 'Model_Building';
Shared_Settings.BatchMode = false; % Sets up batch job when true, or runs immediately when false
Shared_Settings.Submit_Jobs = false; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Shared_Settings.JobSettings.N_Calc = 1; % Number of chained calculations
Shared_Settings.JobSettings.Hours = 6; % Max time for each job (hours)
Shared_Settings.JobSettings.Mins = 0; % Max time for job (minutes)
Shared_Settings.JobSettings.Nodes = 0; % Minimum number of cores to request for calculation.
Shared_Settings.JobSettings.Cores = 12; % Minimum number of cores to request for calculation. Set to -1 for entire node
Shared_Settings.JobSettings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Shared_Settings.JobSettings.SinglePrecision = false; % choose true for single precision mode, false for double
Shared_Settings.JobSettings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
Shared_Settings.JobSettings.MPI_Ranks = 2; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.JobSettings.OMP_Threads = 6; % Set the number of OMP threads per MPI rank
Shared_Settings.JobSettings.npme = 0; % Number of rank assigned to PME
Shared_Settings.JobSettings.dd = [1 1 2]; % Domain decomposition


% Shared calculation parameters
Shared_Settings.Liquid_Interface = true; % When true, creates an system with half STRUCTURE half LIQUID for melting point testing
Shared_Settings.MeltFreezeThreshold = 0.25; % CHANGE in fraction [0,1] OR Number of atoms (1,inf) of liquid/solid required to establish a phase change
Shared_Settings.Optimizer = 'MPSearcher';
Shared_Settings.lb = 0; % K, lower bound on MP search
Shared_Settings.ub = 2200; % K, upper bound on MP search
Shared_Settings.BracketThreshold = 5; % [K] Sets the target bracket for the melting point
Shared_Settings.MinStepSize = 0.25; % [K] Sets the minimum step size for MPsearcher algorithm
Shared_Settings.SlopeThreshold = 1e10; % The change in the % fraction per unit time must be smaller than the absolute value of this threshold for the system to be considered at the melting point. Units of [% Structure Fraction/ps]
Shared_Settings.Liquid_Fraction = 0.50;
Shared_Settings.MaxTDiff = 0.01; % K, maximum change in temperature between points before selecting new initial conditions
Shared_Settings.MaxWarn = 2;
Shared_Settings.Liquid_Test_Time = 50; % ps. simulation time to sample the liquid (second half averaged for enthalpy / volume)
Shared_Settings.Solid_Test_Time = 30; % ps. simulation time to sample the solid (second half averaged for enthalpy / volume)


Exp = Load_Experimental_Data;

% Initial calculation index
idx=0;

%% JC Models EB, EC, ...
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theory = 'JC';
Models = {'EB' 'EC' 'EF' 'EH' 'EI' 'EK' 'EL' 'EN' 'EO' 'EQ' 'ER' 'EU'};
Reps = 1:10;
for jdx = 1:length(Salts)
    Salt = Salts{jdx};
    for kdx = 1:length(Models)
        for ldx = 1:length(Reps)
            Model = [Models{kdx} num2str(Reps(ldx))];

            idx = idx+1;
            Settings_array(idx) = Shared_Settings;
            Settings_array(idx).JobSettings.N_Calc = 1;
            Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
            Settings_array(idx).JobID = [Theory '_Model_' Model]; % An ID that is tacked onto the folder name of all current jobs
            Settings_array(idx).Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
            Settings_array(idx).MDP.RVDW_Cutoff = 1.00; % nm
            Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
            Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
            Settings_array(idx).Cutoff_Buffer = 1.20;
            Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
            Settings_array(idx).c_over_a = 2;
            Settings_array(idx).N_atoms = 2000;
            Settings_array(idx).BracketThreshold = 5; % K
            Settings_array(idx).MinStepSize = 0.25;
            Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
            Settings_array(idx).MeltFreezeThreshold = 0.25;

            Settings_array(idx).MaxWarn = 2;
            Settings_array(idx).Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
            Settings_array(idx).Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
            Settings_array(idx).PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
            Settings_array(idx).InitialMeshSize = 20;
            Settings_array(idx).MeshSizeMultiplier = 5;
            Settings_array(idx).QECompressibility = 1e-7; % sets the compressibility during the system preparation stages
            Settings_array(idx).ScaleInitialLiqDensity = 0.8;
            Settings_array(idx).Delete_Equil = false; % switch to delete temporary calculation folders for finite T calcs
            
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
            Settings_array(idx).Target_T = Exp.(Salt).mp; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            Settings_array(idx).MDP.Initial_T = Exp.(Salt).mp; % Initial termpature at which to generate velocities
            Settings_array(idx).T0 = Exp.(Salt).mp; % K, Initial temperature

            Settings_array(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
            Settings_array(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
            Settings_array(idx).MDP.Fourier_Spacing = 0.12;
            Settings_array(idx).MDP.VerletBT = -1;
        end
    end
end


%% Check for already running jobs
if check_running && isunix
    er = 1;
    jdx = 1;
    while er ~= 0
        [er,out] = system('squeue -u $USER -o %100j | tail -n +2');
        jdx = jdx+1;
        if jdx > 10
            disp('Failed to check current running jobs.')
            out = '';
            break
        end
    end
    
    % Grab the currently running jobs
    out = strtrim(out);
    c_jobs = strtrim(split(out,newline));
    
    % Remove the job number tags
    c_jobs_unique = unique(regexprep(c_jobs,'-([0-9]){3}',''));
else
    c_jobs_unique = {''};
end

%% Prepare batch scripts
[~,~,~,Slurm] = find_home;
calc_cmd = 'matlab -nodisplay -r "Thermal_Properties_At_MP(##INPUTFILE##)" >> ##LOGFILE##';

% Loop through all models and build their submission file
for idx = 1:length(Settings_array)
    if any(idx == skip_calculations)
        continue
    end
    Settings = Settings_array(idx);
    
    % Load the batch script template and remove unnecessary fields
    Batch_Template = MD_Batch_Template(Settings.JobSettings);
    Batch_Template = strrep(Batch_Template,['##PREMIN##' newline],'');
    Batch_Template = strrep(Batch_Template,['##CLEANUP##' newline],'');
    Batch_Template = strrep(Batch_Template,['##EXT1##' newline],'');
    Batch_Template = strrep(Batch_Template,['##EXT2##' newline],'');
    
    [Settings.WorkDir,Settings.JobName,Settings.Full_Model_Name] = GetMDWorkdir(Settings,Settings.JobID);
    OuterDir = Settings.WorkDir;
    if ~exist(Settings.WorkDir, 'dir')
    	% Model not found
    	continue
    end
    Settings.WorkDir = fullfile(OuterDir,'BestPoint_Thermal');
    Batch_Template = strrep(Batch_Template,'##DIRECTORY##',Settings.WorkDir);
    TaskName = [Settings.Salt '_' Settings.JobName];
    
    % Generate and move to the submission directory
    if ~exist(Settings.WorkDir, 'dir')
       mkdir(Settings.WorkDir)
    end
    cd(Settings.WorkDir);
    
    % Check if job is already complete
    FlagCompleteFile = fullfile(Settings.WorkDir,'Thermal_Properties_Complete');
    if check_complete && isfile(FlagCompleteFile)
        disp([TaskName ': Job already completed. Skipping Job Submission.'])
        continue
    end
    
    % Check if Model is already running
    if check_running && ismember(TaskName,c_jobs_unique) 
        disp([TaskName ': Job already Running. Skipping Job Submission.'])
        continue
    end
    
    % Save input file
    inp_file = fullfile(Settings.WorkDir,[Settings.JobName '.inp']);
    save(inp_file,'Settings','-mat')
    
    % Create input command
    calc_cmd_idx = strrep(calc_cmd,'##INPUTFILE##',['''' Settings.JobName '.inp' '''']);
    calc_cmd_idx = strrep(calc_cmd_idx,'##LOGFILE##',[Settings.JobName '.MPlog']);
    
    % Set up job links
    for jdx = 1:Settings.JobSettings.N_Calc
        
        if jdx == 1
            calc_cmd_idx_jdx = calc_cmd_idx;
        else
            EXT1 = ['if [[ ! -f "' FlagCompleteFile '" ]]; then'];
            EXT2 = 'fi';
            calc_cmd_idx_jdx = [EXT1 newline calc_cmd_idx newline EXT2];
        end
        Batch_Template_jdx = strrep(Batch_Template,'##MDRUN##',calc_cmd_idx_jdx);
        
        % Copy over batch template and modify
        TaskName_jdx = [TaskName '-' num2str(jdx,'%03.f')];
        Batch_Template_jdx = strrep(Batch_Template_jdx,'##TASKNAME##',TaskName_jdx);
        Batch_Template_jdx = strrep(Batch_Template_jdx,'##ERROR##',TaskName_jdx);

        % Save job submission script
        subm_file_jdx = fullfile(Settings.WorkDir,[Settings.JobName '-' num2str(jdx,'%03.f') '.subm']);
        fid = fopen(subm_file_jdx,'wt');
        fwrite(fid,regexprep(Batch_Template_jdx,{'\r'}',{''}));
        fclose(fid);

        % Submit job
        if jdx == 1
            if Slurm
                qsub_cmd = ['sbatch ' subm_file_jdx];
            else
                qsub_cmd = ['qsub ' subm_file_jdx];
            end
        else
            if Slurm
                qsub_cmd = ['sbatch --depend=afterany:' PrevJobID ' ' subm_file_jdx ];
            else
                qsub_cmd = ['qsub -W depend=afterany:' PrevJobID ' ' subm_file_jdx ];
            end
        end
        
        disp(['Submitting Job ' TaskName_jdx])
        [errcode,output] = system(qsub_cmd);


        if errcode ~= 0 && ~Slurm
            disp(output);
            error(['Error submitting Job: ' newline qsub_cmd]);
        end

        % Attempt to recover from failed job submissions on SLURM servers
        while errcode ~= 0
            disp(output);
            disp(['Error submitting Job: ' newline qsub_cmd]);
            disp('Retrying.');

            errcode2 = 1;
            while errcode2 ~= 0
                pause(5);

                % Check if the job was actually submit or not
                check_cmd = 'squeue -u $USER -o "%10i %60j"';
                [errcode2,output] = system(check_cmd);
            end

            % Check if the previous job is in the queue
            regcheck = ['([0-9]+) +' TaskName_jdx];
            jobid = regexp(output,regcheck,'once','tokens');

            if isempty(jobid) % Job did not submit successfully
                % Retry submitting
                disp(['Submitting Job ' TaskName_jdx])
                [errcode,output] = system(qsub_cmd);

            else % Job did submit successfully (but failed to output properly)
                output = jobid{1};
                errcode = 0;
            end
        end
        
        PrevJobID = regexp(output,'[0-9]+','match','ONCE');
        pause(0.5);
    end
end