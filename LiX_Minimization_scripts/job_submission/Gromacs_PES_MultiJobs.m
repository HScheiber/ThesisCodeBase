% Gromacs_PES_MultiJobs
JobType = 'Structure_Minimization'; % 'Structure_Minimization'
Optimize_Position = true; % Only meaningful for structure minimization: activates optimization of fractional coordinates
Salts = {'LiF'};
Structures = {'Wurtzite'}; %{'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'};
Models = {'TF'}; % Input model(s) to use: JC, JC3P, JC4P, TF
Scale_Dispersion = 1.0; % Works for both JC and TF
Scale_Epsilon = 1.0; % Works only for JC
Parallel = false;
Hours = 12;
Minutes = 0;
Nodes = 1;
Submit = true; % Set to false for debugging

% Activate GLOST?
Glost_On = true; % activate GLOST: Greedy Launcher Of Small Tasks (does not work with parallel matlab!)
Glost_Nodes = 1; % How many nodes to request for GLOST job

% Cannot run GLOST and parallel simultaneously
if Parallel && Glost_On
    error('Cannot run in Parallel with GLOST')
end

% Make sure the Glost_txt variable does not already exist
clearvars Glost_txt

[~,Servertxt] = system('hostname -s | cut -c 1-3');
Server = strtrim(Servertxt);

if Parallel || Glost_On
    Cores_Per_Node = feature('numcores');
else
    Cores_Per_Node = 1;
end

if ispc 
    Username = 'Trusty';
    Account = '';
    Submission_dir = 'C:\Users\Hayden\Documents\Patey_Lab\GTesting';
elseif strcmpi(Server,'ced') || strcmp(Server,'cdr')
    Username = 'scheiber';
    Account = 'rrg-patey-ad';
    Submission_dir = '/home/scheiber/project/GROMACS_LiHalides';
elseif strcmpi(Server,'gra')
    Username = 'scheiber';
    Account = 'def-patey';
    Submission_dir = '/home/scheiber/project/GROMACS_LiHalides';
elseif strcmpi(Server,'sea') || strcmpi(Server,'pod')
    Username = 'scheiber';
    Account = 'def-patey';
    Submission_dir = '/home/scheiber/project/GROMACS_LiHalides';
elseif strcmpi(Server,'pat')
    Username = 'user';
    Account = '';
    Submission_dir = '/home/user/project/GROMACS_LiHalides';
else
    Username = 'scheiber';
    Account = 'def-patey';
    Submission_dir = '/home/scheiber/project/GROMACS_LiHalides';
    if Glost_On
        warning(['GLOST may not work on this server: ' Server])
    end
end

if Glost_On
    % Generate automatic glost job name
    if Optimize_Position
        poslab = '_Full';
    else
        poslab = '_Cell';
    end
    Labels = strrep(Label_replace(Structures),'_','');
    disper = [num2str(min(Scale_Dispersion),'%2.1f') '-' num2str(max(Scale_Dispersion),'%2.1f')];
    Glost_Name = [Salts{:} '_' Labels{:} '_' Models{:} '_D' disper poslab];
    
    % Create GLOST Jobs folder
    Glost_Directory = [Submission_dir filesep 'GLOST_INPUTS'];
    if ~exist(Glost_Directory,'dir')
        mkdir(Glost_Directory)
    end
    
    % Make submission script
    if strcmpi(Server,'pat')
        cd(Glost_Directory);
        Glost_submission = ['mpiexec glost_launch ' Glost_Name '.joblist'];
    else
        Glost_submission = [ ...
            '#!/bin/bash' newline ...
            '#SBATCH --time=' num2str(Hours) ':0:00' newline ...
            '#SBATCH --nodes=' num2str(Glost_Nodes) newline ...
            '#SBATCH --tasks-per-node=' num2str(Cores_Per_Node) newline ...
            '#SBATCH --cpus-per-task=1' newline ...
            '#SBATCH --mem-per-cpu=MaxMemPerCPU' newline ...
            '#SBATCH --account=' Account newline ...
            '#SBATCH --job-name=' Glost_Name newline ...
            '#SBATCH --error=' Glost_Name '.stde' newline ...
            '#SBATCH --export=ALL' newline ...
            '#SBATCH --exclusive' newline newline ...
            '# set EXE environment' newline ...
            'PATH=$PATH' newline ...
            'export PATH' newline newline ...
            '# Load GLOST module along with the modules required' newline ...
            'module load nixpkgs/16.09  intel/2016.4  openmpi/2.1.1 glost/0.3.1' newline newline ...
            '# Load additional program modules' newline ... 
            'module load matlab/2018a' newline ... 
            'module load gromacs/2018' newline ... 
            'echo "Starting run at: `date`"' newline newline ...
            '# Move to GLOST start directory' newline ...
            'cd ' Glost_Directory newline newline ...        
            '# Run GLOST with the job list' newline ...
            'mpiexec glost_launch ' Glost_Name '.joblist' newline ....
            'grep -rnL ' Glost_Name ' -e "convergence reached after" >> ' Glost_Name filesep 'IncompleteJobs.txt' newline ...
            'echo "Job finished at"' newline ...
            'date' newline ...
            '################### Job Ended ###################' newline ...
            'exit 0'];
        qsub = 'sbatch';
    end
    
    EnvironmentSetupText = ['addpath(''/home/scheiber/.matlab/R2018a''); startup; ' ... 
        JobType '(##SALT##,##STRUCTURE##,##MODEL##,##DISPERSION##,##EPSILON##,##PARALLEL##); ' ... 
        'quit'];
    
elseif strcmpi(Server,'sea') || strcmpi(Server,'pod')  % Orcinus
    TemplateText = ['#!/bin/bash' newline ...
        '#PBS -S /ThesisCodeBase/bash' newline ...
        '#PBS -W x=GRES:MATLAB' newline ...
        '#PBS -l walltime=' num2str(Hours) ':' num2str(Minutes) ':00' newline ...
        '#PBS -l nodes=' num2str(Nodes) ':ppn=' num2str(Cores_Per_Node) newline ...
        '#PBS -V' newline ...
        '#PBS -N ##TASKNAME##' newline ...
        '#PBS -e ##TASKNAME##.stde' newline ...
        '#PBS -o ##TASKNAME##.stdo' newline ...
        '' newline ...
        '# Check on some basics:' newline ...
        'echo "Running on host: " `hostname`' newline ...
        'echo "Changing to directory from which PBS script was submitted."' newline ...
        'cd $PBS_O_WORKDIR' newline ...
        'echo "Current working directory is now: " `pwd`' newline ...
        'echo "Starting MATLAB at `date`"' newline ...
        newline ...
        newline ...
        '# set EXE environment' newline ...
        'unalias matlab' newline ...
        'module load matlab/matlab_2014a' newline ...
        'module load gromacs/5.1.4' newline ...
        'matlab -nodisplay -nojvm -nodesktop -nosplash < ##TASKNAME##.m > ##TASKNAME##.log' newline ...
        'echo "MATLAB run completed at `date`"' newline ...
        'exit 0'];

    EnvironmentSetupText = ['addpath(''/home/scheiber/.matlab/R2014a'') ' newline ... 
        'startup; ' newline ... 
        JobType '(##SALT##,##STRUCTURE##,##MODEL##,##DISPERSION##,##EPSILON##,##PARALLEL##);' newline ... 
        'quit'];
    
    qsub = 'qsub';

elseif strcmpi(Server,'ced') || strcmpi(Server,'cdr') || strcmpi(Server,'gra')
    TemplateText = ['#!/bin/bash' newline ... 
        '#SBATCH --time=' num2str(Hours) ':' num2str(Minutes) ':00' newline ... 
        '#SBATCH --nodes=' num2str(Nodes) newline ... 
        '#SBATCH --tasks-per-node=1' newline ... 
        '#SBATCH --cpus-per-task=' num2str(Cores_Per_Node) newline ... 
        '#SBATCH --mem-per-cpu=MaxMemPerCPU' newline ... 
        '#SBATCH --account=' Account newline ... 
        '#SBATCH --job-name=##TASKNAME##' newline ... 
        '#SBATCH --error=##TASKNAME##.stde' newline ... 
        '#SBATCH --export=ALL' newline ... 
        newline ... 
        '# Check on some basics:' newline ... 
        'echo "Running on host: " `hostname`' newline ... 
        'echo "Changing to directory from which PBS script was submitted."' newline ... 
        'cd $SLURM_SUBMIT_DIR' newline ... 
        'echo "Current working directory is now: " `pwd`' newline ... 
        'echo "Starting MATLAB at `date`"' newline ... 
        newline ... 
        newline ... 
        '# set EXE environment' newline ... 
        'module load matlab/2018a' newline ... 
        'module load gromacs/2018' newline ... 
        'matlab -nodisplay -nodesktop -nosplash < ##TASKNAME##.m > ##TASKNAME##.log' newline ... 
        'echo "MATLAB run completed at `date`"' newline ... 
        'exit 0'];
    
    EnvironmentSetupText = ['addpath(''/home/scheiber/.matlab/R2018a'') ' newline ... 
        'startup; ' newline ... 
        JobType '(##SALT##,##STRUCTURE##,##MODEL##,##DISPERSION##,##EPSILON##,##PARALLEL##);' newline ... 
        'quit'];
    
    qsub = 'sbatch';
elseif strcmpi(Server,'pat')
    EnvironmentSetupText = ['source /home/user/Documents/MATLAB/.matlabrc; ' newline ... 
        JobType '(##SALT##,##STRUCTURE##,##MODEL##,##DISPERSION##,##EPSILON##,##PARALLEL##);' newline ... 
        'quit'];
    qsub = '';
end

if strcmp(JobType,'Structure_Minimization')
    if Parallel || Glost_On
        if Optimize_Position
            EnvironmentSetupText = strrep(EnvironmentSetupText,'##PARALLEL##','true,true');
        else
            EnvironmentSetupText = strrep(EnvironmentSetupText,'##PARALLEL##','true,false');
        end
    else
        if Optimize_Position
            EnvironmentSetupText = strrep(EnvironmentSetupText,'##PARALLEL##','false,true');
        else
            EnvironmentSetupText = strrep(EnvironmentSetupText,'##PARALLEL##','false,false');
        end
    end
else % Not structure minimization
    if Parallel
        EnvironmentSetupText = strrep(EnvironmentSetupText,'##PARALLEL##','true,false'); %#ok<*UNRCH>
    elseif Glost_On
        EnvironmentSetupText = strrep(EnvironmentSetupText,'##PARALLEL##','false,true');
    else
        EnvironmentSetupText = strrep(EnvironmentSetupText,'##PARALLEL##','false,false');
    end
end

for i=1:length(Salts)
    Salt = Salts{i};
    
    % Update Directory
    Current_Directory = fullfile(Submission_dir,Salt);
    
    % Create directory if it does not exist
    if ~exist(Current_Directory,'dir')
        mkdir(Current_Directory)
    end
    
    
    for j=1:length(Structures)
        Structure = Structures{j};
        
        % Update Directory
        Current_Directory = fullfile(Submission_dir,Salt,Structure);

        % Create directory if it does not exist
        if ~exist(Current_Directory,'dir')
            mkdir(Current_Directory)
        end
        
        for k=1:length(Models)
            Model = Models{k};
            
            for l=1:length(Scale_Dispersion)
                
                Dispersion_Scale = Scale_Dispersion(l);
                
                if ismembertol(1.0,Dispersion_Scale,1e-5)
                    Model_Scaled = Model;
                else
                    Model_Scaled = [Model '_D' sprintf('%7.5f',Dispersion_Scale)];
                end
                
                for m=1:length(Scale_Epsilon)
                    Epsilon_Scale = Scale_Epsilon(m);
                
                    if ismembertol(1.0,Scale_Epsilon,1e-5)
                        Model_Scaled2 = Model_Scaled;
                    else
                        Model_Scaled2 = [Model_Scaled '_E' sprintf('%7.5f',Epsilon_Scale)];
                    end
                    
                    TaskName = [Salt '_' Structure '_' Model_Scaled2];
                        
                    if Glost_On
                        
                        % Update job text
                        CurrentEnvironmentText = strrep(EnvironmentSetupText,'##SALT##',...
                            ['{''' Salt '''}']);
                        CurrentEnvironmentText = strrep(CurrentEnvironmentText,'##STRUCTURE##',...
                            ['{''' Structure '''}']);
                        CurrentEnvironmentText = strrep(CurrentEnvironmentText,'##MODEL##',...
                            ['{''' Model '''}']);
                        CurrentEnvironmentText = strrep(CurrentEnvironmentText,'##DISPERSION##',...
                            num2str(Dispersion_Scale));
                        CurrentEnvironmentText = strrep(CurrentEnvironmentText,'##EPSILON##',...
                            num2str(Epsilon_Scale));
                        
                        % Make current job script
                        Script_Name = fullfile(Glost_Directory,[TaskName '.m']);
                        Output_Name = fullfile(Glost_Directory,Glost_Name,[TaskName '.glostlog']);
                        
                        % Save matlab script
                        fid = fopen(Script_Name,'wt');
                        fwrite(fid,CurrentEnvironmentText);
                        fclose(fid);                        
                        
                        % Generate job line
                        Glost_newline = ['export PATH=$PATH; matlab -nodisplay -nodesktop -nosplash < '  ...
                            Script_Name ' > ' Output_Name '; rm -f ' Script_Name];

                        % Add new command to command list
                        if ~exist('Glost_txt','var')
                            Glost_txt = Glost_newline;
                        else
                            Glost_txt = [Glost_txt newline Glost_newline];
                        end
                    else
                        %% Glost not on

                        % Update Directory
                        Current_Directory = fullfile(Submission_dir,Salt,Structure,Model_Scaled2);

                        % Create directory if it does not exist
                        if ~exist(Current_Directory,'dir')
                            mkdir(Current_Directory)
                        end


                        % Move to current directory for job submit
                        cd(Current_Directory)

                        % Add in settings
                        CurrentBatchText = strrep(TemplateText,'##TASKNAME##',TaskName);

                        CurrentEnvironmentText = strrep(EnvironmentSetupText,'##SALT##',...
                            ['{''' Salt '''}']);
                        CurrentEnvironmentText = strrep(CurrentEnvironmentText,'##STRUCTURE##',...
                            ['{''' Structure '''}']);
                        CurrentEnvironmentText = strrep(CurrentEnvironmentText,'##MODEL##',...
                            ['{''' Model '''}']);
                        CurrentEnvironmentText = strrep(CurrentEnvironmentText,'##DISPERSION##',...
                            num2str(Dispersion_Scale));
                        CurrentEnvironmentText = strrep(CurrentEnvironmentText,'##EPSILON##',...
                            num2str(Epsilon_Scale));

                        % Print text into text files
                        Subid = fopen([TaskName '.subm'],'wt');
                        fprintf(Subid,CurrentBatchText);
                        fclose(Subid);

                        Envid = fopen([TaskName '.m'],'wt');
                        fprintf(Envid,CurrentEnvironmentText);
                        fclose(Envid);

                        % Submit job
                        if Submit
                            disp(['Submitting Job: ' TaskName])
                            [status,~] = system([qsub ' ' TaskName '.subm']);

                            if status == 0
                                disp('Success, job sumitted.')
                            else
                                error(['Job: "' TaskName '" failed to submit.'])
                            end
                        else
                            disp(['Test mode, job not submitted: ' TaskName])
                        end
                    end
                end
            end
        end
    end
end

% Run Final Glost command
if Glost_On
    
    cd(Glost_Directory)
    if ~exist(Glost_Name,'dir')
        mkdir(Glost_Name)
    end
    
    % Save jobslist file
    FID = fopen([Glost_Name '.joblist'],'wt');
    fprintf(FID,'%s',Glost_txt);
    fclose(FID);
    
    % Save submission script
    FID = fopen([Glost_Name '.subm'],'w');
    fprintf(FID,'%s',Glost_submission);
    fclose(FID);

    % Submit job
    command = [qsub ' ' Glost_Name '.subm'];
    disp(['Submitting New GLOST Job Set: ' Glost_Name '...']);
    if ~ispc && Submit
        system(command);
        disp('Job Sumitted.');
    else
        disp(command); % for debugging
    end
end

cd(Submission_dir)
