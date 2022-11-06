% Run_Project_Structure_Classifier_Batch

clear;
Settings = Initialize_MD_Settings;
% Load Settings object and set some shared calculation settings
Settings.JobLinks = 20; % Number of jobs to link together.
Settings.Hours = 3; % Max time for each job (hours)
Settings.Mins = 0; % Max time for job (minutes)
Settings.Nodes = 1; % Minimum number of cores to request for calculation.
Settings.Cores = -1; % Minimum number of cores to request for calculation. Set to -1 for entire node
Settings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Settings.SinglePrecision = false; % choose true for single precision mode, false for double
Settings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
Settings.npme = []; % Number of rank assigned to PME
Settings.dd = []; % Domain decomposition
Settings.openMP = false; % Domain decomposition
Continue_Chain = '3246559'; % leave blank to not continue
pyloc = '/arc/software/spack-2021/spack/opt/spack/linux-centos7-skylake_avx512/gcc-9.4.0/miniconda3-4.9.2-osd7mgcccret4btsp42lvasqoymzt4mf/bin/python';

Modules = ['# set EXE environment' newline ...
           'module purge all' newline ...
           'module load Software_Collection/2021' newline ...
           'module load gcc/9.4.0' newline ...
           'module load python/3.8.10' newline ...
           'conda activate ml' newline ...
           'export PYTHONPATH=''/arc/home/haydensc/ThesisCodeBase/LiXStructureDetector''' newline newline];

SubmDir = '/home/haydensc/scratch/Molten_Salts_MD'; % Name of project directory to contain job within the main project folder
Batch_Template = MD_Batch_Template(Settings,'sockeye');

Batch_Template = strrep(Batch_Template,'mpiprocs=32','mpiprocs=1');
Batch_Template = strrep(Batch_Template,'##DIRECTORY##',SubmDir);
Batch_Template = regexprep(Batch_Template,'# set EXE environment.+?\n\n',Modules);

cd(SubmDir)

for idx = 1:Settings.JobLinks
    
    JobName_i = ['StrucPred-' num2str(idx,'%03.f')];
    
    Batch_Text_i = strrep(Batch_Template,'##ERROR##',JobName_i);
    Batch_Text_i = strrep(Batch_Text_i,'##TASKNAME##',JobName_i);
    
    logfile_i = [JobName_i '.log'];
    Batch_Text_i = regexprep(Batch_Text_i,'##PREMIN##.+##CLEANUP##',[pyloc ' Run_Project_Structure_Classifier.py']);

    % Save batch script
    batch_file = fullfile(SubmDir,[JobName_i '.subm']);
    fidBS = fopen(batch_file,'wt');
    fwrite(fidBS,regexprep(Batch_Text_i,{'\r' Settings.wsl},{'' ''}));
    fclose(fidBS);

    
    if idx > 1
        if Settings.slurm
            [~,output] = system([Settings.qsub ' --depend=afterany:' PrevJobID ' ' batch_file]);
        else
            [~,output] = system([Settings.qsub ' -W depend=afterany:' PrevJobID ' ' batch_file]);
        end
    elseif ~isempty(Continue_Chain)
        if Settings.slurm
            [~,output] = system([Settings.qsub ' --depend=afterany:' Continue_Chain ' ' batch_file]);
        else
            [~,output] = system([Settings.qsub ' -W depend=afterany:' Continue_Chain ' ' batch_file]);
        end
    else
        [~,output] = system([Settings.qsub ' ' batch_file]);
    end
    PrevJobID = regexp(output,'[0-9]+','match','ONCE');
    disp(['Job (Link ' num2str(idx) ') submitted.'])
end



