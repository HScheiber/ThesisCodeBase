function [home,project,computer,slurm,BO_Models,qsub,passlog,pipe,wsl,StrSelModels,scratch] = find_home
    if ispc
        home = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase';
        project = 'C:\Users\Hayden\Documents\Patey_Lab';
        scratch = 'C:\Users\Hayden\Documents\Patey_Lab';
        BO_Models = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';
		StrSelModels = 'C:\Users\Hayden\Documents\Patey_Lab\Machine_Learning_Structure_Finder\LiX_Structure_Selector_Optimization';
        computer = lower(getenv('COMPUTERNAME'));
        slurm = false;
        qsub = '';
        passlog = ' ^&^> ';
        pipe = '^|';
        wsl = 'wsl source ~/.bashrc; ';
    elseif isunix
        Servertxt = lower(getenv('HOSTNAME'));
        if isempty(Servertxt)
            [~,Servertxt] = system('echo $HOSTNAME');
            Servertxt = lower(Servertxt);
        end
        server = Servertxt(1:3);
        passlog = ' &> ';
        pipe = '|';
        wsl = '';
        
        if strcmp(server,'ced') || strcmp(server,'cdr')
            home = '/home/scheiber/ThesisCodeBase';     
            project = '/project/6001647/scheiber';
            scratch = '/scratch/scheiber';
            BO_Models = '/home/scheiber/project/Model_Building/Completed';
            computer = 'cedar';
            slurm = true;
            qsub = 'sbatch';
			StrSelModels = '/home/scheiber/project/LiX_Structure_Selector';
        elseif strcmp(server,'nar') || ~isempty(regexp(server,'nc[0-9]','ONCE'))
            home = '/home/scheiber/ThesisCodeBase';     
            project = '/project/6001647/scheiber';
            computer = 'narval';
            scratch = '/lustre07/scratch/scheiber';
            slurm = true;
            BO_Models = '/home/scheiber/project/Model_Building/Completed';
            qsub = 'sbatch';
			StrSelModels = '/home/scheiber/project/LiX_Structure_Selector';
        elseif ~isempty(regexp(server,'se[0-9]','ONCE')) || strcmpi(server,'log')
            home = '/home/haydensc/ThesisCodeBase';
            project = '/scratch/st-gpatey-1/haydensc';
            BO_Models = '/home/haydensc/project/BO_Models';
            scratch = '/scratch/st-gpatey-1/haydensc';
            computer = 'sockeye';
            slurm = false;
            qsub = 'qsub';
			StrSelModels = '/home/haydensc/project/LiX_Structure_Selector';
        elseif strcmp(server,'pat')
            home = '/home/user/ThesisCodeBase';
            project = '/home/user/project';
            scratch = '/home/user/scratch';
            BO_Models = '/home/user/project/BO_Models';
            computer = 'patey';
            slurm = false;
            qsub = 'local';
			StrSelModels = '/home/user/project/LiX_Structure_Selector';
        elseif strcmp(server,'gra')
            home = '/home/scheiber/ThesisCodeBase';     
            project = '/project/6001647/scheiber';
            scratch = '/scratch/scheiber';
            computer = 'graham';
            slurm = true;
            BO_Models = '/home/scheiber/project/Model_Building/Completed';
            qsub = 'sbatch';
			StrSelModels = '/home/scheiber/project/LiX_Structure_Selector';
        elseif strcmp(server,'sea')
            home = '/home/scheiber/ThesisCodeBase';     
            project = '/home/scheiber/project';
            scratch = '/home/scheiber/scratch';
            computer = 'orcinus';
            slurm = false;
            BO_Models = '/home/scheiber/project/BO_Models';
            qsub = 'qsub';
			StrSelModels = '/home/scheiber/project/LiX_Structure_Selector';
        else
            home = '/home/scheiber/ThesisCodeBase';
            project = '/project/6001647/scheiber';
            scratch = '/scratch/scheiber';
            computer = lower(system('hostname'));
            slurm = true;
			BO_Models = '/home/scheiber/project/BO_Models';
            qsub = 'sbatch';
			StrSelModels = '/home/scheiber/project/LiX_Structure_Selector';
        end
    end
end