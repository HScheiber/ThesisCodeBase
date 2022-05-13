function [TemplateText,gmx,gmx_loc,mdrun_opts,postprocess] = MD_Batch_Template(Settings,varargin)

if nargin > 1
    Server = varargin{1};
else
    if ispc
        Server = getenv('COMPUTERNAME');
    else
        [~,Servertxt] = system('hostname -s | cut -c 1-3');
        Server = strtrim(Servertxt);
    end
end

if strcmpi(Server,'ced') || strcmpi(Server,'cdr') % cedar
    Account = 'rrg-patey-ad';
elseif ~isempty(regexp(Server,'se[0-9]','ONCE')) || strcmpi(Server,'log') % sockeye
    Account = 'st-gpatey-1';
else
    Account = 'def-patey';
end

postprocess = '';

% if strcmpi(Server,'sea') || strcmpi(Server,'pod') % Orcinus
%     if Settings.Cores < 0
%         Cores_per_node = 12;
%     else
%         Cores_per_node = Settings.Cores;
%     end
%     
%     if strcmp(Settings.Mempernode,'-1')
%         memline = '';
%     elseif strcmp(Settings.Mempernode,'0')
%         memline = '';
%     else
%         memline = ['#PBS -l pmem=' Settings.Mempernode newline];
%     end
%     
%     if Settings.Nodes > 1 && Settings.openMP
%         mdrun_opts = ['mpiexec -np ' num2str(Settings.Nodes*2) ...
%             ' --npernode $MPI_PPN --mca mpi_paffinity_alone 0 ##MDRUN## -v -notunepme -ntomp $OMP_NUM_THREADS -maxh ' ...
%             num2str(Settings.Hours)];
%         if Settings.SinglePrecision
%             gmx = 'gmx_mpi';
%             gmx_loc = 'gmx';
%         else
%             gmx = 'gmx_mpi_d';
%             gmx_loc = 'gmx_d';
%         end
%         MPI_PPN = 2;
%         OMP_NUM_THREADS = 6;
%     elseif Settings.Nodes > 1 && ~Settings.openMP
%         mdrun_opts = ['mpiexec -machinefile $PBS_NODEFILE -np ' num2str(Settings.Nodes*Cores_per_node) ...
%             ' --npernode $MPI_PPN --mca mpi_paffinity_alone 0 ##MDRUN## -v -notunepme -ntomp $OMP_NUM_THREADS -maxh ' ...
%             num2str(Settings.Hours)];
%         if Settings.SinglePrecision
%             gmx = 'gmx_mpi';
%             gmx_loc = 'gmx';
%         else
%             gmx = 'gmx_mpi_d';
%             gmx_loc = 'gmx_d';
%         end
%         MPI_PPN = Cores_per_node;
%         OMP_NUM_THREADS = 1;
%     else
%         mdrun_opts = ['mpiexec -machinefile $PBS_NODEFILE -np ' num2str(Cores_per_node) ...
%             ' --mca mpi_paffinity_alone 0 ##MDRUN## -v -notunepme -ntomp $OMP_NUM_THREADS -maxh ' ...
%             num2str(Settings.Hours)];
%         if Settings.SinglePrecision
%             gmx = 'gmx_mpi';
%             gmx_loc = 'gmx';
%         else
%             gmx = 'gmx_mpi_d';
%             gmx_loc = 'gmx_d';
%         end
%         MPI_PPN = Cores_per_node;
%         OMP_NUM_THREADS = 1;
%     end
%     
%     TemplateText = ['#!/bin/bash' newline ...
%         '#PBS -S /bin/bash' newline ...
%         '#PBS -l walltime=' num2str(Settings.Hours) ':' num2str(Settings.Mins) ':00' newline ...
%         '#PBS -l nodes=' num2str(Settings.Nodes) ':ppn=' num2str(Cores_per_node) newline ...
%         memline ...
%         '#PBS -V' newline ...
%         '#PBS -N ##TASKNAME##' newline ...
%         '#PBS -e ##ERROR##.stde' newline ...
%         '#PBS -o ##ERROR##.stdo' newline ...
%         '#PBS -l partition=QDR' newline ...
%         newline newline ...
%         '# Check on some basics:' newline ...
%         'echo "Running on host: " `hostname`' newline ...
%         'echo "Changing to directory from which PBS script was submitted."' newline ...
%         'cd ##DIRECTORY##' newline ...
%         'echo "Current working directory is now: " `pwd`' newline ...
%         newline newline ...
%         '# set EXE environment' newline ...
%         'module load gromacs/5.1.4' newline ...
%         'module load matlab/matlab_2015b' newline ...
%         'export MPI_PPN=' num2str(MPI_PPN) newline ...
%         'export OMP_NUM_THREADS=' num2str(OMP_NUM_THREADS) newline ...
%         newline newline ...
%         '# Run Job' newline ...
%         '##PREMIN##' newline ...
%         '##EXT1##' newline ...
%         mdrun_opts newline ...
%         '##EXT2##' newline ...
%         '##CLEANUP##' newline ...
%         'echo "Job completed at `date`"' newline ... 
%         'exit 0'];
if ~isempty(regexp(Server,'se[0-9]','ONCE')) || strcmpi(Server,'log') % sockeye
    
    new_modules = true; % Set this back to true if the bugs in the new modules get fixed
    
    if Settings.Cores < 1
        if Settings.BigNode
            Cores_per_node = 40;
            cpu_arch = ':cpu_arch=cascade';
        else
            Cores_per_node = 32;
            cpu_arch = ':cpu_arch=skylake';
        end
    else
        Cores_per_node = Settings.Cores/max(Settings.Nodes,1);
        cpu_arch = '';
    end
    
    % Deal with MPI Ranks and OMP Threads
    if Settings.MPI_Ranks < 0
        MPI_Ranks_Per_Node = Cores_per_node/Settings.OMP_Threads;
    else
        MPI_Ranks_Per_Node = Settings.MPI_Ranks;
    end
    
    % Deal with memory
    if strcmp(Settings.Mempernode,'-1') || strcmp(Settings.Mempernode,'0')
        memline = '184gb';
    else
        memline = Settings.Mempernode;
    end
    
    % Deal with the executable to call
    if Settings.MPI_Ranks == 1
        if Settings.SinglePrecision
            gmx = 'gmx_mpi';
            gmx_loc = 'gmx_mpi';
        else
            gmx = 'gmx_mpi_d';
            gmx_loc = 'gmx_mpi_d';
        end
    else
        if Settings.SinglePrecision
            gmx = ['mpiexec -machinefile $PBS_NODEFILE -np ' num2str(MPI_Ranks_Per_Node*max(Settings.Nodes,1)) ...
            ' --map-by ppr:' num2str(MPI_Ranks_Per_Node) ':node gmx_mpi'];
            gmx_loc = 'gmx_mpi';
        else
            gmx = ['mpiexec -machinefile $PBS_NODEFILE -np ' num2str(MPI_Ranks_Per_Node*max(Settings.Nodes,1)) ...
            ' --map-by ppr:' num2str(MPI_Ranks_Per_Node) ':node gmx_mpi_d'];
            gmx_loc = 'gmx_mpi_d';
        end
    end

    mdrun_opts = '';
    if Settings.OMP_Threads > 1
        mdrun_opts = [mdrun_opts ' -ntomp ' num2str(Settings.OMP_Threads)];
    end
    if ~Settings.DLB
        mdrun_opts = [mdrun_opts ' -dlb no'];
    end
    if ~Settings.TunePME
        mdrun_opts = [mdrun_opts ' -notunepme'];
    end
    if length(Settings.dd) == 3
        mdrun_opts = [mdrun_opts ' -dd ' regexprep(num2str(Settings.dd),' +',' ')];
    end
    if ~isempty(Settings.npme)
        mdrun_opts = [mdrun_opts ' -npme ' num2str(Settings.npme)];
    end
    if abs(Settings.dds - 0.8) > sqrt(eps)
        mdrun_opts = [mdrun_opts ' -dds ' num2str(Settings.dds)];
    end
    mdrun_opts = [mdrun_opts ' -maxh ' num2str(Settings.Hours)];
    
    if Settings.SinglePrecision
        if new_modules
            dependencies = ['module load Software_Collection/2021' newline ...
                            'module load gcc/9.4.0' newline ...
                            'module load openmpi/4.1.1-cuda11-3' newline ...
                            'module load openblas/0.3.15' newline ...
                            'module load python/3.8.10' newline];
            gmx_module = 'gromacs/2019.2';
            matlab_module = 'matlab/R2021a';
        else
            postprocess =  ['    module load Software_Collection/2021' newline ...
                            '    module load gcc/9.4.0' newline ...
                            '    module load python/3.8.10' newline ...
                            '    module load matlab/R2021a' newline];
            
            dependencies = ['module load gcc/9.1.0' newline ...
                            'module load openmpi/3.1.4' newline ...
                            'module load fftw/3.3.8'];
            gmx_module = 'gromacs/2019.2';
            matlab_module = 'matlab/R2018b';
        end
    else
        if new_modules
            dependencies = ['module load Software_Collection/2021' newline ...
                            'module load gcc/9.4.0' newline ...
                            'module load openmpi/4.1.1-cuda11-3' newline ...
                            'module load openblas/0.3.15' newline ...
                            'module load python/3.8.10' newline];
            
            gmx_module = 'gromacs/2019.6-double';
            matlab_module = 'matlab/R2021a';
        else
            postprocess =  ['    module purge all' newline ...
                            '    module load Software_Collection/2021' newline ...
                            '    module load gcc/9.4.0' newline ...
                            '    module load python/3.8.10' newline ...
                            '    module load matlab/R2021a' newline];
            
            dependencies = ['module load gcc/9.1.0' newline ...
                            'module load openmpi/3.1.4' newline ...
                            'module load fftw/3.3.8'];
            gmx_module = 'gromacs/5.1.4-gr5.1.4';
            matlab_module = 'matlab/R2018b';
        end
    end
    
    TemplateText = ['#!/bin/bash' newline ...
        '#PBS -l walltime=' num2str(Settings.Hours) ':' num2str(Settings.Mins) ':00,' ...
        'select=' num2str(max(Settings.Nodes,1)) ':ncpus=' num2str(Cores_per_node) ':mpiprocs=' num2str(MPI_Ranks_Per_Node) ...
        ':ompthreads=' num2str(Settings.OMP_Threads) ':mem=' memline cpu_arch newline ...
        '#PBS -A ' Account newline ...
        '#PBS -N ##TASKNAME##' newline ...
        '#PBS -e ##ERROR##.stde' newline ...
        '#PBS -o ##ERROR##.stdo' newline ...
        newline newline ...
        '# Check on some basics:' newline ...
        'echo "Running on host: " `hostname`' newline ...
        'echo "Changing to directory from which PBS script was submitted."' newline ...
        'cd ##DIRECTORY##' newline ...
        'echo "Current working directory is now: " `pwd`' newline ...
        newline newline ...
        '# set EXE environment' newline ...
        dependencies newline ...
        'module load ' gmx_module newline ...
        'module load ' matlab_module newline ...
        'export MPI_PPN=' num2str(MPI_Ranks_Per_Node) newline ...
        'mkdir $TMPDIR/.matlab' newline ...
        'export MATLAB_PREFDIR=$TMPDIR/.matlab' newline ...
        newline newline ...
        '# Run Job' newline ...
        '##PREMIN##' newline ...
        '##EXT1##' newline ...
        '##MDRUN##' newline ...
        '##EXT2##' newline ...
        '##CLEANUP##' newline ...
        'echo "Job completed at `date`"' newline ... 
        'exit 0'];
elseif strcmpi(Server,'ced') || strcmpi(Server,'cdr') || strcmpi(Server,'gra') % Cedar and graham
    
    % Number of cores per node
    if strcmpi(Server,'ced') || strcmpi(Server,'cdr') % Cedar
        if Settings.Cores < 1
            if Settings.BigNode
                Cores_per_node = 48;
            else
                Cores_per_node = 32;
            end
        else
            Cores_per_node = Settings.Cores/max(Settings.Nodes,1);
        end
    else
        if Settings.Cores < 1
            Cores_per_node = 32; % Graham or testing pc
        else
            Cores_per_node = Settings.Cores/max(Settings.Nodes,1);
        end
    end
    
    % Deal with MPI Ranks and OMP Threads
    if Settings.MPI_Ranks < 0
        MPI_Ranks_Per_Node = Cores_per_node/Settings.OMP_Threads;
    else
        MPI_Ranks_Per_Node = Settings.MPI_Ranks;
    end
    
    if Settings.Nodes > 0 % Whole node or Multi-node calculation
        nodeline = ['#SBATCH --nodes=' num2str(Settings.Nodes) newline ...
                    '#SBATCH --exclusive' newline];
        tasksline = ['#SBATCH --tasks-per-node=' num2str(MPI_Ranks_Per_Node) newline ...
                     '#SBATCH --cpus-per-task=' num2str(Settings.OMP_Threads) newline];
        
        % Deal with the executable to call
        if Settings.MPI_Ranks == 1
            if Settings.SinglePrecision
                gmx = 'gmx';
                gmx_loc = 'gmx';
            else
                gmx = 'gmx_d';
                gmx_loc = 'gmx_d';
            end
        else
            if Settings.SinglePrecision
                gmx = ['mpiexec -np ' num2str(MPI_Ranks_Per_Node*Settings.Nodes) ' gmx_mpi'];
                gmx_loc = 'gmx';
            else
                gmx = ['mpiexec -np ' num2str(MPI_Ranks_Per_Node*Settings.Nodes) ' gmx_mpi_d'];
                gmx_loc = 'gmx_d';
            end
        end
        
    else % Partial node calculation
        nodeline = '';
        tasksline = ['#SBATCH --tasks=' num2str(MPI_Ranks_Per_Node) newline ...
                     '#SBATCH --cpus-per-task=' num2str(Settings.OMP_Threads) newline];
        
        % Deal with the executable to call
        if Settings.MPI_Ranks == 1
            if Settings.SinglePrecision
                gmx = 'gmx';
                gmx_loc = 'gmx';
            else
                gmx = 'gmx_d';
                gmx_loc = 'gmx_d';
            end
        else
            if Settings.SinglePrecision
                gmx = ['mpiexec -np ' num2str(MPI_Ranks_Per_Node) ' gmx_mpi'];
                gmx_loc = 'gmx';
            else
                gmx = ['mpiexec -np ' num2str(MPI_Ranks_Per_Node) ' gmx_mpi_d'];
                gmx_loc = 'gmx_d';
            end
        end
    end
    
    mdrun_opts = '';
    if Settings.OMP_Threads > 1
        mdrun_opts = [mdrun_opts ' -ntomp ' num2str(Settings.OMP_Threads)];
    end
    if ~Settings.DLB
        mdrun_opts = [mdrun_opts ' -dlb no'];
    end
    if ~Settings.TunePME
        mdrun_opts = [mdrun_opts ' -notunepme'];
    end
    if length(Settings.dd) == 3
        mdrun_opts = [mdrun_opts ' -dd ' regexprep(num2str(Settings.dd),' +',' ')];
    end
    if ~isempty(Settings.npme)
        mdrun_opts = [mdrun_opts ' -npme ' num2str(Settings.npme)];
    end
    if abs(Settings.dds - 0.8) > sqrt(eps)
        mdrun_opts = [mdrun_opts ' -dds ' num2str(Settings.dds)];
    end
    mdrun_opts = [mdrun_opts ' -maxh ' num2str(Settings.Hours)];
    
    % Deal with memory allocation
    if strcmp(Settings.Mempernode,'-1') % Default
        memline = '';
    elseif strcmp(Settings.Mempernode,'0') % Max memory per CPU
        memline = ['#SBATCH --mem-per-cpu=3800M' newline];
    else % Custom memory per node
        memline = ['#SBATCH --mem=' Settings.Mempernode newline];
    end
    
    TemplateText = ['#!/bin/bash' newline ... 
        '#SBATCH --time=' num2str(Settings.Hours) ':' num2str(Settings.Mins) ':00' newline ... 
        nodeline ...
        tasksline ... 
        memline ...
        '#SBATCH --account=' Account newline ... 
        '#SBATCH --job-name=##TASKNAME##' newline ... 
        '#SBATCH --output=##ERROR##.stdo' newline ... 
        '#SBATCH --error=##ERROR##.stde' newline ... 
        '#SBATCH --export=ALL' newline ... 
        newline newline ... 
        '# Check on some basics:' newline ... 
        'echo "Running on host: " `hostname`' newline ... 
        'echo "Changing to directory from which PBS script was submitted."' newline ... 
        'cd ##DIRECTORY##' newline ... 
        'echo "Current working directory is now: " `pwd`' newline ... 
        newline newline ... 
        '# Load modules' newline ...
        'module purge --all' newline ...
        'source $HOME/ml/bin/activate' newline ...
        'module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3' newline ...
        'module load gromacs-plumed/2019.6' newline ...
        'module load tbb/2020.2 python/3.8' newline ...
        'module load matlab/2021a.1' newline ...
        newline newline...
        '# Set variables' newline ...
        'export OMP_NUM_THREADS=' num2str(Settings.OMP_Threads) newline ...
        'mkdir $SLURM_TMPDIR/.matlab' newline ...
        'export MATLAB_PREFDIR=$SLURM_TMPDIR/.matlab' newline ...
        newline newline ...
        '# Run Job' newline ...
        '##PREMIN##' newline ...
        '##EXT1##' newline ...
        '##MDRUN##' newline ...
        '##EXT2##' newline ...
        '##CLEANUP##' newline ...
        'echo "Job completed at `date`"' newline ... 
        'exit 0'];
else
    if ispc
        gmx = 'wsl source ~/.bashrc; gmx_d';
    else
        if Settings.SinglePrecision
            gmx = 'gmx';
        else
            gmx = 'gmx_d';
        end
    end
    gmx_loc = gmx;
    
    if Settings.Cores < 1
        Cores_per_node = feature('numcores');
    else
        Cores_per_node = Settings.Cores;
    end
    
    % Deal with MPI Ranks and OMP Threads
    if Settings.MPI_Ranks < 0
        MPI_Ranks_Per_Node = Cores_per_node/Settings.OMP_Threads;
    else
        MPI_Ranks_Per_Node = Settings.MPI_Ranks;
    end
    
    mdrun_opts = ' -v';
    if MPI_Ranks_Per_Node > 1
        mdrun_opts = [mdrun_opts ' -ntmpi ' num2str(MPI_Ranks_Per_Node)];
    end
    if Settings.OMP_Threads > 1
        mdrun_opts = [mdrun_opts ' -ntomp ' num2str(Settings.OMP_Threads)];
    end
    if ~Settings.DLB
        mdrun_opts = [mdrun_opts ' -dlb no'];
    end
    if ~Settings.TunePME
        mdrun_opts = [mdrun_opts ' -notunepme'];
    end
    if length(Settings.dd) == 3
        mdrun_opts = [mdrun_opts ' -dd ' regexprep(num2str(Settings.dd),' +',' ')];
    end
    if ~isempty(Settings.npme)
        mdrun_opts = [mdrun_opts ' -npme ' num2str(Settings.npme)];
    end
    if abs(Settings.dds - 0.8) > sqrt(eps)
        mdrun_opts = [mdrun_opts ' -dds ' num2str(Settings.dds)];
    end
    
    [TemplateText,~,~] = MD_Batch_Template(Settings,'ced');
end

end