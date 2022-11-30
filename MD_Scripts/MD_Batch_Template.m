function [TemplateText,Settings] = MD_Batch_Template(Settings,varargin)

if nargin > 1
    Server = varargin{1};
else
    [~,~,Server] = find_home;
end
Settings.postprocess = '';

switch Server
case 'sockeye' 
    Account = 'st-gpatey-1';
    new_modules = true; % Set this back to true if the bugs in the new modules get fixed
    Settings.grompp = 'grompp';
    Settings.mdrun = 'mdrun';
    Settings.g_energy = 'energy';
    Settings.g_msd = 'msd';
    Settings.trjconv = 'trjconv';
    Settings.g_check = 'check';
    Settings.genconf = 'genconf';
    Settings.g_density = 'density';
    Settings.editconf = 'editconf';
    Settings.insert_molecules = 'insert-molecules';
    Settings.convert_tpr = 'convert-tpr';
    Settings.g_traj = 'traj';
    
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
            Settings.gmx = 'gmx_mpi ';
            Settings.gmx_loc = 'gmx_mpi ';
        else
            Settings.gmx = 'gmx_mpi_d ';
            Settings.gmx_loc = 'gmx_mpi_d ';
        end
    else
        if Settings.SinglePrecision
            Settings.gmx = ['mpiexec -machinefile $PBS_NODEFILE -np ' num2str(MPI_Ranks_Per_Node*max(Settings.Nodes,1)) ...
            ' --map-by ppr:' num2str(MPI_Ranks_Per_Node) ':node gmx_mpi '];
            Settings.gmx_loc = 'gmx_mpi ';
        else
            Settings.gmx = ['mpiexec -machinefile $PBS_NODEFILE -np ' num2str(MPI_Ranks_Per_Node*max(Settings.Nodes,1)) ...
            ' --map-by ppr:' num2str(MPI_Ranks_Per_Node) ':node gmx_mpi_d '];
            Settings.gmx_loc = 'gmx_mpi_d ';
        end
    end

    Settings.mdrun_opts = '';
    if Settings.OMP_Threads > 1
        Settings.mdrun_opts = [Settings.mdrun_opts ' -ntomp ' num2str(Settings.OMP_Threads)];
    end
    if ~Settings.DLB
        Settings.mdrun_opts = [Settings.mdrun_opts ' -dlb no'];
    end
    if ~Settings.TunePME
        Settings.mdrun_opts = [Settings.mdrun_opts ' -notunepme'];
    end
    if length(Settings.dd) == 3
        Settings.mdrun_opts = [Settings.mdrun_opts ' -dd ' regexprep(num2str(Settings.dd),' +',' ')];
    end
    if ~isempty(Settings.npme)
        Settings.mdrun_opts = [Settings.mdrun_opts ' -npme ' num2str(Settings.npme)];
    end
    if abs(Settings.dds - 0.8) > sqrt(eps)
        Settings.mdrun_opts = [Settings.mdrun_opts ' -dds ' num2str(Settings.dds)];
    end
    Settings.mdrun_opts = [Settings.mdrun_opts ' -maxh ' num2str(Settings.Hours)];
    
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
            Settings.postprocess =  ['    module load Software_Collection/2021' newline ...
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
            Settings.postprocess =  ['    module purge all' newline ...
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

case {'cedar' 'graham' 'narval'} % Cedar, graham, and narval
    
    % Number of cores per node
    switch Server
    case 'cedar'
        Account = 'rrg-patey-ad';
        if Settings.Cores < 1
            if Settings.BigNode
                Cores_per_node = 48;
            else
                Cores_per_node = 32;
            end
        else
            Cores_per_node = Settings.Cores/max(Settings.Nodes,1);
        end
    case 'graham'
        Account = 'def-patey';
        if Settings.Cores < 1
            Cores_per_node = 32; % Graham
        else
            Cores_per_node = Settings.Cores/max(Settings.Nodes,1);
        end
    case 'narval'
        Account = 'def-patey';
        if Settings.Cores < 1
            Cores_per_node = 64; % Narval
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
    
    if strcmp(Settings.gmx_version,'4.6.7')
        gmx_old_version = true;
    elseif strcmp(Settings.gmx_version,'2019')
        gmx_old_version = false;
    else
        gmx_old_version = Settings.Polarization;
    end
    
    Settings.mdrun_opts = '';
    if Settings.Nodes > 0 % Whole node or Multi-node calculation
        nodeline = ['#SBATCH --nodes=' num2str(Settings.Nodes) newline ...
                    '#SBATCH --exclusive' newline];
        tasksline = ['#SBATCH --tasks-per-node=' num2str(MPI_Ranks_Per_Node) newline ...
                     '#SBATCH --cpus-per-task=' num2str(Settings.OMP_Threads) newline];
        
        % Deal with the executable to call
        if gmx_old_version % Use version 4.6.7
            Settings.gmx = ['mpiexec -np ' num2str(MPI_Ranks_Per_Node*Settings.Nodes) ' '];
            Settings.gmx_loc = '';
        elseif Settings.MPI_Ranks == 1
            if Settings.SinglePrecision
                Settings.gmx = 'gmx ';
                Settings.gmx_loc = 'gmx ';
            else
                Settings.gmx = 'gmx_d ';
                Settings.gmx_loc = 'gmx_d ';
            end
            Settings.mdrun_opts = [Settings.mdrun_opts ' -ntmpi 1'];
        else
            if Settings.SinglePrecision
                Settings.gmx = ['mpiexec -np ' num2str(MPI_Ranks_Per_Node*Settings.Nodes) ' gmx_mpi '];
                Settings.gmx_loc = 'gmx ';
            else
                Settings.gmx = ['mpiexec -np ' num2str(MPI_Ranks_Per_Node*Settings.Nodes) ' gmx_mpi_d '];
                Settings.gmx_loc = 'gmx_d ';
            end
        end
        
    else % Partial node calculation
        nodeline = '';
        tasksline = ['#SBATCH --tasks=' num2str(MPI_Ranks_Per_Node) newline ...
                     '#SBATCH --cpus-per-task=' num2str(Settings.OMP_Threads) newline];
        
        % Deal with the executable to call
        if gmx_old_version % Use version 4.6.7
            Settings.gmx = ['mpiexec -np ' num2str(MPI_Ranks_Per_Node) ' '];
            Settings.gmx_loc = '';
        elseif Settings.MPI_Ranks == 1
            if Settings.SinglePrecision
                Settings.gmx = 'gmx ';
                Settings.gmx_loc = 'gmx ';
            else
                Settings.gmx = 'gmx_d ';
                Settings.gmx_loc = 'gmx_d ';
            end
            Settings.mdrun_opts = [Settings.mdrun_opts ' -ntmpi 1'];
        else
            if Settings.SinglePrecision
                Settings.gmx = ['mpiexec -np ' num2str(MPI_Ranks_Per_Node) ' gmx_mpi '];
                Settings.gmx_loc = 'gmx ';
            else
                Settings.gmx = ['mpiexec -np ' num2str(MPI_Ranks_Per_Node) ' gmx_mpi_d '];
                Settings.gmx_loc = 'gmx_d ';
            end
        end
    end
    
    if gmx_old_version
        Settings.grompp = 'grompp_mpi_d';
        Settings.mdrun = 'mdrun_mpi_d';
        Settings.g_energy = 'g_energy_mpi_d';
        Settings.g_msd = 'g_msd_mpi_d';
        Settings.trjconv = 'trjconv_mpi_d';
        Settings.g_check = 'gmxcheck_mpi_d';
        Settings.genconf = 'genconf_mpi_d';
        Settings.g_density = 'g_density_mpi_d';
        Settings.editconf = 'editconf_mpi_d';
        Settings.g_traj = 'g_traj_mpi_d';
        Settings.convert_tpr = 'tpbconv_mpi_d';
        if Settings.SinglePrecision
            Settings.insert_molecules = 'gmx_d insert-molecules';
        else
            Settings.insert_molecules = 'gmx insert-molecules';
        end
        Settings.mdrun_opts = [Settings.mdrun_opts ' -pd'];
    else
        Settings.grompp = 'grompp';
        Settings.mdrun = 'mdrun';
        Settings.g_energy = 'energy';
        Settings.g_msd = 'msd';
        Settings.trjconv = 'trjconv';
        Settings.g_check = 'check';
        Settings.genconf = 'genconf';
        Settings.g_density = 'density';
        Settings.editconf = 'editconf';
        Settings.insert_molecules = 'insert-molecules';
        Settings.convert_tpr = 'convert-tpr';
        Settings.g_traj = 'traj';
    end
    if Settings.OMP_Threads > 1
        Settings.mdrun_opts = [Settings.mdrun_opts ' -ntomp ' num2str(Settings.OMP_Threads)];
    end
    if ~Settings.DLB
        Settings.mdrun_opts = [Settings.mdrun_opts ' -dlb no'];
    end
    if ~Settings.TunePME
        Settings.mdrun_opts = [Settings.mdrun_opts ' -notunepme'];
    end
    if length(Settings.dd) == 3
        Settings.mdrun_opts = [Settings.mdrun_opts ' -dd ' regexprep(num2str(Settings.dd),' +',' ')];
    end
    if ~isempty(Settings.npme)
        Settings.mdrun_opts = [Settings.mdrun_opts ' -npme ' num2str(Settings.npme)];
    end
    if abs(Settings.dds - 0.8) > sqrt(eps)
        Settings.mdrun_opts = [Settings.mdrun_opts ' -dds ' num2str(Settings.dds)];
    end
    Settings.mdrun_opts = [Settings.mdrun_opts ' -maxh ' num2str(Settings.Hours)];
    
    % Deal with memory allocation
    if strcmp(Settings.Mempernode,'-1') % Default
        memline = '';
    elseif strcmp(Settings.Mempernode,'0') % Max memory per CPU
        memline = ['#SBATCH --mem-per-cpu=3800M' newline];
    else % Custom memory per node
        memline = ['#SBATCH --mem=' Settings.Mempernode newline];
    end
    
    if gmx_old_version
        Gromacs_Module = ['module load gromacs-plumed/2019.6' newline ...
                          'module load fftw/3.3.8' newline ...
                          'source /home/scheiber/.local/bin/GMXRC' newline];
    else
        Gromacs_Module = ['module load gromacs-plumed/2019.6' newline];
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
        Gromacs_Module ...
        'module load tbb/2020.2 python/3.8' newline ...
        'module load matlab/2022a' newline ...
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
otherwise
    if ispc
        Settings.gmx = 'wsl source ~/.bashrc; gmx_d ';
    else
        if Settings.SinglePrecision
            Settings.gmx = 'gmx ';
        else
            Settings.gmx = 'gmx_d ';
        end
    end
    Settings.gmx_loc = Settings.gmx;
    Settings.grompp = 'grompp';
    Settings.mdrun = 'mdrun';
    Settings.g_energy = 'energy';
    Settings.g_msd = 'msd';
    Settings.trjconv = 'trjconv';
    Settings.g_check = 'check';
    Settings.genconf = 'genconf';
    Settings.g_density = 'density';
    Settings.editconf = 'editconf';
    Settings.insert_molecules = 'insert-molecules';
    Settings.convert_tpr = 'convert-tpr';
    Settings.g_traj = 'traj';
    
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
    
    Settings.mdrun_opts = ' -v';
    Settings.mdrun_opts = [Settings.mdrun_opts ' -ntmpi ' num2str(MPI_Ranks_Per_Node)];
    Settings.mdrun_opts = [Settings.mdrun_opts ' -ntomp ' num2str(Settings.OMP_Threads)];
    if MPI_Ranks_Per_Node == 1 && Settings.OMP_Threads == 1
        Settings.mdrun_opts = [Settings.mdrun_opts ' -pin on'];
    end
    
    if ~Settings.DLB
        Settings.mdrun_opts = [Settings.mdrun_opts ' -dlb no'];
    end
    if ~Settings.TunePME
        Settings.mdrun_opts = [Settings.mdrun_opts ' -notunepme'];
    end
    if length(Settings.dd) == 3
        Settings.mdrun_opts = [Settings.mdrun_opts ' -dd ' regexprep(num2str(Settings.dd),' +',' ')];
    end
    if ~isempty(Settings.npme)
        Settings.mdrun_opts = [Settings.mdrun_opts ' -npme ' num2str(Settings.npme)];
    end
    if abs(Settings.dds - 0.8) > sqrt(eps)
        Settings.mdrun_opts = [Settings.mdrun_opts ' -dds ' num2str(Settings.dds)];
    end
    
    [TemplateText,~] = MD_Batch_Template(Settings,'cedar');
end

end