% Submit_deriv_free_opt
%% List of theories %%
% 'PBE' 'optB88' 'optB88-vdW' 'PBE-D3' 
% 'PBE-rVV10' 'PBE-rVV10L' 'PBEsol' 
% 'PBEsol-rVV10' 'optPBE' 'optPBE-vdW'
% 'B86r' 'rev-vdW-DF2' 'PW86r' 'vdW-DF2' 
% 'BEEF-vdW' 'SG4' 'SG4-rVV10m' 'SCAN' 
% 'SCAN-rVV10' 'TMTPSS' 'TMTPSS-rVV10L'
% 'HF' 'None'
%% Calculation Parameters
Salts = {'LiBr'};
Structures = {'BetaBeO'}; %  'Wurtzite' 'NiAs'
Theories = {'TMTPSS-rVV10L'};
Basis_Experiment = 'Sapporo-QZP';
Cutoff = 800;
vdw_cutoff = 600; % Only applies when using vdW or rVV10 corrections
kpoints = 7; % does not apply to hybrids
Eps_default = 1e-15;
Eps_pgf_orb = 1e-20;
Eps_SCF = 1e-5;
finer_grid = true; % Use finer grid for the XC
Hours = 48; % number of hours to give calculation
MPIcores = 9; % Number of MPI cores to use for calculation
OMPcores = 1; % Number of OMP threads to use for calculation
Cores_total = 32; % Total number of CPU cores to request
restart = true; % Checks if job complete, if not then finds restart file and restarts the job
skip_running = true; % checks if a job is already running and skips if it is.
Increase_Radial_Grid = 0; % Increase radial grid size by this amount
Large_memory = true;

%% Script begins here
[home,project,computer] = find_home;
CP2K_Scripts = [home filesep 'CP2K_scripts'];
switch computer
    case 'unbearabull'
        pdir = pwd;
        scheduler = '_slurm';
    case 'cedar'
        pdir = [project filesep 'CP2K_Calcs' filesep Basis_Experiment];
        scheduler = '_slurm';
    case 'sockeye'
        pdir = [project filesep 'CP2K_Calcs' filesep Basis_Experiment];
        scheduler = '_pbs';
    case 'patey'
        pdir = [project filesep 'CP2K_Calcs' filesep Basis_Experiment];
        scheduler = 'local';
    otherwise
        pdir = [project filesep 'CP2K_Calcs' filesep Basis_Experiment];
        scheduler = '_slurm';
end

% Load the submission template
if strcmp(scheduler,'local')
    subtext = '';
else
    subtext = fileread([CP2K_Scripts filesep 'deriv_free_submission_template' scheduler '.subm']);

    % Add in some global calculation parameters to the submission script
    subtext = strrep(subtext,'##HOURS##',num2str(Hours));
    subtext = strrep(subtext,'##CORESTOT##',num2str(Cores_total));

    if Large_memory
        subtext = strrep(subtext,'##MEMORY##','16000M');
    else
        subtext = strrep(subtext,'##MEMORY##','3800M');
    end
end

% Add in global calculation parameters to the matlab optimizer script
calctext = fileread([CP2K_Scripts filesep 'CP2K_deriv_free_optimizer.m']);
calctext = strrep(calctext,'##MPICORES##',num2str(MPIcores));
calctext = strrep(calctext,'##OMPCORES##',num2str(OMPcores));
calctext = strrep(calctext,'##CUTOFF##',num2str(Cutoff));
calctext = strrep(calctext,'##VDWCUTOFF##',num2str(vdw_cutoff));
calctext = strrep(calctext,'##KPOINTS##',num2str(kpoints));
calctext = strrep(calctext,'##BASIS_EXPERIMENT##',Basis_Experiment);

calctext = strrep(calctext,'##EPS_DEFAULT##',num2str(Eps_default));
calctext = strrep(calctext,'##EPS_PGF_ORB##',num2str(Eps_pgf_orb));
calctext = strrep(calctext,'##EPS_SCF##',num2str(Eps_SCF));

if finer_grid
    calctext = strrep(calctext,'##FINERGRID##','true');
else
    calctext = strrep(calctext,'##FINERGRID##','false');
end


for idx = 1:length(Theories)
    Theory = Theories{idx};
    
    for jdx = 1:length(Structures)
        Structure = Structures{jdx};
        
        for kdx = 1:length(Salts)
            Salt = Salts{kdx};
            
            % Make sure to use
            switch Salt
                case {'LiBr' 'LiI' 'LiAt'}
                    Theory_Dir = [Theory '-DKH'];
                otherwise
                    Theory_Dir = Theory;
            end
            
            % Submission directory
            cdir = [pdir filesep Theory_Dir filesep Structure filesep Salt];
            if ~exist(cdir, 'dir')
               mkdir(cdir)
            end
            
            % Job name
            JobName = [Theory_Dir '_' Structure '_' Salt];
            
            % Copy over the submission template
            subtext_copy = subtext;
            calctext_copy = calctext;
            
            % Add in the job name and location to submission template
            subtext_copy = strrep(subtext_copy,'##JOBNAME##',strrep(JobName,'-','_'));
            subtext_copy = strrep(subtext_copy,'##CDIR##',cdir);
            
            % Prepare the job input text
            calctext_copy = strrep(calctext_copy,'##SALT##',Salt);
            calctext_copy = strrep(calctext_copy,'##STRUCTURE##',Structure);
            calctext_copy = strrep(calctext_copy,'##THEORY##',Theory);
            calctext_copy = strrep(calctext_copy,'##CDIR##',cdir);
            
            % Save input and submission files to the working directory
            inp_text_filename = fullfile(cdir,[strrep(JobName,'-','_') '.m']);
            fid = fopen(inp_text_filename,'w');
            fprintf(fid,'%s',calctext_copy);
            fclose(fid);
            
            if ~strcmp(scheduler,'local')
                submit_text_filename = fullfile(cdir,[JobName '.subm']);
                fid = fopen(submit_text_filename,'w');
                fprintf(fid,'%s',subtext_copy);
                fclose(fid);
            end
            
            % Submit the batch job
            cd(cdir)
            disp(['submitting ' JobName])
            if ~ispc
                switch scheduler
                    case '_pbs'
                        system(['qsub ' submit_text_filename]);
                    case '_slurm'
                        system(['sbatch ' submit_text_filename]);
                    case 'local'
                        dname = [JobName '.optlog'];
                        diary(dname)
                        evalc(strrep(JobName,'-','_'))
                        diary off
                end
                pause(0.2)
            end
        end
    end
end
cd(pdir)