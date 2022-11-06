% Collect_MP_Data

%% Initialize directories of input and output data
Settings = Initialize_MD_Settings;
Settings.Project_Directory_Name = 'Alexandria_Polarized_Model_MPs';
DataSetName = 'Alexandria_Melting_Point_Data.mat';
DataKeyword = '';
ProjectDir = fullfile(Settings.project,Settings.Project_Directory_Name);
SaveDataDir = fullfile(Settings.home,'data',DataSetName);
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
Renew_all_data = true;

% Load previously collected data
if isfile(SaveDataDir) && ~Renew_all_data
    D = load(SaveDataDir,'Data');
    Data = D.Data;
else
    Data = struct();
end
Warn_summary = '';
Missing_data = 0;

% List all basis experiment directories
disp('Beginning Search Through Raw Data...')
Experiment = Load_Experimental_Data;

for Salt_idx = 1:length(Salts)
    Salt = Salts{Salt_idx};
    SaltDataDir = fullfile(ProjectDir,Salt);
    disp('************************************')
    disp(['Current Salt: ' Salt])
    disp('************************************')
    Warn_summary = [Warn_summary newline '************************************' ...
        newline Salt newline '************************************']; %#ok<*AGROW>
    
    % List all inner directories
    d_jobs = dir(SaltDataDir);
    d_jobs = d_jobs([d_jobs(:).isdir]);
    d_jobs = d_jobs(~ismember({d_jobs(:).name},{'.','..'}));
    
    if isfield(Data,Salt)
        prev_jobs = get_all_MP_job_names(Data.(Salt));
    else
        prev_jobs = {};
    end
    
    [~,d_jobs_idx] = setdiff(strrep({d_jobs.name},'-','_'),prev_jobs);
    d_jobs = d_jobs(d_jobs_idx);

    for job_idx = 1:length(d_jobs)
        JobName = d_jobs(job_idx).name;
        if ~contains(JobName,DataKeyword)
            continue
        end
        JobDir = fullfile(d_jobs(job_idx).folder,JobName);
        
        Theory_info = regexp(JobName,'(.+?)_Alexandria','tokens','once');
        Theory = Theory_info{1};
        disp(['Current Job: ' Salt ' ' Theory])
        
        % Search for MPResults and Density profile files 
        ResultsFile = fullfile(JobDir,[JobName '_MPResults.mat']);
        InputFile = fullfile(JobDir,[JobName '.inp']);
        DensityProfileFile = fullfile(JobDir,[JobName '_Density_Profile.xvg']);
        
        if ~isfile(ResultsFile)
            Warn_summary = [Warn_summary newline Salt ' ' JobName ': Job not complete (no results file found).'];
            Missing_data = Missing_data + 1;
            continue
        elseif ~isfile(InputFile)
            Warn_summary = [Warn_summary newline Salt ' ' JobName ': Missing input file.'];
            Missing_data = Missing_data + 1;
        end
        
        try
            dat = load(ResultsFile,'T_dat');
            T_dat = dat.T_dat;
        catch
            Warn_summary = [Warn_summary newline Salt ' ' JobName ': Cannot load results file, check for file corruption.'];
            continue
        end
        
        try
            dat = load(InputFile,'Settings','-mat');
            Settings = dat.Settings;
            Structure = Settings.Structure;
        catch
            Warn_summary = [Warn_summary newline Salt ' ' JobName ': Cannot load input file, check for file corruption.'];
            continue
        end
        
        Tm = T_dat.T_Trace(~T_dat.Freeze_Trace & ~T_dat.Melt_Trace);
        MP_confirmed = true;
        if length(Tm) > 1
            Tms = sort(Tm);
            idx = floor(length(Tms)/2);
            Tm = Tms(idx);
            
        elseif isempty(Tm)
            % This is reached when the algorithm bracketed the MP less than
            % the desired range, but did not find a simulation that
            % maintained solid-liquid coexistence.
            if isfolder(fullfile(JobDir,['Tm_' num2str(T_dat.T,'%.4f')]))
                Tm = T_dat.T; % Check for old convention (slower) which was the midpoint
            else
                Tm = T_dat.dT(1); % By new convention we use lower bound as "estimated melting point"
            end
            
            MP_confirmed = false;
        end
        Tm_dir =  fullfile(JobDir,['Tm_' num2str(Tm,'%.4f')]);
        
        % Gather some basic calculation data
        if ~isfolder(Tm_dir) && sum(~T_dat.Freeze_Trace & ~T_dat.Melt_Trace) > 1
            Tmold = T_dat.T_Trace(~T_dat.Freeze_Trace & ~T_dat.Melt_Trace);
            Tmold = Tmold(end);
            Tm_dir2 =  fullfile(JobDir,['Tm_' num2str(Tmold,'%.4f')]);
            
            if ~isfolder(Tm_dir2)
                Warn_summary = [Warn_summary newline Salt ' ' JobName ': Cannot find folder: Tm_' num2str(Tm,'%.4f')];
                continue
            else
                Tm = Tmold;
                Tm_dir = Tm_dir2;
            end
        elseif ~isfolder(Tm_dir)
            Warn_summary = [Warn_summary newline Salt ' ' JobName ': Cannot find folder: Tm_' num2str(Tm,'%.4f')];
            continue
        end
        
        % Find the final geometry
        d_tm = dir(fullfile(Tm_dir,'*.gro'));
        Final_geom_file = d_tm(contains({d_tm.name},'OutConf'));

        if length(Final_geom_file) > 1
            % If multiple gro files exist, choose most recent
            Times = zeros(1,length(Final_geom_file));
            for y = 1:length(Final_geom_file)
                Times(y) = datenum(Final_geom_file(y).date,...
                    'dd-mmm-yyyy HH:MM:SS');
            end
            [~,Idx] = max(Times);
            Final_geom_file = Final_geom_file(Idx);
        end

        if isempty(Final_geom_file)
            Warn_summary = [Warn_summary newline Salt ' ' JobName ': No output gro file found in subfolder: Tm_' num2str(Tm,'%.4f')];
            continue
        end
        
        try
            Final_geom = load_gro_file(fullfile(Final_geom_file.folder,Final_geom_file.name));
        catch
            Warn_summary = [Warn_summary newline Salt ' ' JobName ': Problem loading final gro file, check for file corruption.'];
            continue
        end
        N_atoms_total = Final_geom.N_atoms;
        N_atoms_liquid = ceil((N_atoms_total/2)*Settings.Liquid_Fraction)*2;
        N_atoms_solid = N_atoms_total - N_atoms_liquid;

        % Create a data structure and add the data
        Data.(Salt).(Theory).(Structure) = T_dat;
        Data.(Salt).(Theory).(Structure).N_Total = N_atoms_total;
        Data.(Salt).(Theory).(Structure).N_Solid = N_atoms_solid;
        Data.(Salt).(Theory).(Structure).N_Liquid = N_atoms_liquid;
        Data.(Salt).(Theory).(Structure).a = norm(Final_geom.a_vec)*10; % Angstroms, final size of box in a-direction
        Data.(Salt).(Theory).(Structure).b = norm(Final_geom.b_vec)*10; % Angstroms
        Data.(Salt).(Theory).(Structure).c = norm(Final_geom.c_vec)*10; % Angstroms
        Data.(Salt).(Theory).(Structure).Tm = Tm;
        Data.(Salt).(Theory).(Structure).IsBubble = Settings.GenCluster;
        Data.(Salt).(Theory).(Structure).MP_confirmed = MP_confirmed;
        Data.(Salt).(Theory).(Structure).JobName = JobName;
        
        % Difference between experimental and calculated melting points
        Data.(Salt).(Theory).(Structure).dT_Exp = Tm - Experiment.(Salt).mp;
        
        % Load in the density profile if it exists
        if isfile(DensityProfileFile)
            dat = import_xvg(DensityProfileFile);
            Data.(Salt).(Theory).(Structure).DensityProfile = dat;
        else
            Warn_summary = [Warn_summary newline Salt ' ' JobName ': No density profile found.'];
            Data.(Salt).(Theory).(Structure).DensityProfile = [];
        end
        
        % If the calculation is a interface, find the initial length of the solid portion of the box
        Sol_info_file = fullfile(Tm_dir,[JobName '_' num2str(Tm,'%.4f') '_SolInfo.mat']);
        if isfile(Sol_info_file)
            dat = load(Sol_info_file,'Solid_Geometry');
            Solid_Geometry = dat.Solid_Geometry;
            Data.(Salt).(Theory).(Structure).a_vec_sol = Solid_Geometry.a_vec;
            Data.(Salt).(Theory).(Structure).b_vec_sol = Solid_Geometry.b_vec;
            Data.(Salt).(Theory).(Structure).c_vec_sol = Solid_Geometry.c_vec;
            Data.(Salt).(Theory).(Structure).c_vec_tot = Solid_Geometry.c_vec_tot;
            Data.(Salt).(Theory).(Structure).Liq_L     = Solid_Geometry.L;
            if isfield(Solid_Geometry,'Liquid_Vol')
                Data.(Salt).(Theory).(Structure).Liquid_Vol= Solid_Geometry.Liquid_Vol;
                Data.(Salt).(Theory).(Structure).Solid_Vol = Solid_Geometry.Solid_Vol;
                Data.(Salt).(Theory).(Structure).Solid_R   = Solid_Geometry.Solid_R;
                Data.(Salt).(Theory).(Structure).Total_Vol = Solid_Geometry.Solid_Vol + Solid_Geometry.Liquid_Vol;
            else
                V_Liquid = Cell_Volume(Solid_Geometry.a_vec,Solid_Geometry.b_vec,Solid_Geometry.c_vec_tot - Solid_Geometry.c_vec);
                V_Solid = Cell_Volume(Solid_Geometry.a_vec,Solid_Geometry.b_vec,Solid_Geometry.c_vec);
                Data.(Salt).(Theory).(Structure).Liquid_Vol= V_Liquid;
                Data.(Salt).(Theory).(Structure).Solid_Vol = V_Solid;
                Data.(Salt).(Theory).(Structure).Solid_R   = ( 3*V_Solid/(4*pi) )^(1/3);
                Data.(Salt).(Theory).(Structure).Total_Vol = V_Liquid + V_Solid;
            end
        else
            MinJobInpFile = fullfile(JobDir,'Minimization','TempJobInfo.mat');
            if isfile(MinJobInpFile)
                dat = load(MinJobInpFile);
                
                Molar_Volume = Cell_Volume(dat.Geometry.Transform(1,:)*dat.Geometry.a,...
                    dat.Geometry.Transform(2,:)*dat.Geometry.b,...
                    dat.Geometry.Transform(3,:)*dat.Geometry.c)/dat.Geometry.N; % Ang^3/atom
                Solid_Vol = Molar_Volume*N_atoms_solid/(10^3); % nm^3
                Solid_R = ( 3*Solid_Vol/(4*pi) )^(1/3); % nm
                Total_Vol = Cell_Volume(Final_geom.a_vec,Final_geom.b_vec,Final_geom.c_vec); % nm^3
                
                Data.(Salt).(Theory).(Structure).a_vec_sol = dat.a_vec; % nm
                Data.(Salt).(Theory).(Structure).b_vec_sol = dat.b_vec; % nm
                Data.(Salt).(Theory).(Structure).c_vec_sol = dat.c_vec; % nm
                Data.(Salt).(Theory).(Structure).c_vec_tot = Final_geom.c_vec;
                Data.(Salt).(Theory).(Structure).Liq_L     = norm(Final_geom.c_vec) - norm(dat.c_vec);
                Data.(Salt).(Theory).(Structure).Liquid_Vol= Total_Vol - Solid_Vol;
                Data.(Salt).(Theory).(Structure).Solid_Vol = Solid_Vol;
                Data.(Salt).(Theory).(Structure).Solid_R   = Solid_R;   % radius of solid (if it is a cluster)
                Data.(Salt).(Theory).(Structure).Total_Vol = Total_Vol; % nm^3
            else
                Data.(Salt).(Theory).(Structure).a_vec_sol = nan; % nm
                Data.(Salt).(Theory).(Structure).b_vec_sol = nan; % nm
                Data.(Salt).(Theory).(Structure).c_vec_sol = nan; % nm
                Data.(Salt).(Theory).(Structure).c_vec_tot = nan;
                Data.(Salt).(Theory).(Structure).Liq_L     = nan;
                Data.(Salt).(Theory).(Structure).Liquid_Vol= nan;
                Data.(Salt).(Theory).(Structure).Solid_Vol = nan;
                Data.(Salt).(Theory).(Structure).Solid_R   = nan;
                Data.(Salt).(Theory).(Structure).Total_Vol = nan;
            end
        end
        if Settings.GenCluster
            Data.(Salt).(Theory).(Structure).c_over_a = 1;
        else
            Data.(Salt).(Theory).(Structure).c_over_a = Settings.c_over_a;
        end
        
        Tm_dt = T_dat.dT;
        dT_Exp = Data.(Salt).(Theory).(Structure).dT_Exp;
        if length(Tm_dt) == 1
            if Tm_dt < Tm
                Lower_Bound = Tm_dt;
                Upper_Bound = [];
            else
                Lower_Bound = [];
                Upper_Bound = Tm_dt;
            end
        elseif isempty(Tm_dt)
            Warn_summary = [Warn_summary newline Salt ' ' JobName ': No bounds on melting point.'];
            Lower_Bound = [];
            Upper_Bound = [];
        else
            Lower_Bound = Tm_dt(1);
            Upper_Bound = Tm_dt(2);
        end
    end
        
    Warn_summary = [Warn_summary newline '************************************'];
end
disp('Finished collecting Data.')


disp('Finished Calculaing Differences in Lattice Energy Between Structures.')
disp('Saving Data.')

% Save the dataset
save(SaveDataDir,'Data')
disp('Finished Saving Data! Job Complete.')
disp('************************************')
disp('Summary of missing data.')
disp('************************************')
disp(['Number of missing data points: ' num2str(Missing_data)])
disp(Warn_summary)
