% Collect_MP_Data

%% Initialize directories of input and output data
Settings = Initialize_MD_Settings;
Settings.Project_Directory_Name = 'Melting_Point_Studies';
DataSetName = 'Melting_Point_Data.mat';
DataKeyword = '';
ProjectDir = fullfile(Settings.project,Settings.Project_Directory_Name);
SaveDataDir = fullfile(Settings.home,'data',DataSetName);
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
Renew_all_data = false;

fileID = fopen(fullfile(ProjectDir,'Melting_Point_Data.csv'),'wt');
fprintf(fileID, [strjoin(repmat({'%s'},1,14),',') '\n'],...
    'Salt',...
    'Theory',...
    'Model',...
    'Job Name',...
    'Tm [K]',...
    'Tm Lower Bound [K]',...
    'Tm Upper Bound [K]',...
    'dTm Experiment [K]',...
    'N Total',...
    'N Liquid',...
    'N Solid',...
    'Box length a [A]',...
    'Box length b [A]',...
    'Box length c [A]',...    
    'Solid a [A]',...
    'Solid b [A]',...
    'Solid c [A]');

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

    % Put experimental data on the first line
    Out = {Salt,...
          'Experiment',...
          '',...
          '',...
          Experiment.(Salt).mp,...
          Experiment.(Salt).mp-Experiment.(Salt).dmp,...
          Experiment.(Salt).mp+Experiment.(Salt).dmp,...
          0,...
          '',...
          '',...
          '',...
          '',...
          '',...
          '',...
          '',...
          '',...
          ''};
    fprintf(fileID, '%s,%s,%s,%s,%f,%f,%f,%f,%s,%s,%s,%s,%s,%s\n',Out{:});
    
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
        
        Theory_info = regexp(JobName,'(.+?)_[A-Z]_(JC|JC3P|JC4P|JCSD|TF|BH)(_Model_.+?)*_(NPT|NVT|NVE)','tokens','once');
        JobNametxt = ['MP_' strrep(strrep(Theory_info{1},'-','_'),'.','_')];
        Theory = Theory_info{2};
        if ~isempty(Theory_info{3})
            Model = regexp(Theory_info{3},'_Model_(.+)','tokens','once');
            Model = Model{1};
            Modeltxt = [Theory '_Model_' Model];
        else
            Model = '';
            Modeltxt = Theory;
        end
        disp(['Current Job: ' Modeltxt ' / ' strrep(JobNametxt,'MP_','')])
        
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
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt) = T_dat;
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt).N_Total = N_atoms_total;
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt).N_Solid = N_atoms_solid;
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt).N_Liquid = N_atoms_liquid;
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt).a = norm(Final_geom.a_vec)*10; % Angstroms, final size of box in a-direction
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt).b = norm(Final_geom.b_vec)*10; % Angstroms
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt).c = norm(Final_geom.c_vec)*10; % Angstroms
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Tm = Tm;
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt).IsBubble = Settings.GenCluster;
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt).MP_confirmed = MP_confirmed;
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt).JobName = JobName;
        
        % Difference between experimental and calculated melting points
        Data.(Salt).(Modeltxt).(Structure).(JobNametxt).dT_Exp = Tm - Experiment.(Salt).mp;
        
        % Load in the density profile if it exists
        if isfile(DensityProfileFile)
            dat = import_xvg(DensityProfileFile);
            Data.(Salt).(Modeltxt).(Structure).(JobNametxt).DensityProfile = dat;
        else
            Warn_summary = [Warn_summary newline Salt ' ' JobName ': No density profile found.'];
            Data.(Salt).(Modeltxt).(Structure).(JobNametxt).DensityProfile = [];
        end
        
        % If the calculation is a interface, find the initial length of the solid portion of the box
        Sol_info_file = fullfile(Tm_dir,[JobName '_' num2str(Tm,'%.4f') '_SolInfo.mat']);
        if isfile(Sol_info_file)
            dat = load(Sol_info_file,'Solid_Geometry');
            Solid_Geometry = dat.Solid_Geometry;
            Data.(Salt).(Modeltxt).(Structure).(JobNametxt).a_vec_sol = Solid_Geometry.a_vec;
            Data.(Salt).(Modeltxt).(Structure).(JobNametxt).b_vec_sol = Solid_Geometry.b_vec;
            Data.(Salt).(Modeltxt).(Structure).(JobNametxt).c_vec_sol = Solid_Geometry.c_vec;
            Data.(Salt).(Modeltxt).(Structure).(JobNametxt).c_vec_tot = Solid_Geometry.c_vec_tot;
            Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Liq_L     = Solid_Geometry.L;
            if isfield(Solid_Geometry,'Liquid_Vol')
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Liquid_Vol= Solid_Geometry.Liquid_Vol;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Solid_Vol = Solid_Geometry.Solid_Vol;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Solid_R   = Solid_Geometry.Solid_R;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Total_Vol = Solid_Geometry.Solid_Vol + Solid_Geometry.Liquid_Vol;
            else
                V_Liquid = Cell_Volume(Solid_Geometry.a_vec,Solid_Geometry.b_vec,Solid_Geometry.c_vec_tot - Solid_Geometry.c_vec);
                V_Solid = Cell_Volume(Solid_Geometry.a_vec,Solid_Geometry.b_vec,Solid_Geometry.c_vec);
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Liquid_Vol= V_Liquid;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Solid_Vol = V_Solid;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Solid_R   = ( 3*V_Solid/(4*pi) )^(1/3);
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Total_Vol = V_Liquid + V_Solid;
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
                
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).a_vec_sol = dat.a_vec; % nm
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).b_vec_sol = dat.b_vec; % nm
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).c_vec_sol = dat.c_vec; % nm
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).c_vec_tot = Final_geom.c_vec;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Liq_L     = norm(Final_geom.c_vec) - norm(dat.c_vec);
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Liquid_Vol= Total_Vol - Solid_Vol;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Solid_Vol = Solid_Vol;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Solid_R   = Solid_R;   % radius of solid (if it is a cluster)
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Total_Vol = Total_Vol; % nm^3
            else
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).a_vec_sol = nan; % nm
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).b_vec_sol = nan; % nm
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).c_vec_sol = nan; % nm
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).c_vec_tot = nan;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Liq_L     = nan;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Liquid_Vol= nan;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Solid_Vol = nan;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Solid_R   = nan;
                Data.(Salt).(Modeltxt).(Structure).(JobNametxt).Total_Vol = nan;
            end
        end
        if Settings.GenCluster
            Data.(Salt).(Modeltxt).(Structure).(JobNametxt).c_over_a = 1;
        else
            Data.(Salt).(Modeltxt).(Structure).(JobNametxt).c_over_a = Settings.c_over_a;
        end
        
        Tm_dt = T_dat.dT;
        dT_Exp = Data.(Salt).(Modeltxt).(Structure).(JobNametxt).dT_Exp;
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
        
        Out = {Salt,...
              Theory,...
              Model,...
              JobName,...
              Tm,...
              Lower_Bound,...
              Upper_Bound,...
              dT_Exp,...
              N_atoms_total,...
              N_atoms_liquid,...
              N_atoms_solid,...
              Data.(Salt).(Modeltxt).(Structure).(JobNametxt).a,...
              Data.(Salt).(Modeltxt).(Structure).(JobNametxt).b,...
              Data.(Salt).(Modeltxt).(Structure).(JobNametxt).c,...
              norm(Data.(Salt).(Modeltxt).(Structure).(JobNametxt).a_vec_sol),...
              norm(Data.(Salt).(Modeltxt).(Structure).(JobNametxt).b_vec_sol),...
              norm(Data.(Salt).(Modeltxt).(Structure).(JobNametxt).c_vec_sol)};
        fprintf(fileID, '%s,%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',Out{:});
    end
        
    Warn_summary = [Warn_summary newline '************************************'];
end
disp('Finished collecting Data.')


disp('Finished Calculaing Differences in Lattice Energy Between Structures.')
disp('Saving Data.')
fclose(fileID);

% Save the dataset
save(SaveDataDir,'Data')
disp('Finished Saving Data! Job Complete.')
disp('************************************')
disp('Summary of missing data.')
disp('************************************')
disp(['Number of missing data points: ' num2str(Missing_data)])
disp(Warn_summary)
