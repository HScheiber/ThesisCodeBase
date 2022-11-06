function [feval,fderiv,User_data] = Melting_Point_Check(T,Settings)

    % Load the data trace
    T_dat = load(Settings.CurrentTFile,'T_dat').T_dat;
    if ~isfield(T_dat,'Alt_Structure')
        T_dat.Alt_Structure = false;
    end
    
    if istable(T)
        T = T{1,1};
    end
    
    % Check if the current point has already been calculated (i.e. during a restart)
    if any(abs(T_dat.T_Trace - T) < sqrt(eps))
        [~,idx] = min(abs(T_dat.T_Trace - T));
        f = T_dat.f_Trace(idx);
        df = T_dat.df_Trace(idx);
        if f <= 0 && Settings.Continue
            % Continue a calculation from the previously determined melting point with new updated settings
        else
            feval = f; % function evaluation
            fderiv = df; % function derivative
            User_data = T_dat; % user data
            return
        end
    end
    
    if Settings.Verbose
        disp(newline)
        disp(repmat('*',1,40));
        disp(['Current Temperature: ' num2str(T,'%.4f') ' K']);
        disp(repmat('*',1,40));
        disp(newline)
    end
    
    % Make the new work directory
    RefGeomDir = Settings.RefGeomDir;
    OuterDir = Settings.WorkDir;
    WorkDir = fullfile(Settings.WorkDir,['T_' num2str(T,'%.4f')]);
    if ~exist(WorkDir,'dir')
        mkdir(WorkDir)
    end
    
    % Dynamically set the CheckTime to the closest divisor of MaxCheckTime for better efficiency
    CheckTime_Trace = 1./T_dat.f_Trace; % ps
    T_diffs = T - T_dat.T_Trace; % K
    
    T_diffs_below = T_diffs(T_diffs < 0);
    CheckTime_below = CheckTime_Trace(T_diffs < 0);
    [T_closest_below,bel_idx] = max(T_diffs_below);
    CheckTime_closest_below = CheckTime_below(bel_idx);
    if CheckTime_closest_below < 0
        CheckTime_closest_below = Settings.MaxCheckTime/2;
    end
    
    T_diffs_above = T_diffs(T_diffs > 0);
    CheckTime_above = CheckTime_Trace(T_diffs > 0);
    [T_closest_above,abv_idx] = min(T_diffs_above);
    CheckTime_closest_above = CheckTime_above(abv_idx);
    if CheckTime_closest_above < 0
        CheckTime_closest_above = Settings.MaxCheckTime/2;
    end
    
    if isempty(CheckTime_closest_below) || isempty(CheckTime_closest_above)
        CheckTime = Settings.CheckTime; % ps
    else
        w_bel = 1 - abs(T_closest_below)/(abs(T_closest_below) + abs(T_closest_above));
        w_abv = 1 - abs(T_closest_above)/(abs(T_closest_below) + abs(T_closest_above));
        CheckTime_Guess = w_bel*CheckTime_closest_below + w_abv*CheckTime_closest_above; % ps
        
        MaxCheckTimeDivisors = alldivisors(Settings.MaxCheckTime);
        [~,idx] = min(abs(MaxCheckTimeDivisors - CheckTime_Guess));
        CheckTime = max(MaxCheckTimeDivisors(idx),Settings.CheckTime); % ps
    end
    
    % Check if the current temperature is too far from the reference temperature used to built the initial conditions
    Strucure_Ref_File = fullfile(RefGeomDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '.' Settings.CoordType]);
    Sol_Ref_File = fullfile(RefGeomDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '_SolInfo.mat']);
    if abs(T_dat.T_ref - T) >= Settings.MaxTDiff || ~isfile(Strucure_Ref_File) || ~isfile(Sol_Ref_File)
        
        % Search for previous initial conditions with a suitable temperature
        MakeNewConfig = true;
        if any(abs(T_dat.Ref_Density_Trace - T) < Settings.MaxTDiff)
            [~,idx] = min(abs(T_dat.Ref_Density_Trace - T));
            T_dat.T_ref = T_dat.Ref_Density_Trace(idx);
            Strucure_Ref_File = fullfile(RefGeomDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '.' Settings.CoordType]);
            Sol_Ref_File = fullfile(RefGeomDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '_SolInfo.mat']);
            Strucure_In_File = fullfile(WorkDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '.' Settings.CoordType]);
            Sol_In_File = fullfile(WorkDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '_SolInfo.mat']);
            
            if isfile(Strucure_Ref_File) && isfile(Sol_Ref_File)
                MakeNewConfig = false;
            end
        end
        
        if MakeNewConfig
            if Settings.Verbose
                disp('Updating initial system configuration...')
            end
            % Make a new reference T
            T_dat.T_ref = T;
            Settings.Target_T = T_dat.T_ref; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            Settings.MDP.Initial_T = T_dat.T_ref; % Initial termpature at which to generate velocities

            % Reset the initial configuration
            Output = Setup_LiX_Simulation(Settings);
            
            % Check if the output is aborted
            if Output.Aborted
                if Output.SolidMelted
                    
                    time_to_phase_change = 1;
                    f = 1/time_to_phase_change;
                    df = 5000/(time_to_phase_change);
                    T_dat.T_Melt_Trace = [T_dat.T_Freeze_Trace T];
                    
                    % Update the current T_mp lower error bound
                    T_dat.dT(2) = T;
                    T_dat.df_bracket(2) = df;
                    T_dat.Ref_Density_Trace = [T_dat.Ref_Density_Trace T_dat.T_ref];
                    T_dat.T = T;
                    T_dat.T_Trace = [T_dat.T_Trace T];
                    T_dat.f_Trace = [T_dat.f_Trace f];
                    T_dat.df_Trace = [T_dat.df_Trace df];
                    T_dat.Freeze_Trace = [T_dat.Freeze_Trace false];
                    T_dat.Melt_Trace = [T_dat.Melt_Trace true];
                    T_dat.Alt_Structure = false;
                    
                elseif Output.LiquidFroze || Output.LiquidAmorphous

                    time_to_phase_change = 1;
                    f = 1/time_to_phase_change;
                    df = -5000/(time_to_phase_change);
                    T_dat.T_Freeze_Trace = [T_dat.T_Freeze_Trace T];
                    
                    % Update the current T_mp lower error bound
                    T_dat.dT(1) = T;
                    T_dat.df_bracket(1) = df;
                    T_dat.Ref_Density_Trace = [T_dat.Ref_Density_Trace T_dat.T_ref];
                    T_dat.T = T;
                    T_dat.T_Trace = [T_dat.T_Trace T];
                    T_dat.f_Trace = [T_dat.f_Trace f];
                    T_dat.df_Trace = [T_dat.df_Trace df];
                    T_dat.Freeze_Trace = [T_dat.Freeze_Trace true];
                    T_dat.Melt_Trace = [T_dat.Melt_Trace false];
                    T_dat.Alt_Structure = false;
                    
                else
                    % Solid or liquid crystallized to an unwanted structure
                    % homogeneously, or calculation failed due to bad potential
                    % In this case, abort the calculation...
                    f = -1;
                    df = 0;
                    T_dat.Alt_Structure = true;
                    T_dat.T = T;
                    T_dat.T_Trace = [T_dat.T_Trace T];
                    T_dat.f_Trace = [T_dat.f_Trace f];
                    T_dat.df_Trace = [T_dat.df_Trace df];
                    T_dat.Freeze_Trace = [T_dat.Freeze_Trace false];
                    T_dat.Melt_Trace = [T_dat.Melt_Trace false];
                end
                
                while true
                    try
                        copyfile(Settings.CurrentTFile,Settings.PrevTFile)
                        save(Settings.CurrentTFile,'T_dat')
                        break
                    catch
                        pause(10)
                    end
                end
                feval = f; % function evaluation
                fderiv = df; % function derivative
                User_data = T_dat; % user data
                return
            else
                % Rename the gro and sol.mat file at the given temperature
                Strucure_Ref_File_old = fullfile(OuterDir,[Settings.JobName '.' Settings.CoordType]);
                Strucure_Ref_File = fullfile(RefGeomDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '.' Settings.CoordType]);
                movefile(Strucure_Ref_File_old, Strucure_Ref_File)
                % Move the partial files from minimization
                Sol_Ref_File_old = fullfile(OuterDir,[Settings.JobName '_SolInfo.mat']);
                Sol_Ref_File = fullfile(RefGeomDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '_SolInfo.mat']);
                movefile(Sol_Ref_File_old, Sol_Ref_File)
                Strucure_In_File = fullfile(WorkDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '.' Settings.CoordType]);
                Sol_In_File = fullfile(WorkDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '_SolInfo.mat']);

                % Update reference T file
                T_dat.Ref_Density_Trace = [T_dat.Ref_Density_Trace T_dat.T_ref];
                while true
                    try
                        copyfile(Settings.CurrentTFile,Settings.PrevTFile)
                        save(Settings.CurrentTFile,'T_dat')
                        break
                    catch
                        pause(10)
                    end
                end
                if Settings.Verbose
                    disp('Successfully updated initial configuration.')
                end
            end
            
        elseif Settings.Manual_Box
            Settings.SuperCellFile = Strucure_Ref_File;
            UpdateTopology(Settings)
        end
    else
        Strucure_In_File = fullfile(WorkDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '.' Settings.CoordType]);
        Sol_In_File = fullfile(WorkDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '_SolInfo.mat']);
        
        if Settings.Manual_Box
            Settings.SuperCellFile = Strucure_Ref_File;
            UpdateTopology(Settings)
        end
    end
    if Settings.Verbose
        disp(['Intial configuration reference density: T = ' num2str(T_dat.T_ref,'%.4f') ' K.']);
    end
    
    % Look for a checkpoint file
    cpt_check = dir(fullfile(WorkDir,[Settings.JobName '_*.cpt']));
    
    % Delete any corrupted checkpoint files
    for idx = length(cpt_check):-1:1
        cptch = fullfile(WorkDir,cpt_check(idx).name);
        [errc,~] = system([Settings.gmx_loc Settings.g_check ' -f ' windows2unix(cptch)]);
        if errc ~= 0
            if Settings.Verbose
                disp(['Detected and deleted corrupted checkpoint file: ' cptch])
            end
            delete(cptch)
            prevcpt = strrep(cptch,'.cpt','_prev.cpt');
            if isfile(prevcpt)
                movefile(prevcpt,cptch)
            end
            cpt_check = dir(fullfile(WorkDir,[Settings.JobName '_*.cpt']));
        end
    end
    
    if length(cpt_check) > 1
        % Remove any "prev" checkpoints
        if sum(~contains({cpt_check.name},'prev')) > 0
            cpt_check = cpt_check(~contains({cpt_check.name},'prev'));
        end
        [~,idx] = max([cpt_check(:).datenum]);
        cpt_check = cpt_check(idx);
    end
    
    % Checkpoint file found. Ensure it is usable.
    if ~isempty(cpt_check)
        ContinueFromCheckPoint = true;
        
        % Grab the current checkpoint number (text in format NNN)
        cpt_number = regexp(cpt_check.name,'_([0-9]{3})(_step[0-9]+)*(_prev)*.cpt','once','tokens');
        
        % Grab the step size from the mdp file and update the CheckTime to match the previous one
        MDP_in_File = fullfile(WorkDir,[Settings.JobName '.mdp']);
        MDP_txt = fileread(MDP_in_File);
        nstepsre = regexp(MDP_txt,'nsteps *= *(.+?) *;','tokens','once');
        stepsizere = regexp(MDP_txt,'dt *= *(.+?) *;','tokens','once');
        PrevCheckTime = str2double(nstepsre{1})*str2double(stepsizere{1});
        
        if PrevCheckTime < Settings.CheckTime
            ContinueFromCheckPoint = false;  % this should not be reached unless the settings are changed
        end
        
        Trajectory_File = fullfile(WorkDir,[Settings.JobName '.trr']);
        gmx_check_cmd = [Settings.gmx_loc Settings.g_check ' -f ' windows2unix(Trajectory_File)];
        [errcode,outchk] = system(gmx_check_cmd);
        if errcode ~= 0
            ContinueFromCheckPoint = false;
        end
        
        ftime = regexp(outchk,'Coords +([0-9]|\.|e|E)+ +([0-9]|\.|e|E)+','tokens','once');
        if isempty(ftime)
            ContinueFromCheckPoint = false;
        else
            telpse = (str2double(ftime{1})-1)*str2double(ftime{2}); % Elapsed time in ps
            if telpse < Settings.CheckTime
                ContinueFromCheckPoint = false;
            end
        end
    else
        ContinueFromCheckPoint = false;
    end
    
    % If usable checkpoint file found...
    if ContinueFromCheckPoint
        if Settings.Verbose
            disp(['Usable checkpoint file found for T = ' num2str(T,'%.4f')])
            disp(['Previously, simulation stopped at ' num2str(telpse,'%.1f') ' ps'])
        end
        CheckTime = max(PrevCheckTime,Settings.CheckTime);
        
        % Maximum number of time segments
        max_steps = num2str(Settings.MaxCheckTime/CheckTime,'%03.f');
        
        % Define some gmx files
        Log_File = fullfile(WorkDir,[Settings.JobName '.log']);
        Energy_file = fullfile(WorkDir,[Settings.JobName '.edr']);
        Trajectory_File = fullfile(WorkDir,[Settings.JobName '.trr']);
        Structure_Out_File = fullfile(WorkDir,[Settings.JobName '_OutConf_' cpt_number{1} '.' Settings.CoordType]);
        CheckPoint_File = fullfile(WorkDir,cpt_check.name);
        Traj_Conf_File = fullfile(WorkDir,[Settings.JobName '_' cpt_number{1} '.tpr']);
        Table_File = fullfile(OuterDir,[Settings.JobName '_Table.xvg']);
        
        % Check if the structure OUT file exists, indicating the previous segment finished
        if isfile(Structure_Out_File)
            if Settings.Verbose
                disp(['Step ' cpt_number{1} '/' max_steps ' previously completed.'])
            end
        else
            % If not, continue the current segment of the simulation
            mdrun_command = [Settings.gmx Settings.mdrun ' -s ' windows2unix(Traj_Conf_File) ...
                ' -o ' windows2unix(Trajectory_File) ' -g ' windows2unix(Log_File) ...
                ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Structure_Out_File) ...
                ' -cpi ' windows2unix(CheckPoint_File) ' -cpo ' windows2unix(CheckPoint_File) ...
                Settings.mdrun_opts];
            
            if isfile(Table_File) % If tabulated potential required
                mdrun_command = [mdrun_command ' -table ' windows2unix(Table_File)];
            end

            MDtimer = tic;
            if Settings.Verbose
                disp(['Continuing MD simulation segment ' cpt_number{1} '/' max_steps ...
                    ' (' num2str(str2double(cpt_number{1})*CheckTime,'%.0f') '/' num2str(Settings.MaxCheckTime,'%.0f') ' ps).'])
            end
            
            if Settings.slurm
                mdrun_command = regexprep(mdrun_command,' -maxh ([0-9]|\.)+','','once');
                HoursRemain = JobTimeRemaining;
                mdrun_command = [mdrun_command ' -maxh ' num2str(HoursRemain)];
            else
                mdrun_command = regexprep(mdrun_command,' -maxh ([0-9]|\.)+','','once');
            end
            
            [errcode,~] = system(mdrun_command);
            if errcode ~= 0
                if Settings.Verbose
                    disp(['MD simulation segment ' cpt_number{1} '/' max_steps ...
                        ' (' num2str(str2double(cpt_number{1})*CheckTime,'%.0f') '/' num2str(Settings.MaxCheckTime,'%.0f') ...
                        ' ps) failed to complete. Time elapsed: ' datestr(seconds(toc(MDtimer)),'HH:MM:SS')])
                end
            elseif ~isfile(Structure_Out_File)
                if Settings.Verbose
                    disp('Calculation time is almost up, ending MATLAB session now.')
                end
                exit
            end
            if Settings.Verbose && errcode == 0
                disp(['MD simulation segment ' cpt_number{1} '/' max_steps ...
                    ' (' num2str(CheckTime,'%.0f') '/' num2str(Settings.MaxCheckTime,'%.0f') ...
                    ' ps) complete. Time elapsed: ' datestr(seconds(toc(MDtimer)),'HH:MM:SS')])
            end
        end
        
        if Settings.Verbose
            disp('Checking melting/freezing status...')
        end
        CheckStructureTimer = tic;
        try
            % Check volume has not exploded
            VolCheck_File = fullfile(WorkDir,[Settings.JobName '_001_VolCheck.xvg']);
            VolCheck_Log_File = fullfile(WorkDir,[Settings.JobName '_001_VolCheck.log']);
            
            % Check energy options
            gmx_command = [Settings.wsl 'echo "0" ' Settings.pipe ...
                ' ' strrep(Settings.gmx_loc,Settings.wsl,'') Settings.g_energy ...
                ' -f ' windows2unix(Energy_file) ' -s ' windows2unix(Traj_Conf_File)];
            [~,outpt] = system(gmx_command);
            en_opts = regexp(outpt,'-+\n.+','match','once');
            En_set = '';
            En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Volume','tokens','once'))];
            En_set = [En_set ' 0'];
            En_set = regexprep(En_set,' +',' ');
            
            % Grab data from results
            VolCheck_command = [Settings.wsl 'echo ' En_set ' ' Settings.pipe ...
                ' ' strrep(Settings.gmx_loc,Settings.wsl,'') Settings.g_energy ...
                ' -f ' windows2unix(Energy_file) ' -o ' windows2unix(VolCheck_File) ...
                ' -s ' windows2unix(Traj_Conf_File) ' ' Settings.passlog windows2unix(VolCheck_Log_File)];
            [err,~] = system(VolCheck_command);
            if err ~= 0
                error('Failed to collect energy data.')
            end
            
            VolCheck_Data = import_xvg(VolCheck_File);
            dev_from_init = VolCheck_Data(:,2)./VolCheck_Data(1,2);
            max_dev = max(abs(dev_from_init));
            
            if max_dev >= 1.5
                if Settings.Verbose
                    disp(['System volume exploded (' num2str(max_dev) ...
                        ' deviation in volume detected). Aborting calculation.'])
                end
                
                f = -1;
                df = 0;
                T_dat.Alt_Structure = true;
                T_dat.T = T;
                T_dat.T_Trace = [T_dat.T_Trace T];
                T_dat.f_Trace = [T_dat.f_Trace f];
                T_dat.df_Trace = [T_dat.df_Trace df];
                T_dat.Freeze_Trace = [T_dat.Freeze_Trace false];
                T_dat.Melt_Trace = [T_dat.Melt_Trace false];
                while true
                    try
                        copyfile(Settings.CurrentTFile,Settings.PrevTFile)
                        save(Settings.CurrentTFile,'T_dat')
                        break
                    catch
                        pause(10)
                    end
                end
                feval = -1; % function evaluation
                fderiv = 0; % function derivative
                User_data = T_dat; % user data
                if Settings.Verbose
                    disp('System blow up. This potential may be unusable!')
                    disp('Aborting Melting Point calculation.')
                end
                return
            else
                delete(VolCheck_File)
                delete(VolCheck_Log_File)
            end
            
            PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(WorkDir, Settings.Salt, ...
                pyargs('SystemName',Settings.JobName,...
                'RefStructure',Settings.RefStructure,... % 'InitialRefFrac',Settings.Liquid_Fraction,...
                'RefChangeThreshold',Settings.MeltFreezeThreshold,...
                'FileType',Settings.CoordType,...
                'ML_TimeLength',0,...
                'ML_TimeStep',0,...
                'T_Ref',T_dat.T_ref));
            Froze = logical(PyOut{1});
            Melted = logical(PyOut{2});
            Froze_alt = logical(PyOut{6});
            IsNotComplete = ~Froze && ~Melted;
        catch
            IsNotComplete = true;
            Froze_alt = false;
        end
        ext_idx = str2double(cpt_number{1}) + 1; % increment step
        
        if IsNotComplete && ~Froze_alt && errcode ~= 0
            
            try % Clean up
                [~,~] = system([Settings.wsl 'find ' windows2unix(WorkDir) ' -iname "#*#" ^| xargs rm']);
            catch me
                disp(me.message)
            end
            
            if Settings.MDP.dt/2 >= Settings.MinTimeStep
                if Settings.Verbose
                    disp('Simulation failed. Restarting with reduced time step.')
                end
                
                % Look for and delete all checkpoint files
                cpt_check = dir(fullfile(WorkDir,[Settings.JobName '_*.cpt']));
                for idx = 1:length(cpt_check)
                    CheckPoint_File = fullfile(WorkDir,cpt_check(idx).name);
                    delete(CheckPoint_File)
                end
                
                Settings.MDP.dt = Settings.MDP.dt/2;
                Settings.Output_Coords = Settings.Output_Coords*2;
                [feval,fderiv,User_data] = Melting_Point_Check(T,Settings);
                return
            else
                f = -1;
                df = 0;
                T_dat.Alt_Structure = true;
                T_dat.T = T;
                T_dat.T_Trace = [T_dat.T_Trace T];
                T_dat.f_Trace = [T_dat.f_Trace f];
                T_dat.df_Trace = [T_dat.df_Trace df];
                T_dat.Freeze_Trace = [T_dat.Freeze_Trace false];
                T_dat.Melt_Trace = [T_dat.Melt_Trace false];
                while true
                    try
                        copyfile(Settings.CurrentTFile,Settings.PrevTFile)
                        save(Settings.CurrentTFile,'T_dat')
                        break
                    catch
                        pause(10)
                    end
                end
                feval = -1; % function evaluation
                fderiv = 0; % function derivative
                User_data = T_dat; % user data
                if Settings.Verbose
                    disp('Possible system blow up. This potential may be unusable!')
                    disp('Aborting Melting Point calculation.')
                end
                return
            end
        end
        
    else % No checkpoint file found
        
        % Number of MD steps to perform before ending each segment
        MD_nsteps = CheckTime/Settings.MDP.dt;
        
        % Maximum number of time segments
        max_steps = num2str(Settings.MaxCheckTime/CheckTime,'%03.f');
        
        % Copy over the structure, mdp, and top files
        while true
            try
                copyfile(Strucure_Ref_File,Strucure_In_File)
                copyfile(Sol_Ref_File,Sol_In_File)
                break
            catch
                pause(10)
            end
        end
        Topology_File = fullfile(OuterDir,[Settings.JobName '.top']);
        MDP_in_File = fullfile(WorkDir,[Settings.JobName '.mdp']);
        while true
            try
                copyfile(fullfile(OuterDir,[Settings.JobName '.mdp']),MDP_in_File)
                break
            catch
                pause(10)
            end
        end

        Traj_Conf_File = fullfile(WorkDir,[Settings.JobName '_001.tpr']);
        MDP_out_File = fullfile(WorkDir,[Settings.JobName '_out.mdp']);
        GromppLog_File = fullfile(WorkDir,[Settings.JobName '_Grompplog.log']);
        Table_File = fullfile(OuterDir,[Settings.JobName '_Table.xvg']);

        % Modify the temperature/number of steps of the mdp file and resave
        MDP_txt = fileread(MDP_in_File);
        MDP_txt = regexprep(MDP_txt,'(ref-t *= *).+?;',['$1' pad(num2str(T),18) ';']);
        MDP_txt = regexprep(MDP_txt,'(gen-temp *= *).+?;',['$1' pad(num2str(T),18) ';']);
        MDP_txt = regexprep(MDP_txt,'(nsteps *= *).+?;',['$1' pad(num2str(MD_nsteps),18) ';']);
        MDP_txt = regexprep(MDP_txt,'(nstxout += *)(.+?)( *);',['$1' num2str(round(Settings.Output_Coords)) '$3;']);
        MDP_txt = regexprep(MDP_txt,'(dt += *)(.+?)( *);',['$1' num2str(Settings.MDP.dt) '$3;']);
        
        fidMDP = fopen(MDP_in_File,'wt');
        fwrite(fidMDP,regexprep(MDP_txt,'\r',''));
        fclose(fidMDP);
        
        % Only applies to polarizable models
        ndx_filename = fullfile(WorkDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '.ndx']);
        ndx_add = add_polarization_shells(Settings,Strucure_In_File,...
            'ndx_filename',ndx_filename,'add_shells',false);
        
        % Run gmx grompp
        GROMPP_command = [Settings.gmx_loc Settings.grompp ' -c ' windows2unix(Strucure_In_File) ...
            ' -f ' windows2unix(MDP_in_File) ' -p ' windows2unix(Topology_File) ...
            ' -o ' windows2unix(Traj_Conf_File) ' -po ' windows2unix(MDP_out_File) ...
            ndx_add ' -maxwarn ' num2str(Settings.MaxWarn) Settings.passlog windows2unix(GromppLog_File)];
        [errcode,~] = system(GROMPP_command);

        % Catch error in grompp
        if errcode ~= 0
            error(['Error running GROMPP. Problem command: ' newline GROMPP_command]);
        end

        % Run the first segment of the simulation
        Log_File = fullfile(WorkDir,[Settings.JobName '.log']);
        Energy_file = fullfile(WorkDir,[Settings.JobName '.edr']);
        Trajectory_File = fullfile(WorkDir,[Settings.JobName '.trr']);
        Structure_Out_File = fullfile(WorkDir,[Settings.JobName '_OutConf_001.' Settings.CoordType]);
        CheckPoint_File = fullfile(WorkDir,[Settings.JobName '_001.cpt']);
        mdrun_command = [Settings.gmx Settings.mdrun ' -s ' windows2unix(Traj_Conf_File) ...
            ' -o ' windows2unix(Trajectory_File) ' -g ' windows2unix(Log_File) ...
            ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Structure_Out_File) ...
            ' -cpo ' windows2unix(CheckPoint_File) Settings.mdrun_opts];
        
        if isfile(Table_File) % If tabulated potential required
            mdrun_command = [mdrun_command ' -table ' windows2unix(Table_File)];
        end
        
        MDtimer = tic;
        if Settings.Verbose
            disp(['Beginning MD simulation segment 001/' max_steps ...
                ' (' num2str(CheckTime,'%.0f') '/' num2str(Settings.MaxCheckTime,'%.0f') ' ps).'])
        end
        
        if Settings.slurm
            mdrun_command = regexprep(mdrun_command,' -maxh ([0-9]|\.)+','','once');
            HoursRemain = JobTimeRemaining;
            mdrun_command = [mdrun_command ' -maxh ' num2str(HoursRemain)];
        else
            mdrun_command = regexprep(mdrun_command,' -maxh ([0-9]|\.)+','','once');
        end

        [errcode,outp] = system(mdrun_command);
        if errcode ~= 0
            try % Clean up
                [~,~] = system([Settings.wsl 'find ' windows2unix(WorkDir) ' -iname "#*#" ^| xargs rm']);
            catch me
                disp(me.message)
            end
            
            if Settings.MDP.dt/2 >= Settings.MinTimeStep
                if Settings.Verbose
                    disp('Simulation failed. Restarting with reduced time step.')
                end
                if isfile(CheckPoint_File)
                    delete(CheckPoint_File)
                end
                Settings.MDP.dt = Settings.MDP.dt/2;
                Settings.Output_Coords = Settings.Output_Coords*2;
                [feval,fderiv,User_data] = Melting_Point_Check(T,Settings);
                return
            else
                if Settings.Verbose
                    disp(['MD simulation segment 001/' max_steps ...
                        ' (' num2str(CheckTime,'%.0f') '/' num2str(Settings.MaxCheckTime,'%.0f') ...
                        ' ps) failed to complete. Time elapsed: ' datestr(seconds(toc(MDtimer)),'HH:MM:SS')])
                end

                % Check the final frame of the simulation chunk
                if Settings.Verbose
                    disp('Checking melting/freezing status...')
                end
                try
                    PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(WorkDir, Settings.Salt, ...
                        pyargs('SystemName',Settings.JobName,...
                        'RefStructure',Settings.RefStructure,... % 'InitialRefFrac',Settings.Liquid_Fraction,...
                        'RefChangeThreshold',Settings.MeltFreezeThreshold,...
                        'FileType',Settings.CoordType,...
                        'ML_TimeLength',0,...
                        'ML_TimeStep',0,...
                        'T_Ref',T_dat.T_ref));
                    Froze = logical(PyOut{1});
                    Melted = logical(PyOut{2});
                    Froze_alt = logical(PyOut{6});
                    IsNotComplete = ~Froze && ~Melted;
                catch
                    IsNotComplete = true;
                    Froze_alt = false;
                end
                
                if IsNotComplete && ~Froze_alt
                    f = -1;
                    df = 0;
                    T_dat.Alt_Structure = true;
                    T_dat.T = T;
                    T_dat.T_Trace = [T_dat.T_Trace T];
                    T_dat.f_Trace = [T_dat.f_Trace f];
                    T_dat.df_Trace = [T_dat.df_Trace df];
                    T_dat.Freeze_Trace = [T_dat.Freeze_Trace false];
                    T_dat.Melt_Trace = [T_dat.Melt_Trace false];
                    while true
                        try
                            copyfile(Settings.CurrentTFile,Settings.PrevTFile)
                            save(Settings.CurrentTFile,'T_dat')
                            break
                        catch
                            pause(10)
                        end
                    end
                    feval = -1; % function evaluation
                    fderiv = 0; % function derivative
                    User_data = T_dat; % user data
                    if Settings.Verbose
                        disp(outp);
                        disp('Possible system blow up. This potential may be unusable!')
                        disp('Aborting Melting Point calculation.')
                    end
                    return
                end
            end
        elseif ~isfile(Structure_Out_File)
            if Settings.Verbose
                disp('Calculation time is almost up, ending MATLAB session now.')
            end
            exit
        else
            if Settings.Verbose
                disp(['MD simulation segment 001/' max_steps ...
                    ' (' num2str(CheckTime,'%.0f') '/' num2str(Settings.MaxCheckTime,'%.0f') ...
                    ' ps) complete. Time elapsed: ' datestr(seconds(toc(MDtimer)),'HH:MM:SS')])
            end
            
            % Check the final frame of the simulation chunk
            if Settings.Verbose
                disp('Checking melting/freezing status...')
            end
            CheckStructureTimer = tic;
            
            % Check volume has not exploded
            VolCheck_File = fullfile(WorkDir,[Settings.JobName '_001_VolCheck.xvg']);
            VolCheck_Log_File = fullfile(WorkDir,[Settings.JobName '_001_VolCheck.log']);
            
            % Check energy options
            gmx_command = [Settings.wsl 'echo "0" ' Settings.pipe ...
                ' ' strrep(Settings.gmx_loc,Settings.wsl,'') Settings.g_energy ...
                ' -f ' windows2unix(Energy_file) ' -s ' windows2unix(Traj_Conf_File)];
            [~,outpt] = system(gmx_command);
            en_opts = regexp(outpt,'-+\n.+','match','once');
            En_set = '';
            En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Volume','tokens','once'))];
            En_set = [En_set ' 0'];
            En_set = regexprep(En_set,' +',' ');
            
            % Grab data from results
            VolCheck_command = [Settings.wsl 'echo ' En_set ' ' Settings.pipe ...
                ' ' strrep(Settings.gmx_loc,Settings.wsl,'') Settings.g_energy ...
                ' -f ' windows2unix(Energy_file) ' -o ' windows2unix(VolCheck_File) ...
                ' -s ' windows2unix(Traj_Conf_File) ' ' Settings.passlog windows2unix(VolCheck_Log_File)];
            [err,~] = system(VolCheck_command);
            if err ~= 0
                error('Failed to collect energy data.')
            end
            
            VolCheck_Data = import_xvg(VolCheck_File);
            dev_from_init = VolCheck_Data(:,2)./VolCheck_Data(1,2);
            max_dev = max(abs(dev_from_init));
            
            if max_dev >= 1.5
                if Settings.Verbose
                    disp(['System volume exploded (' num2str(max_dev) ...
                        ' deviation in volume detected). Aborting calculation.'])
                end
                
                f = -1;
                df = 0;
                T_dat.Alt_Structure = true;
                T_dat.T = T;
                T_dat.T_Trace = [T_dat.T_Trace T];
                T_dat.f_Trace = [T_dat.f_Trace f];
                T_dat.df_Trace = [T_dat.df_Trace df];
                T_dat.Freeze_Trace = [T_dat.Freeze_Trace false];
                T_dat.Melt_Trace = [T_dat.Melt_Trace false];
                while true
                    try
                        copyfile(Settings.CurrentTFile,Settings.PrevTFile)
                        save(Settings.CurrentTFile,'T_dat')
                        break
                    catch
                        pause(10)
                    end
                end
                feval = -1; % function evaluation
                fderiv = 0; % function derivative
                User_data = T_dat; % user data
                if Settings.Verbose
                    disp('System blow up. This potential may be unusable!')
                    disp('Aborting Melting Point calculation.')
                end
                return
            else
                delete(VolCheck_File)
                delete(VolCheck_Log_File)
            end
            
            PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(WorkDir, Settings.Salt, ...
                pyargs('SystemName',Settings.JobName,...
                'RefStructure',Settings.RefStructure,... % 'InitialRefFrac',Settings.Liquid_Fraction,...
                'RefChangeThreshold',Settings.MeltFreezeThreshold,...
                'FileType',Settings.CoordType,...
                'ML_TimeLength',0,...
                'ML_TimeStep',0,...
                'T_Ref',T_dat.T_ref,...
                'CheckFullTrajectory',true,...
                'SaveTrajectory',true,...
                'SavePredictionsImage',true,...
                'TimePerFrame',1));
            Froze = logical(PyOut{1});
            Melted = logical(PyOut{2});
            Froze_alt = logical(PyOut{6});

            IsNotComplete = ~Froze && ~Melted;
            ext_idx = 2;
        end
    end
    
    while IsNotComplete && ~Froze_alt && (CheckTime*ext_idx <= Settings.MaxCheckTime) 
        if Settings.Verbose
            disp(['Simulation inconclusive, extending simulation to step ' num2str(ext_idx,'%03.f') '/' max_steps ...
                ' (' num2str(CheckTime*ext_idx,'%.0f') '/' num2str(Settings.MaxCheckTime,'%.0f') ' ps).' ...
                ' Time elapsed: ' datestr(seconds(toc(CheckStructureTimer)),'HH:MM:SS')])
        end
        
        MDtimer = tic;
        % If simulation did not find an answer after 1 time chunk, keep trying. Update the input files.
        Traj_Conf_File_idx =  fullfile(WorkDir,[Settings.JobName '_' num2str(ext_idx,'%03.f') '.tpr']);
        Checkpoint_File_idx = fullfile(WorkDir,[Settings.JobName '_' num2str(ext_idx,'%03.f') '.cpt']);
        Structure_Out_File_idx = fullfile(WorkDir,[Settings.JobName '_OutConf_' num2str(ext_idx,'%03.f') '.' Settings.CoordType]);

        % Extend simulation by converting tpr file with gmx convert-tpr
        convert_tpr_command = [Settings.gmx_loc Settings.convert_tpr ' -s ' windows2unix(Traj_Conf_File) ...
            ' -o ' windows2unix(Traj_Conf_File_idx) ' -extend ' num2str(CheckTime)];

        [errcode,~] = system(convert_tpr_command);
        % Catch error in grompp
        if errcode ~= 0
            error(['Error running convert-tpr. Problem command: ' newline convert_tpr_command]);
        end

        % Set up next gmx mdrun command
        mdrun_command = [Settings.gmx Settings.mdrun ' -s ' windows2unix(Traj_Conf_File_idx) ...
            ' -o ' windows2unix(Trajectory_File) ' -g ' windows2unix(Log_File) ...
            ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Structure_Out_File_idx) ...
            ' -cpi ' windows2unix(CheckPoint_File) ' -cpo ' windows2unix(Checkpoint_File_idx) ...
            Settings.mdrun_opts ' -append'];
        
        if isfile(Table_File) % If tabulated potential required
            mdrun_command = [mdrun_command ' -table ' windows2unix(Table_File)]; %#ok<AGROW>
        end
        
        if Settings.slurm
            mdrun_command = regexprep(mdrun_command,' -maxh ([0-9]|\.)+','','once');
            HoursRemain = JobTimeRemaining;
            mdrun_command = [mdrun_command ' -maxh ' num2str(HoursRemain)]; %#ok<AGROW>
        else
            mdrun_command = regexprep(mdrun_command,' -maxh ([0-9]|\.)+','','once');
        end
        
        [errcode,~] = system(mdrun_command);
        if errcode == 0
            if Settings.Verbose
                disp(['MD simulation segment ' num2str(ext_idx,'%03.f') '/' max_steps ...
                    ' (' num2str(CheckTime*ext_idx,'%.0f') '/' num2str(Settings.MaxCheckTime,'%.0f') ...
                    ' ps) complete. Time elapsed: ' datestr(seconds(toc(MDtimer)),'HH:MM:SS')])
            end
        elseif ~isfile(Structure_Out_File)
            if Settings.Verbose
                disp('Calculation time is almost up, ending MATLAB session now.')
            end
            exit
        else
            if Settings.Verbose
                disp(['MD simulation segment ' num2str(ext_idx,'%03.f') '/' max_steps ...
                    ' (' num2str(CheckTime*ext_idx,'%.0f') '/' num2str(Settings.MaxCheckTime,'%.0f') ...
                    ' ps) failed to complete. Time elapsed: ' datestr(seconds(toc(MDtimer)),'HH:MM:SS')])
            end
        end

        if Settings.Verbose
            disp('Checking melting/freezing status...')
        end
        CheckStructureTimer = tic;
        try            
            % Check volume has not exploded
            VolCheck_File = fullfile(WorkDir,[Settings.JobName '_' num2str(ext_idx,'%03.f') '_VolCheck.xvg']);
            VolCheck_Log_File = fullfile(WorkDir,[Settings.JobName '_' num2str(ext_idx,'%03.f') '_VolCheck.log']);
            
            % Check energy options
            gmx_command = [Settings.wsl 'echo "0" ' Settings.pipe ...
                ' ' strrep(Settings.gmx_loc,Settings.wsl,'') Settings.g_energy ...
                ' -f ' windows2unix(Energy_file) ' -s ' windows2unix(Traj_Conf_File_idx)];
            [~,outpt] = system(gmx_command);
            en_opts = regexp(outpt,'-+\n.+','match','once');
            En_set = '';
            En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Volume','tokens','once'))]; %#ok<AGROW>
            En_set = [En_set ' 0']; %#ok<AGROW>
            En_set = regexprep(En_set,' +',' ');
            
            % Grab data from results
            VolCheck_command = [Settings.wsl 'echo ' En_set ' ' Settings.pipe ...
                ' ' strrep(Settings.gmx_loc,Settings.wsl,'') Settings.g_energy ...
                ' -f ' windows2unix(Energy_file) ' -o ' windows2unix(VolCheck_File) ...
                ' -s ' windows2unix(Traj_Conf_File_idx) ' ' Settings.passlog windows2unix(VolCheck_Log_File)];
            [err,~] = system(VolCheck_command);
            if err ~= 0
                error('Failed to collect energy data.')
            end
            
            VolCheck_Data = import_xvg(VolCheck_File);
            dev_from_init = VolCheck_Data(:,2)./VolCheck_Data(1,2);
            max_dev = max(abs(dev_from_init));
            
            if max_dev >= 1.5
                if Settings.Verbose
                    disp(['System volume exploded (' num2str(max_dev) ...
                        ' deviation in volume detected). Aborting calculation.'])
                end
                
                f = -1;
                df = 0;
                T_dat.Alt_Structure = true;
                T_dat.T = T;
                T_dat.T_Trace = [T_dat.T_Trace T];
                T_dat.f_Trace = [T_dat.f_Trace f];
                T_dat.df_Trace = [T_dat.df_Trace df];
                T_dat.Freeze_Trace = [T_dat.Freeze_Trace false];
                T_dat.Melt_Trace = [T_dat.Melt_Trace false];
                while true
                    try
                        copyfile(Settings.CurrentTFile,Settings.PrevTFile)
                        save(Settings.CurrentTFile,'T_dat')
                        break
                    catch
                        pause(10)
                    end
                end
                feval = -1; % function evaluation
                fderiv = 0; % function derivative
                User_data = T_dat; % user data
                if Settings.Verbose
                    disp('System blow up. This potential may be unusable!')
                    disp('Aborting Melting Point calculation.')
                end
                return
            else
                delete(VolCheck_File)
                delete(VolCheck_Log_File)
            end
            
            % Check the final frame of the simulation chunk for phase change
            PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(WorkDir, Settings.Salt, ...
                pyargs('SystemName',Settings.JobName,...
                'RefStructure',Settings.RefStructure,... % 'InitialRefFrac',Settings.Liquid_Fraction,...
                'RefChangeThreshold',Settings.MeltFreezeThreshold,...
                'FileType',Settings.CoordType,...
                'ML_TimeLength',0,...
                'ML_TimeStep',0,...
                'T_Ref',T_dat.T_ref));
            Froze = logical(PyOut{1});
            Melted = logical(PyOut{2});
            Froze_alt = logical(PyOut{6});
            IsNotComplete = ~Froze && ~Melted;

        catch
            IsNotComplete = true;
            Froze_alt = false;
        end
        % Update and increment
        CheckPoint_File = Checkpoint_File_idx;
        Structure_Out_File = Structure_Out_File_idx;
        Traj_Conf_File = Traj_Conf_File_idx;
        ext_idx = ext_idx+1;
        
%         % Delete previous cpt and tpr files
%         if Settings.Delete_Backups
%             delete(CheckPoint_File)
%             delete(Traj_Conf_File)
%             delete(Structure_Out_File)
%         end
        
        if IsNotComplete && ~Froze_alt && errcode ~= 0
            try % Clean up
                [~,~] = system([Settings.wsl 'find ' windows2unix(WorkDir) ' -iname "#*#" ^| xargs rm']);
            catch me
                disp(me.message)
            end
            if Settings.MDP.dt/2 >= Settings.MinTimeStep
                if Settings.Verbose
                    disp('Simulation failed with inconclusive result. Restarting with reduced time step.')
                end

                % Look for and delete all checkpoint files
                cpt_check = dir(fullfile(WorkDir,[Settings.JobName '_*.cpt']));
                for idx = 1:length(cpt_check)
                    CheckPoint_File = fullfile(WorkDir,cpt_check(idx).name);
                    delete(CheckPoint_File)
                end

                Settings.MDP.dt = Settings.MDP.dt/2;
                Settings.Output_Coords = Settings.Output_Coords*2;
                [feval,fderiv,User_data] = Melting_Point_Check(T,Settings);
                return
            else
                f = -1;
                df = 0;
                T_dat.Alt_Structure = true;
                T_dat.T = T;
                T_dat.T_Trace = [T_dat.T_Trace T];
                T_dat.f_Trace = [T_dat.f_Trace f];
                T_dat.df_Trace = [T_dat.df_Trace df];
                T_dat.Freeze_Trace = [T_dat.Freeze_Trace false];
                T_dat.Melt_Trace = [T_dat.Melt_Trace false];
                while true
                    try
                        copyfile(Settings.CurrentTFile,Settings.PrevTFile)
                        save(Settings.CurrentTFile,'T_dat')
                        break
                    catch
                        pause(10)
                    end
                end
                feval = -1; % function evaluation
                fderiv = 0; % function derivative
                User_data = T_dat; % user data
                if Settings.Verbose
                    disp('Possible system blow up. This potential may be unusable!')
                    disp('Aborting Melting Point calculation.')
                end
                return
            end
        end
    end

    % Once complete, calculate the function and "gradient"
    if Settings.Verbose
        disp('Checking full trajectory for time-to-phase-change...')
    end
    CheckStructureTimer = tic;
    PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(WorkDir, Settings.Salt, ...
        pyargs('SystemName',Settings.JobName,...
        'RefStructure',Settings.RefStructure,...
        'CheckFullTrajectory',true,...
        'SaveTrajectory',Settings.SaveTrajectory,...
        'SavePredictionsImage',Settings.SavePredictionsImage,... % 'InitialRefFrac',Settings.Liquid_Fraction,...
        'RefChangeThreshold',Settings.MeltFreezeThreshold,...
        'SlopeThreshold',Settings.SlopeThreshold,...
        'FileType',Settings.CoordType,...
        'ML_TimeLength',0,...
        'ML_TimeStep',0,...
        'T_Ref',T_dat.T_ref,...
        'T',T,...
        'TimePerFrame',10,...
        'SaveTrajectoryAux',2));
    if Settings.Verbose
        disp(['Completed trajectory check. Time elapsed: ' ...
            datestr(seconds(toc(CheckStructureTimer)),'HH:MM:SS')])
    end
    Froze = logical(PyOut{1});
    Melted = logical(PyOut{2});
    time_to_phase_change = PyOut{3};
    Froze_alt = logical(PyOut{6});
    
    if Froze_alt
        if time_to_phase_change < 1
            time_to_phase_change = 1;
        end
        f = -1;
        df = 0;
        T_dat.Alt_Structure = true;
        if Settings.Verbose
            disp(['Detected system freezing into unexpected structure at ' num2str(time_to_phase_change) ' ps.'])
            disp('Aborting Melting Point calculation.')
        end
    elseif ~Froze && ~Melted
        f = -1;
        df = 0;
        
        if Settings.Verbose
            disp(repmat('*',1,40));
            disp(['System has not melted or froze after max time of ' num2str(Settings.MaxCheckTime) ' ps. Indeterminate point found.'])
            disp(['Indeterminate point found at T = ' num2str(T,'%.4f') ' K. Error bounds: Tm = ' ...
                num2str(T_dat.dT(1),'%.4f') ' - ' num2str(T_dat.dT(2),'%.4f') ' K.']);
            disp(repmat('*',1,40));
        end
    elseif Froze
        if time_to_phase_change < 1
            time_to_phase_change = 1;
        end
        f = 1/time_to_phase_change;
        df = -5000/(time_to_phase_change);
        T_dat.T_Freeze_Trace = [T_dat.T_Freeze_Trace T];
        
        % Update the current T_mp lower error bound
        if T > T_dat.dT(1) && ~Settings.IgnoreBounds
            T_dat.dT(1) = T;
            T_dat.df_bracket(1) = df;
        end
        if Settings.Verbose
            disp(['Detected system freezing at ' num2str(time_to_phase_change) ' ps.'])
        end
    elseif Melted
        if time_to_phase_change < 1
            time_to_phase_change = 1;
        end
        f = 1/time_to_phase_change;
        df = 5000/(time_to_phase_change);
        T_dat.T_Melt_Trace = [T_dat.T_Melt_Trace T];
        
        % Update the current T_mp upper error bound
        if T < T_dat.dT(2) && ~Settings.IgnoreBounds
            T_dat.dT(2) = T;
            T_dat.df_bracket(2) = df;
        end
        if Settings.Verbose
            disp(['Detected system melting at ' num2str(time_to_phase_change) ' ps.'])
        end
    end
    
    % Save data to file before returning
    T_dat.T = T;
    T_dat.T_Trace = [T_dat.T_Trace T];
    T_dat.f_Trace = [T_dat.f_Trace f];
    T_dat.df_Trace = [T_dat.df_Trace df];
    T_dat.Freeze_Trace = [T_dat.Freeze_Trace Froze];
    T_dat.Melt_Trace = [T_dat.Melt_Trace Melted];
    while true
        try
            copyfile(Settings.CurrentTFile,Settings.PrevTFile)
            save(Settings.CurrentTFile,'T_dat')
            break
        catch
            pause(10)
        end
    end
    
    % Delete all but final cpt, tpr, and OutConf files
    if Settings.Delete_Backups
        
        % Delete checkpoint files
        cpt_check = dir(fullfile(WorkDir,[Settings.JobName '_*.cpt']));
        cpt_numbers = regexp({cpt_check.name},'_([0-9]{3})(_step[0-9]+)*(_prev)*.cpt','once','tokens');
        max_cpt = max(str2double(cellfun(@(v)v(1),cpt_numbers)));
        max_cpt_str = num2str(max_cpt,'%03.f');
        for idx = 1:length(cpt_check)
            if ~strcmp(cpt_numbers{idx}{1},max_cpt_str)
                delete(fullfile(cpt_check(idx).folder,cpt_check(idx).name))
            end
        end
        
        % Delete tpr files
        tpr_check = dir(fullfile(WorkDir,[Settings.JobName '_*.tpr']));
        tpr_numbers = regexp({tpr_check.name},'_([0-9]{3}).tpr','once','tokens');
        max_tpr = max(str2double(cellfun(@(v)v(1),tpr_numbers)));
        max_tpr_str = num2str(max_tpr,'%03.f');
        for idx = 1:length(tpr_check)
            if ~strcmp(tpr_numbers{idx}{1},max_tpr_str)
                delete(fullfile(tpr_check(idx).folder,tpr_check(idx).name))
            end
        end
        
        % Delete OutConf files
        gro_check = dir(fullfile(WorkDir,[Settings.JobName '_OutConf*.' Settings.CoordType]));
        gro_numbers = regexp({gro_check.name},['_([0-9]{3}).' Settings.CoordType],'once','tokens');
        max_gro = max(str2double(cellfun(@(v)v(1),gro_numbers)));
        max_gro_str = num2str(max_gro,'%03.f');
        for idx = 1:length(gro_check)
            if ~strcmp(gro_numbers{idx}{1},max_gro_str)
                delete(fullfile(gro_check(idx).folder,gro_check(idx).name))
            end
        end
        
        % Delete any backup files
        system([Settings.wsl 'find ' windows2unix(WorkDir) ' -iname "#*#" -delete']);
    end
    
    feval = f; % function evaluation
    fderiv = df; % function derivative
    User_data = T_dat; % user data 
end
