function Output = Setup_LiX_Simulation(Settings)

%% Sanity checks
if ~strcmpi(Settings.Annealing,'no')
    if Settings.Annealing_Times(end)/1000 ~= Settings.MDP.Trajectory_Time
        Settings.MDP.Trajectory_Time = Settings.Annealing_Times(end)/1000;
        if Settings.Verbose
            disp(['Warning: total trajectory time differant than final simulated annealing time point. Changing trajectory to ' num2str(Settings.MDP.Trajectory_Time) ' ns.']);
        end
    end
    
    if Settings.Annealing_Temps(1) ~= Settings.MDP.Initial_T || Settings.Annealing_Temps(1) ~= Settings.Target_T
        Settings.MDP.Initial_T = Settings.Annealing_Temps(1);
        Settings.Target_T = Settings.Annealing_Temps(1);
        if Settings.Verbose
            disp(['Warning: Initial temperature reset to match First Annealing Temp (T = ' num2str(Settings.MDP.Initial_T) ')']);
        end
    end
end

% Initialize Output
Output.StructureChange = false;
Output.SolidMelted = false;
Output.LiquidFroze = false;
Output.LiquidAmorphous = false;
Output.Aborted = false;

% Turns off lattice parameter expansion when false
if ~Settings.Expand_LP
    Settings.Expand_a_SC = 1;
    Settings.Expand_b_SC = 1;
    Settings.Expand_c_SC = 1;
    Settings.Expand_a = 1;
    Settings.Expand_b = 1;
    Settings.Expand_c = 1;
end

%% Preliminary calculation parameters

% Initialize cluster size to 0, it will be updated to a finite number as needed
Settings.Cluster_N = 0;

% Load Default Geometry info for Salt/Structure
if ~isfield(Settings,'Geometry') || ~strcmp(Settings.Geometry.Structure,Settings.Structure)
    Settings.Geometry = Default_Crystal(Settings,'Center_Coordinates',true);
end

% Get Metal and Halide info
[Settings.Metal,Settings.Halide] = Separate_Metal_Halide(Settings.Salt);
Metal_Info = elements('Sym',Settings.Metal);
Halide_Info = elements('Sym',Settings.Halide);

if ~isfield(Settings,'WorkDir')
    [Settings.WorkDir,Settings.JobName,Settings.Full_Model_Name] = GetMDWorkdir(Settings);
end
if Settings.Verbose
    disp(['Current Job: ' Settings.Salt ' ' Settings.JobName])
end

% Create the working directory if it does not exist
if ~exist(Settings.WorkDir,'dir')
    mkdir(Settings.WorkDir)
end

% Load Model parameters
% if ~isempty(Settings.Model)
%     [Settings,~] = Load_Model_Params(Settings);
% end

% Load Batch script (if applicable) settings and gromacs stuff
[Batch_Template,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts,Settings.postprocess] = MD_Batch_Template(Settings.JobSettings);

% Calculate the largest cutoff distance
Longest_Cutoff = max([Settings.MDP.RList_Cutoff Settings.MDP.RCoulomb_Cutoff Settings.MDP.RVDW_Cutoff]); % nm

% Calculate number of steps
MD_nsteps = Settings.MDP.Trajectory_Time*1000/Settings.MDP.dt; % Number of steps to perform before ending (should be 0 for single point energy calculations)

% Check if parameters can be input without a table;
Settings.Table_Req = IsGmxTableRequired(Settings);

% If required, step size of tabulated potentials in nm
if Settings.JobSettings.SinglePrecision
    Settings.Table_StepSize = 0.002;
else
    Settings.Table_StepSize = 0.0005;
end

% Grompp log filename
Settings.GrompLog_File = fullfile(Settings.WorkDir,[Settings.JobName '_Grompplog.log']);

% Trajectory parameter filename (gmx grompp output)
Settings.Traj_Conf_File = fullfile(Settings.WorkDir,[Settings.JobName '.tpr']);

% Load MDP template
MDP_Template = fileread(fullfile(Settings.home,'templates','Gromacs_Templates',...
'MDP_MD.template'));

%% Topology parameters

% Topology filename and directory for output
Settings.Topology_File = fullfile(Settings.WorkDir,[Settings.JobName '.top']);

% Load Topology template
Settings.Topology_Text = fileread(fullfile(Settings.home,'templates','Gromacs_Templates',...
'Topology.template'));

% Add in global parameters to Topology template
Settings.Topology_Text = strrep(Settings.Topology_Text,'##GENPAIRS##',Settings.Top_gen_pairs);
Settings.Topology_Text = strrep(Settings.Topology_Text,'##FUDGELJ##',num2str(Settings.Top_fudgeLJ));
Settings.Topology_Text = strrep(Settings.Topology_Text,'##FUDGEQQ##',num2str(Settings.Top_fudgeQQ));

% Insert element info into topology template
Settings.Topology_Text = strrep(Settings.Topology_Text,'##MET##',pad(Settings.Metal,2));
Settings.Topology_Text = strrep(Settings.Topology_Text,'##METZ##',pad(num2str(Metal_Info.atomic_number),3));
Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMASS##',pad(num2str(Metal_Info.atomic_mass),7));
Settings.Topology_Text = strrep(Settings.Topology_Text,'##MCHRG##',pad(num2str(Settings.S.Q),2));

Settings.Topology_Text = strrep(Settings.Topology_Text,'##HAL##',pad(Settings.Halide,2));
Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALZ##',pad(num2str(Halide_Info.atomic_number),3));
Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALMASS##',pad(num2str(Halide_Info.atomic_mass),7));
Settings.Topology_Text = strrep(Settings.Topology_Text,'##XCHRG##',pad(num2str(-Settings.S.Q),2));

if Settings.Table_Req

    % Define the function type as 1 (required for custom functions)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##NBFUNC##','1');

    % Define the combination rules (Lorenz-berthelot)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##COMBR##','1');

    % Define all the parameters as 1.0 (already included in potentials)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETC##',pad('1.0',10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALC##',pad('1.0',10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALC##',pad('1.0',10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETA##','1.0');
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALA##','1.0');
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALA##','1.0');

    % Generate tables of the potential
    if strcmp(Settings.Theory,'TF')
        [U_MX, U_MM, U_XX] = TF_Potential_Generator(Settings);
    elseif strcmp(Settings.Theory,'BH')
        [U_MX, U_MM, U_XX] = BH_Potential_Generator(Settings);
    elseif contains(Settings.Theory,'JC')
        switch Settings.Theory
            case 'JC'
                Settings.WaterModel = 'SPC/E';
            case 'JC3P'
                Settings.WaterModel = 'TIP3P';
            case 'JC4P'
                Settings.WaterModel = 'TIP4PEW';
            case 'JCSD'
                Settings.WaterModel = 'SD';
        end
        [U_MX, U_MM, U_XX] = JC_Potential_Generator(Settings);
    else
        error(['Warning: Unknown theory type: "' Settings.Theory '".'])
    end
    
    TableName = [Settings.JobName '_Table'];
    TableFile_MX = fullfile(Settings.WorkDir,[TableName '.xvg']);
    TableFile_MM = fullfile(Settings.WorkDir,[TableName '_' Settings.Metal '_' Settings.Metal '.xvg']);
    TableFile_XX = fullfile(Settings.WorkDir,[TableName '_' Settings.Halide '_' Settings.Halide '.xvg']);

    % Save tables into current directory
    fidMX = fopen(TableFile_MX,'wt');
    fwrite(fidMX,regexprep(U_MX,'\r',''));
    fclose(fidMX);

    fidMM = fopen(TableFile_MM,'wt');
    fwrite(fidMM,regexprep(U_MM,'\r',''));
    fclose(fidMM);

    fidXX = fopen(TableFile_XX,'wt');
    fwrite(fidXX,regexprep(U_XX,'\r',''));
    fclose(fidXX);
    
    % Modify the MDP file
    MDP_Template = strrep(MDP_Template,'##VDWTYPE##',pad('user',18));
    MDP_Template = strrep(MDP_Template,'##CUTOFF##',pad('group',18));
    MDP_Template = regexprep(MDP_Template,'ewald-rtol-lj.+?\n','');
    MDP_Template = regexprep(MDP_Template,'lj-pme-comb-rule.+?\n','');
    MDP_Template = regexprep(MDP_Template,'verlet-buffer-tolerance.+?\n','');
    MDP_Template = strrep(MDP_Template,'##RLIST##',pad(num2str(Settings.MDP.RList_Cutoff),18));
    MDP_Template = strrep(MDP_Template,'##RCOULOMB##',pad(num2str(Settings.MDP.RCoulomb_Cutoff),18));
    MDP_Template = strrep(MDP_Template,'##RVDW##',pad(num2str(Settings.MDP.RVDW_Cutoff),18));
    MDP_Template = strrep(MDP_Template,'##VDWMOD##',pad('none',18));

elseif contains(Settings.Theory,'JC')
    switch Settings.Theory
        case 'JC'
            Settings.WaterModel = 'SPC/E';
        case 'JC3P'
            Settings.WaterModel = 'TIP3P';
        case 'JC4P'
            Settings.WaterModel = 'TIP4PEW';
        case 'JCSD'
            Settings.WaterModel = 'SD';
    end

    TableFile_MX = '';

    % Definte the function type as 1 (LJ)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##NBFUNC##','1');

    % Define the combination rules (Lorenz-berthelot in sigma-epsilon form)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##COMBR##','2');

    % Get JC parameters
    [MX_JC_Param,MM_JC_Param,XX_JC_Param] = JC_Potential_Parameters(Settings);
    
    % Add parameters to topology text
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETC##',pad(num2str(MM_JC_Param.sigma,'%10.8e'),10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALC##',pad(num2str(XX_JC_Param.sigma,'%10.8e'),10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALC##',pad(num2str(MX_JC_Param.sigma,'%10.8e'),10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETA##',num2str(MM_JC_Param.epsilon,'%10.8e'));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALA##',num2str(XX_JC_Param.epsilon,'%10.8e'));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALA##',num2str(MX_JC_Param.epsilon,'%10.8e'));

    % Modify the MDP file
    MDP_Template = strrep(MDP_Template,'##VDWTYPE##',pad(Settings.MDP.VDWType,18));
    MDP_Template = strrep(MDP_Template,'##CUTOFF##',pad(Settings.MDP.CutOffScheme,18));
    MDP_Template = regexprep(MDP_Template,'energygrp-table.+?\n','');
    MDP_Template = regexprep(MDP_Template,'ewald-rtol-lj.+?\n','');
    MDP_Template = regexprep(MDP_Template,'lj-pme-comb-rule.+?\n','');
    MDP_Template = strrep(MDP_Template,'##RLIST##',pad(num2str(Settings.MDP.RList_Cutoff),18));
    MDP_Template = strrep(MDP_Template,'##RCOULOMB##',pad(num2str(Settings.MDP.RCoulomb_Cutoff),18));
    MDP_Template = strrep(MDP_Template,'##RVDW##',pad(num2str(Settings.MDP.RVDW_Cutoff),18));
    MDP_Template = strrep(MDP_Template,'##VDWMOD##',pad(Settings.MDP.vdw_modifier,18));

    % Add in Verlet Settings
    if strcmp(Settings.MDP.CutOffScheme,'Verlet')
        MDP_Template = strrep(MDP_Template,'##VerletBT##',pad(num2str(Settings.MDP.VerletBT),18));
    else
        MDP_Template = regexprep(MDP_Template,'verlet-buffer-tolerance.+?\n','');
    end
elseif contains(Settings.Theory,'BH')
    
    TableFile_MX = '';
    
    % Definte the function type as 2 (Buckingham)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##NBFUNC##','2');
    
    % Define the combination rule as 1 (Buckingham only has 1 comb rule)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##COMBR##','1');
    
    % Get BH parameters (cross terms are pre-computed using my combining rules)
    [U_MX,U_MM,U_XX] = BH_Potential_Parameters(Settings);
    
    % Add parameters to topology text
    % For BH potentials, parameter are B*exp(-alpha*r) + C/r^6
    % Parameter order is B alpha C
    Settings.Topology_Text = strrep(Settings.Topology_Text,'ptype  C          A','ptype   a              b           c6');
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETC##',[num2str(U_MM.B,'%10.8e') ' ' num2str(U_MM.alpha,'%10.8e')]);
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETA##',pad(num2str(U_MM.C,'%10.8e'),10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALC##',[num2str(U_XX.B,'%10.8e') ' ' num2str(U_XX.alpha,'%10.8e')]);
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALA##',pad(num2str(U_XX.C,'%10.8e'),10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALC##',[num2str(U_MX.B,'%10.8e') ' ' num2str(U_MX.alpha,'%10.8e')]);
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALA##',pad(num2str(U_MX.C,'%10.8e'),10));

    % Modify the MDP file
    MDP_Template = strrep(MDP_Template,'##VDWTYPE##',pad(Settings.MDP.VDWType,18));
    MDP_Template = strrep(MDP_Template,'##CUTOFF##',pad('group',18));
    MDP_Template = regexprep(MDP_Template,'energygrp-table.+?\n','');
    MDP_Template = regexprep(MDP_Template,'ewald-rtol-lj.+?\n','');
    MDP_Template = regexprep(MDP_Template,'lj-pme-comb-rule.+?\n','');
    MDP_Template = strrep(MDP_Template,'##RLIST##',pad(num2str(Settings.MDP.RList_Cutoff),18));
    MDP_Template = strrep(MDP_Template,'##RCOULOMB##',pad(num2str(Settings.MDP.RCoulomb_Cutoff),18));
    MDP_Template = strrep(MDP_Template,'##RVDW##',pad(num2str(Settings.MDP.RVDW_Cutoff),18));
    MDP_Template = strrep(MDP_Template,'##VDWMOD##',pad(Settings.MDP.vdw_modifier,18));
    MDP_Template = regexprep(MDP_Template,'verlet-buffer-tolerance.+?\n','');
    
else
    error(['Warning: Unknown theory type: "' Settings.Theory '".'])
end

%% Molecular dynamics parameters

% Add in parameters to MDP template
MDP_Template = strrep(MDP_Template,'##NSTEPS##',pad(num2str(MD_nsteps),18));
MDP_Template = strrep(MDP_Template,'##INTEGR##',pad(Settings.MDP.integrator,18));
MDP_Template = strrep(MDP_Template,'##TIMEST##',pad(num2str(Settings.MDP.dt),18));
MDP_Template = strrep(MDP_Template,'##CONTINUE##',pad(Settings.MDP.continuation,18));
MDP_Template = strrep(MDP_Template,'##REFTINIT##',pad(num2str(Settings.MDP.Initial_T),18));
MDP_Template = strrep(MDP_Template,'##COULOMB##',pad(Settings.MDP.CoulombType,18));
MDP_Template = strrep(MDP_Template,'##FOURIER##',pad(num2str(Settings.MDP.Fourier_Spacing),18));
MDP_Template = strrep(MDP_Template,'##PMEORDER##',pad(num2str(Settings.MDP.PME_Order),18));
MDP_Template = strrep(MDP_Template,'##EWALDTOL##',pad(num2str(Settings.MDP.Ewald_rtol),18));
MDP_Template = strrep(MDP_Template,'##LISTUPDATE##',pad(num2str(Settings.Update_NeighbourList),18));
MDP_Template = strrep(MDP_Template,'##POSOUT##',pad(num2str(Settings.Output_Coords),18));
MDP_Template = strrep(MDP_Template,'##POSOUTCOMP##',pad(num2str(Settings.Output_Coords_Compressed),18));
MDP_Template = strrep(MDP_Template,'##VELOUT##',pad(num2str(Settings.Output_Velocity),18));
MDP_Template = strrep(MDP_Template,'##FORCEOUT##',pad(num2str(Settings.Output_Forces),18));
MDP_Template = strrep(MDP_Template,'##ENOUT##',pad(num2str(Settings.Output_Energies),18));
MDP_Template = strrep(MDP_Template,'##CALCE##',pad(num2str(Settings.Calc_Energies),18));
MDP_Template = strrep(MDP_Template,'##MET##',pad(Settings.Metal,3));
MDP_Template = strrep(MDP_Template,'##HAL##',pad(Settings.Halide,3));

% Thermostat settings
if strcmpi(Settings.Thermostat,'no')
    % Remove thermostat inputs from template
    MDP_Template = strrep(MDP_Template,'##THERMOSTAT##',pad(Settings.Thermostat,18));
    MDP_Template = regexprep(MDP_Template,'tc-grps +.+?\n','');
    MDP_Template = regexprep(MDP_Template,'tau-t +.+?\n','');
    MDP_Template = regexprep(MDP_Template,'nsttcouple +.+?\n','');
    MDP_Template = regexprep(MDP_Template,'ref-t +.+?\n','');

else
    TC = '';
    TT = '';
    for x=1:Settings.MDP.Num_Groups
        TC = [num2str(Settings.Time_Constant_T) ' ' TC];
        TT = [num2str(Settings.Target_T) ' ' TT];
    end
    TC = pad(TC,18);
    TT = pad(TT,18);
    MDP_Template = strrep(MDP_Template,'##THERMOSTAT##',pad(Settings.Thermostat,18));
    MDP_Template = strrep(MDP_Template,'##TTIMECONST##',TC);
    MDP_Template = strrep(MDP_Template,'##NSTTCOUPLE##',pad(num2str(Settings.Nsttcouple),18));
    MDP_Template = strrep(MDP_Template,'##REFT##',TT);
end

% Barostat settings
if strcmpi(Settings.Barostat,'no')
    % Remove barostat inputs from template
    MDP_Template = strrep(MDP_Template,'##BAROSTAT##',pad(Settings.Barostat,18));
    MDP_Template = regexprep(MDP_Template,'pcoupltype +.+?\n','');
    MDP_Template = regexprep(MDP_Template,'tau-p +.+?\n','');
    MDP_Template = regexprep(MDP_Template,'compressibility +.+?\n','');
    MDP_Template = regexprep(MDP_Template,'ref-p +.+?\n','');
    MDP_Template = regexprep(MDP_Template,'nstpcouple +.+?\n','');
    
else
    % Get box compressibility based on the salt
    if Settings.Liquid_Interface && Settings.GenCluster % Automatically set liquid + cluster simulation to isotropic (solid) due to system shape.
        Settings.Isotropy = 'isotropic';
        Settings.Target_P = Settings.Target_P(1); % Bar
    end
    
    Compressibility = Get_Alkali_Halide_Compressibility(Settings.Salt,...
        'Isotropy',Settings.Isotropy,'Molten',Settings.UseMoltenCompressibility,...
        'ScaleFactor',Settings.ScaleCompressibility);
    
    MDP_Template = strrep(MDP_Template,'##BAROSTAT##',pad(Settings.Barostat,18));
    MDP_Template = strrep(MDP_Template,'##ISOTROPY##',pad(num2str(Settings.Isotropy),18));
    switch lower(Settings.Isotropy)
        case 'isotropic'
            Npres = 1;
        case 'semiisotropic'
            Npres = 2;
        case 'anisotropic'
            Npres = 6;
        case 'surface-tension'
            Npres = 2;
    end
    % Check to make sure length of input pressure vector is correct
    if length(Settings.Target_P) < Npres
        if Settings.Verbose
            disp('Warning: Not enough values of target pressure given, extending first value to undefined pressure.')
        end
        Settings.Target_P(end+1:Npres) = Settings.Target_P(1);
    elseif length(Settings.Target_P) > Npres
        if Settings.Verbose
            disp('Warning: Too many values of target pressure given, excluding extra values.')
        end
        Settings.Target_P(Npres+1:end) = [];
    end
    MDP_Template = strrep(MDP_Template,'##COMPRESS##',pad(num2str(Compressibility),18));
    MDP_Template = strrep(MDP_Template,'##REFP##',pad(num2str(Settings.Target_P),18));

    MDP_Template = strrep(MDP_Template,'##PTIMECONST##',pad(num2str(Settings.Time_Constant_P),18));
    MDP_Template = strrep(MDP_Template,'##NSTPCOUPLE##',pad(num2str(Settings.Nstpcouple),18));
    
end

if Settings.MDP.Disp_Correction && ~Settings.Table_Req
    MDP_Template = [MDP_Template newline newline ...
        '; Long-range dispersion correction' newline ...
        'DispCorr                 = EnerPres          ; apply long range dispersion corrections for Energy and pressure'];
elseif Settings.MDP.Disp_Correction && Settings.MDP.Disp_Correction_Tables
    disp('Warning: enabling long-range dispersion correction for tabulated potential!')
    MDP_Template = [MDP_Template newline newline ...
        '; Long-range dispersion correction' newline ...
        'DispCorr                 = EnerPres          ; apply long range dispersion corrections for Energy and pressure'];
elseif Settings.MDP.Disp_Correction
    if Settings.Verbose
        disp('Disabling long-range dispersion correction as this is not compatible with tables as implemented here.')
    end
end

% Simulated annealing
if strcmpi(Settings.Annealing,'no')
    % Remove annealing inputs from MDP
    MDP_Template = strrep(MDP_Template,'##ANNEALING##',pad(Settings.Annealing,18));
    MDP_Template = regexprep(MDP_Template,'annealing-npoints +.+?\n','');
    MDP_Template = regexprep(MDP_Template,'annealing-time +.+?\n','');
    MDP_Template = regexprep(MDP_Template,'annealing-temp +.+?\n','');
else                   
    if length(Settings.Annealing_Times) > 1
        AT = num2str(Settings.Annealing_Temps(1),'%10.0f');
        At = num2str(Settings.Annealing_Times(1),'%10.0f');
        for x = 2:length(Settings.Annealing_Times)
            AT = [AT '_' num2str(Settings.Annealing_Temps(x),'%10.0f')];
            At = [At '_' num2str(Settings.Annealing_Times(x),'%10.0f')];
        end
    else
        AT = num2str(Settings.Annealing_Temps,'%10.0f');
        At = num2str(Settings.Annealing_Times,'%10.0f');
    end

    % Prepare MDP inputs
    Anneal_types = pad(repmat([Settings.Annealing ' '],1,Settings.MDP.Num_Groups),18);
    Anneal_points = pad(repmat([num2str(length(Settings.Annealing_Times)) ' '],1,Settings.MDP.Num_Groups),18);
    Anneal_Temps = pad(repmat([strrep(AT,'_',' ') ' '],1,Settings.MDP.Num_Groups),18);
    Anneal_Times = pad(repmat([strrep(At,'_',' ') ' '],1,Settings.MDP.Num_Groups),19);

    % Update MDP template
    MDP_Template = strrep(MDP_Template,'##ANNEALING##',Anneal_types);
    MDP_Template = strrep(MDP_Template,'##ANNEALPNTS##',Anneal_points);
    MDP_Template = strrep(MDP_Template,'##ANNEALTIMES##',Anneal_Times);
    MDP_Template = strrep(MDP_Template,'##ANNEALTEMPS##',Anneal_Temps);
end

% If pbc off, add comm_mode = angular
if strcmp(Settings.pbc,'off')
    comm_txt = 'comm_mode                = Angular           ; Remove center of mass translational and rotational velocity around the center of mass';
    MDP_Template = regexprep(MDP_Template,'(pbc += )xyz(.+?)\n',['$1no $2' newline comm_txt]);
end

% Generate an MDP title
MDP_Title = [Settings.Structure ' ' Settings.Salt ' with ' Settings.Full_Model_Name ...
    ' model, ' Settings.Thermostat ...
    ' thermostat, and ' Settings.Barostat ' barostat.'];
MDP_Template = strrep(MDP_Template,'##TITLE##',MDP_Title);

% Create name for mdp output file
Settings.MDP_out_File = fullfile(Settings.WorkDir,[Settings.JobName '_out.mdp']);

% Save MDP file
Settings.MDP_in_File = fullfile(Settings.WorkDir,[Settings.JobName '.mdp']);
fidMDP = fopen(Settings.MDP_in_File,'wt');
fwrite(fidMDP,regexprep(MDP_Template,'\r',''));
fclose(fidMDP);

%% System geometry parameters

% Unit Cell Filename
Settings.UnitCellFile = fullfile(Settings.WorkDir,[Settings.JobName '_UnitCell.' Settings.CoordType]);

% Supercell Filename
Settings.SuperCellFile = fullfile(Settings.WorkDir,[Settings.JobName '.' Settings.CoordType]);

% Find minimum lattice parameter for this salt/structure/model (or use initial ones)
if Settings.Find_Min_Params && ~Settings.Skip_Minimization
    [Settings,Found_DataMatch] = FindMinLatticeParam(Settings,...
        'Find_Similar_Params',Settings.Find_Similar_Params,'Center_Coordinates',true);
elseif Settings.Skip_Minimization || strcmpi(Settings.Structure,'liquid') || strcmpi(Settings.Structure,'previous')
    Found_DataMatch = true;
else
    Found_DataMatch = false;
end

Premintxt = '';
% Structure Template filename
DoGeomEdit = true;
if strcmpi(Settings.Structure,'previous')
    [Settings.Coordinate_File,Settings.Geometry] = GetStructureFromTraj(Settings.Prev_geom_loc,Settings.Salt,...
        Settings.WorkDir,Settings.Prev_geom_time,Settings.Geometry,Settings.JobName,Settings.gmx_loc,Settings.CoordType);
    
elseif strcmpi(Settings.Structure,'liquid')
    
    MinDir = fullfile(Settings.WorkDir,'Minimization');
    if ~exist(MinDir,'dir')
        mkdir(MinDir);
    end
    JobFile = fullfile(MinDir,'TempJobInfo.mat');
    save(JobFile);
    Premintxt = ['matlab -batch "Create_Liquid_Box(''' MinDir ''')"' newline];
    DoGeomEdit = false;
else
    if ~Settings.Use_Conv_cell && contained_in_cell(Settings.Structure,{'Rocksalt' 'Sphalerite' 'FiveFive'})
        Settings.Coordinate_File = fullfile(Settings.home,'templates',[upper(Settings.CoordType) '_Templates'],...
            [Settings.Structure '_Primitive.' Settings.CoordType]);
    else
        Settings.Coordinate_File = fullfile(Settings.home,'templates',[upper(Settings.CoordType) '_Templates'],...
            [Settings.Structure '.' Settings.CoordType]);
    end
end

if DoGeomEdit

    % Load Coordinates text
    Coordinate_Text = fileread(Settings.Coordinate_File);

    % Add Metal and Halide symbols
    if strcmp(Settings.CoordType,'gro')
        Met = pad(Settings.Metal,2,'left');
        Hal = pad(Settings.Halide,2,'left');
    elseif strcmp(Settings.CoordType,'g96')
        Met = pad(Settings.Metal,2,'right');
        Hal = pad(Settings.Halide,2,'right');
    end
    Coordinate_Text = strrep(Coordinate_Text,'##MET##',Met);
    Coordinate_Text = strrep(Coordinate_Text,'##HAL##',Hal);

    % Initial guess size of supercell (should be 1 for liquids)
    [N_Supercell_a,N_Supercell_b,N_Supercell_c_tot] = SuperCellSize(Settings,Settings.Geometry);
    Sol_fraction = 1 - Settings.Liquid_Fraction;
    if Settings.Liquid_Interface && ~Settings.GenCluster
        N_Supercell_c = ceil(N_Supercell_c_tot*Sol_fraction);
    else
        N_Supercell_c = N_Supercell_c_tot;
    end

    % Expand the volume
    if Settings.Expand_LP && Settings.Thermal_Solid
        Exp_Coeff = LiX_Thermal_Expansion(Settings);
    else
        Exp_Coeff = 1;
    end

    Settings.Geometry.a = Settings.Geometry.a*Settings.Expand_a*Exp_Coeff;
    Settings.Geometry.b = Settings.Geometry.b*Settings.Expand_b*Exp_Coeff;
    Settings.Geometry.c = Settings.Geometry.c*Settings.Expand_c*Exp_Coeff;
    if isfield(Settings.Geometry,'xyz')
        Settings.Geometry = Update_xyz(Settings.Geometry);
    end
    
    % Manual box options
    if Settings.Manual_Box && Settings.GenCluster
        rho_sol = Settings.Geometry.NF/( Cell_Volume(Settings.Geometry.Transform(1,:).*Settings.Expand_a_SC*Settings.Geometry.a...
            ,Settings.Geometry.Transform(2,:).*Settings.Expand_b_SC*Settings.Geometry.b,...
            Settings.Geometry.Transform(3,:).*Settings.Expand_c_SC*Settings.Geometry.c)/(10^3) ); % molecules/nm^3
        Mols_Solid = round(( (4/3)*pi*Settings.Manual_Radius^3 )*rho_sol); % molecules
        Vol_Solid = Mols_Solid/rho_sol; % nm^3
        
        if Settings.Liquid_Interface
            rho_liq = Get_LiX_Liquid_Density(Settings); % molecules/nm^3
            F = Settings.Liquid_Fraction;
            
            N_total = 2*round(Mols_Solid/(1-F));
            N_Solid = 2*Mols_Solid;
            N_Liquid = N_total - N_Solid;
            
            Vol_Liquid = (N_Liquid/2)/rho_liq; % nm^3
            
            Settings.Vol_Total = Vol_Solid + Vol_Liquid; % nm^3
            
            L = Settings.Vol_Total^(1/3); % nm
            a_vec = L.*[1 0 0];
            b_vec = L.*[0 1 0];
            c_vec = L.*[0 0 1];
            
            L_solid = Settings.Manual_Radius*2.1; % nm
            a_vec_sol = L_solid.*[1 0 0];
            
            N_Supercell_a = ceil(norm(a_vec_sol*10)/Settings.Geometry.a);
            N_Supercell_b = N_Supercell_a;
            N_Supercell_c = N_Supercell_a;
            
            
            % Update the liquid fraction
            Settings.Liquid_Fraction = N_Liquid/N_total;
            Settings.Cluster_Frac = N_Solid/N_total;
        else
            a_vec = 2.1*Settings.Manual_Radius.*[1 0 0]; % nm
            b_vec = 2.1*Settings.Manual_Radius.*[0 1 0]; % nm
            c_vec = 2.1*Settings.Manual_Radius.*[0 0 1]; % nm
            
            Settings.Vol_Total = Cell_Volume(a_vec,b_vec,c_vec); % nm^3
            
            N_Supercell_a = round(norm(a_vec*10)/Settings.Geometry.a);
            N_Supercell_b = round(norm(b_vec*10)/Settings.Geometry.b);
            N_Supercell_c = round(norm(c_vec*10)/Settings.Geometry.c);
            
            N_total = Mols_Solid*2;
        end
        
        LatticeLength = min([Settings.Geometry.Skew_a*norm(a_vec) ...
                            Settings.Geometry.Skew_b*norm(b_vec) ...
                            Settings.Geometry.Skew_c*norm(c_vec)]);
        
        if LatticeLength/2 <= Longest_Cutoff*Settings.Cutoff_Buffer
            error('Selected Manual Box is has at least one lattice dimension too short for the given potential cut-off (including buffer).')
        end
        
    elseif Settings.Manual_Box && ~Settings.GenCluster
        N_Supercell_a = round((Settings.Manual_Box_a*10)/Settings.Geometry.a);
        N_Supercell_b = round((Settings.Manual_Box_a*10)/Settings.Geometry.b);
        N_Supercell_c = round((Settings.Manual_Box_c*10)/Settings.Geometry.c);
        
        if Settings.Liquid_Interface
            rho_liq = Get_LiX_Liquid_Density(Settings); % molecules/nm^3
            rho_sol = Settings.Geometry.NF/( Cell_Volume(Settings.Geometry.Transform(1,:).*Settings.Expand_a_SC*Settings.Geometry.a...
                ,Settings.Geometry.Transform(2,:).*Settings.Expand_b_SC*Settings.Geometry.b,...
                Settings.Geometry.Transform(3,:).*Settings.Expand_c_SC*Settings.Geometry.c)/(10^3) ); % molecules/nm^3
            F = Settings.Liquid_Fraction;
            
            c_sol = Settings.Manual_Box_c/( (F/(1-F))*(rho_sol/rho_liq) + 1 ); % nm
            N_Supercell_c = round((c_sol*10)/Settings.Geometry.c);
            c_sol_actual = N_Supercell_c*Settings.Geometry.c/10; % nm, Solid length to the nearest unit cell
            c_liq = Settings.Manual_Box_c - c_sol_actual;
            
            % Check to make sure supercell is large enough for the cutoffs by finding the shortest lattice parameter
            a_vec = Settings.Geometry.Transform(1,:).*Settings.Expand_a_SC*Settings.Geometry.a*N_Supercell_a/10; % nm
            b_vec = Settings.Geometry.Transform(2,:).*Settings.Expand_b_SC*Settings.Geometry.b*N_Supercell_b/10; % nm
            c_vec = Settings.Geometry.Transform(3,:).*Settings.Manual_Box_c; % nm
            
            N_Solid = N_Supercell_a*N_Supercell_b*N_Supercell_c*Settings.Geometry.N; % number of atoms
            N_Liquid = round(norm(cross(a_vec,b_vec))*c_liq*rho_liq)*2; % number of atoms
            N_total = N_Solid + N_Liquid;
            
            % Update the liquid fraction
            Settings.Liquid_Fraction = N_Liquid/N_total;
            
        else
            % Check to make sure supercell is large enough for the cutoffs by finding the shortest lattice parameter
            a_vec = Settings.Geometry.Transform(1,:).*Settings.Expand_a_SC*Settings.Geometry.a*N_Supercell_a/10; % nm
            b_vec = Settings.Geometry.Transform(2,:).*Settings.Expand_b_SC*Settings.Geometry.b*N_Supercell_b/10; % nm
            c_vec = Settings.Geometry.Transform(3,:).*Settings.Expand_c_SC*Settings.Geometry.c*N_Supercell_c/10; % nm
            N_total = (N_Supercell_a*N_Supercell_b*N_Supercell_c)*Settings.Geometry.N;
        end
        
        LatticeLength = min([Settings.Geometry.Skew_a*norm(a_vec) ...
                            Settings.Geometry.Skew_b*norm(b_vec) ...
                            Settings.Geometry.Skew_c*norm(c_vec)]);
        
        if LatticeLength/2 <= Longest_Cutoff*Settings.Cutoff_Buffer
            error('Selected Manual Box is has at least one lattice dimension too short for the given potential cut-off (including buffer).')
        end
    else
        
        % Check to make sure supercell is large enough for the cutoffs by finding the shortest lattice parameter
        a_vec = Settings.Geometry.Transform(1,:).*Settings.Expand_a_SC*Settings.Geometry.a*N_Supercell_a/10; % nm
        b_vec = Settings.Geometry.Transform(2,:).*Settings.Expand_b_SC*Settings.Geometry.b*N_Supercell_b/10; % nm
        c_vec = Settings.Geometry.Transform(3,:).*Settings.Expand_c_SC*Settings.Geometry.c*N_Supercell_c_tot/10; % nm
        
        if Settings.GenCluster % When generating clusters, reset the box vectors to cubic
            L = Cell_Volume(a_vec,b_vec,c_vec)^(1/3); % nm
            a_vec = L.*[1 0 0];
            b_vec = L.*[0 1 0];
            c_vec = L.*[0 0 1];
        end
        
        % Shortest lattice length in nm
        LatticeLength = min([Settings.Geometry.Skew_a*norm(a_vec) ...
                            Settings.Geometry.Skew_b*norm(b_vec) ...
                            Settings.Geometry.Skew_c*norm(c_vec)]);

        if LatticeLength/2 <= Longest_Cutoff*Settings.Cutoff_Buffer

            old_atnum = (N_Supercell_a*N_Supercell_b*N_Supercell_c_tot)*(Settings.Geometry.N);

            N_Supercell_a = ceil(2*(Longest_Cutoff*Settings.Cutoff_Buffer)/(Settings.Geometry.Skew_a*Settings.Geometry.a/10));
            N_Supercell_b = ceil(2*(Longest_Cutoff*Settings.Cutoff_Buffer)/(Settings.Geometry.Skew_b*Settings.Geometry.b/10));
            N_Supercell_c_tot = ceil(2*(Longest_Cutoff*Settings.Cutoff_Buffer)/(Settings.Geometry.Skew_c*Settings.Geometry.c/10));

            a_current = N_Supercell_a*Settings.Geometry.a/10;
            c_current = N_Supercell_c_tot*Settings.Geometry.c/10;

            if Settings.c_over_a > c_current/a_current
                % Expand c to maintain the user-selected c/a ratio
                N_Supercell_c_tot = ceil((Settings.Geometry.a/Settings.Geometry.c)*N_Supercell_a*Settings.c_over_a);
            elseif Settings.c_over_a < c_current/a_current
                % Expand a to maintain the user-selected c/a ratio
                N_Supercell_a = ceil((Settings.Geometry.c/Settings.Geometry.a)*N_Supercell_c_tot/Settings.c_over_a);
                N_Supercell_b = ceil((Settings.Geometry.c/Settings.Geometry.b)*N_Supercell_c_tot/Settings.c_over_a);
            end

            if Settings.Liquid_Interface && ~Settings.GenCluster
                N_Supercell_c = round(N_Supercell_c_tot*Sol_fraction);
            else
                N_Supercell_c = N_Supercell_c_tot;
            end
            if Settings.Verbose
                disp(['Warning: With ' num2str(old_atnum) ...
                    ' atoms, the cut-off length is longer than half the shortest box vector or longer than the smallest box diagonal element.'])
                disp(['Expanding the box to ' num2str((N_Supercell_a*N_Supercell_b*N_Supercell_c_tot)*(Settings.Geometry.N)) ' atoms.'])
            end
        end
        
        N_Solid = N_Supercell_a*N_Supercell_b*N_Supercell_c*Settings.Geometry.N; % number of atoms
        N_Liquid = 2*round((N_Solid/2)*(1/Sol_fraction - 1)); % number of atoms
        N_total = N_Solid + N_Liquid;
    end
    
    % Convert to Na x Nb x Nc supercell
    Na = num2str(N_Supercell_a);
    Nb = num2str(N_Supercell_b);
    Nc = num2str(N_Supercell_c);

    %% Geometry editing and adding coordinates
    if Found_DataMatch

        % Add number of unit cells to topology file
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##N##x##N##x##N##',[Na 'x' Nb 'x' Nc]);

        % Edit topology title
        if strcmpi(Settings.Structure,'liquid') || strcmpi(Settings.Structure,'previous')
            Settings.Topology_Text = strrep(Settings.Topology_Text,'##GEOM##',[Settings.Structure ' Cell']);
            Coordinate_Text = SaveGroFile(Settings.UnitCellFile,Settings.Geometry,false);
        elseif Settings.Liquid_Interface
            Settings.Topology_Text = strrep(Settings.Topology_Text,'##GEOM##','Liquid-Crystal Interface');
            Coordinate_Text = AddCartesianCoord(Coordinate_Text,Settings.Geometry,1,false,Settings.CoordType); %input coordinates
        else
            Settings.Topology_Text = strrep(Settings.Topology_Text,'##GEOM##','Crystal');
            Coordinate_Text = AddCartesianCoord(Coordinate_Text,Settings.Geometry,1,false,Settings.CoordType); %input coordinates
        end

        if sum([N_Supercell_a N_Supercell_b N_Supercell_c]) > 3
            % Save unit cell .gro file into current directory
            fid = fopen(Settings.UnitCellFile,'wt');
            fwrite(fid,regexprep(Coordinate_Text,'\r',''));
            fclose(fid);
            
            Supercell_command = [Settings.gmx_loc ' genconf -f ' windows2unix(Settings.UnitCellFile) ...
                 ' -o ' windows2unix(Settings.SuperCellFile) ' -nbox ' Na ' ' Nb ' ' Nc];
            
            if isfile(Settings.SuperCellFile)
                delete(Settings.SuperCellFile)
            end
            [errcode,output] = system(Supercell_command);
            
            if errcode ~= 0
                if Settings.Verbose
                    disp(output);
                end
                error(['Error creating supercell with genconf. Problem command: ' newline Supercell_command]);
            end
        else
            % Save super cell file into current directory
            fid = fopen(Settings.SuperCellFile,'wt');
            fwrite(fid,regexprep(Coordinate_Text,'\r',''));
            fclose(fid);
        end

        % build a cluster if requested
        if Settings.GenCluster
            Settings.Cluster_N = ceil(N_total*Settings.Cluster_Frac); % atoms
            if mod(Settings.Cluster_N,2) > 0
                Settings.Cluster_N = Settings.Cluster_N + 1; % Ensure number of atoms is even
            end
            
            Build_Cluster(Settings)
            if Settings.Liquid_Interface
                N_Liquid = Settings.Liquid_Fraction*Settings.Cluster_N/(1-Settings.Liquid_Fraction);
                N_Cell = Settings.Cluster_N + N_Liquid; %#ok<NASGU>
                N_total = Settings.Cluster_N + N_Liquid;
            else
                N_Cell = Settings.Cluster_N; %#ok<NASGU>
                N_total = Settings.Cluster_N;
            end
        end

        % Expand supercell by requested amount
        if Settings.Expand_a_SC ~= 1 || Settings.Expand_b_SC ~= 1 || Settings.Expand_c_SC ~= 1
            a_sc = num2str(Settings.Expand_a_SC*Settings.Geometry.a*N_Supercell_a/10,'%10.8e'); % supercell a length in nm
            b_sc = num2str(Settings.Expand_b_SC*Settings.Geometry.b*N_Supercell_b/10,'%10.8e'); % supercell b length in nm
            c_sc = num2str(Settings.Expand_c_SC*Settings.Geometry.c*N_Supercell_c/10,'%10.8e'); % supercell c length in nm

            bc = num2str(Settings.Geometry.alpha,'%10.4e');
            ac = num2str(Settings.Geometry.beta,'%10.4e');
            ab = num2str(Settings.Geometry.gamma,'%10.4e');
            
            if isfile(Settings.SuperCellFile)
                delete(Settings.SuperCellFile)
            end
            Expand_command = [Settings.gmx_loc ' editconf -f ' windows2unix(Settings.SuperCellFile) ...
                 ' -o ' windows2unix(Settings.SuperCellFile) ' -box ' a_sc ' ' b_sc ' ' c_sc ...
                 ' -angles ' bc ' ' ac ' ' ab];
            [errcode,output] = system(Expand_command);

            if errcode ~= 0
                if Settings.Verbose
                    disp(output);
                end
                error(['Error expanding supercell with editconf. Problem command: ' newline Expand_command]);
            end
        end

        if Settings.Liquid_Interface
            MinDir = fullfile(Settings.WorkDir,'Minimization');
            if ~exist(MinDir,'dir')
                mkdir(MinDir);
            end
            JobFile = fullfile(MinDir,'TempJobInfo.mat');
            save(JobFile);
            Premintxt = ['matlab -batch "MD_Preminimization(''' MinDir ''')"' newline];
        else
            % Save number of atoms into .mat file
            NumberFile = fullfile(Settings.WorkDir,[Settings.JobName '.mat']);
            N_Cell = Settings.Geometry.N;
            save(NumberFile,'N_total','N_Cell')

            % Generate topology file
            Atomlist = copy_atom_order(Settings.SuperCellFile);
            Settings.Topology_Text = strrep(Settings.Topology_Text,'##LATOMS##',Atomlist);
            fidTOP = fopen(Settings.Topology_File,'wt');
            fwrite(fidTOP,regexprep(Settings.Topology_Text,'\r',''));
            fclose(fidTOP);

            GROMPP_command = [Settings.gmx_loc ' grompp -c ' windows2unix(Settings.SuperCellFile) ...
                ' -f ' windows2unix(Settings.MDP_in_File) ' -p ' windows2unix(Settings.Topology_File) ...
                ' -o ' windows2unix(Settings.Traj_Conf_File) ' -po ' windows2unix(Settings.MDP_out_File) ...
                ' -maxwarn ' num2str(Settings.MaxWarn) Settings.passlog windows2unix(Settings.GrompLog_File)];
            [errcode,~] = system(GROMPP_command);
            
            % Catch error in grompp
            if errcode ~= 0
                error(['Error running GROMPP. Problem command: ' newline GROMPP_command]);
            end
        end
    else
        MinDir = fullfile(Settings.WorkDir,'Minimization');
        if ~exist(MinDir,'dir')
            mkdir(MinDir);
        end
        JobFile = fullfile(MinDir,'TempJobInfo.mat');
        save(JobFile);
        Premintxt = ['matlab -batch "MD_Preminimization(''' MinDir ''')"' newline];
    end
end

if Settings.PreEquilibration > 0
    MinDir = fullfile(Settings.WorkDir,'Minimization');
    if ~exist(MinDir,'dir')
        mkdir(MinDir);
    end
    JobFile = fullfile(MinDir,'TempJobInfo.mat');
    if ~isfile(JobFile)
        save(JobFile);
    end
    Premintxt = [Premintxt 'matlab -batch "MD_PreEquilibrate(''' MinDir ''')"' newline];
end

%% Check for previous checkpoint files
Files = dir(Settings.WorkDir);

% Find only files with .cpt extention
CptIndex = cellfun(@(a)~isempty(a),arrayfun(@(a) regexp(a,'\.cpt$'),{Files.name}));
Files_Of_Interest = Files(CptIndex);

sFOI = size(Files_Of_Interest,1);

if sFOI > 1
    % If multiple cpt files exist, choose most recent
    Times = zeros(1,sFOI);
    for y = 1:sFOI
        Times(y) = datenum(Files_Of_Interest(y).date,...
            'dd-mmm-yyyy HH:MM:SS');
    end
    [~,Idx] = max(Times);
    Files_Of_Interest = Files_Of_Interest(Idx);
end

if sFOI == 0
    CheckPoint_In = '';
else
    CheckPoint_In = [' -cpi ' fullfile(Settings.WorkDir,Files_Of_Interest.name)];
    disp(['Starting from previous checkpoint file: ' Files_Of_Interest.name]);
end


%% Prepare mdrun and cleanup commands
Log_File = fullfile(Settings.WorkDir,[Settings.JobName '.log']);
Energy_file = fullfile(Settings.WorkDir,[Settings.JobName '.edr']);
Trajectory_File = fullfile(Settings.WorkDir,[Settings.JobName '.trr']);
ConfOut_File = fullfile(Settings.WorkDir,[Settings.JobName '_OutConf.' Settings.CoordType]);
CheckPoint_File = fullfile(Settings.WorkDir,[Settings.JobName '.cpt']);
mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(Settings.Traj_Conf_File) ...
    ' -o ' windows2unix(Trajectory_File) ' -g ' windows2unix(Log_File) ...
    ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(ConfOut_File) ...
    ' -cpo ' windows2unix(CheckPoint_File) Settings.mdrun_opts];

if Settings.Table_Req % If tabulated potential required
    mdrun_command = [mdrun_command ' -table ' windows2unix(TableFile_MX)];
end

% Build the cleanup command
cleanup_command = '';    
if Settings.Delete_subms
    cleanup_command = [cleanup_command Settings.wsl 'find ' windows2unix(Settings.WorkDir) ...
        ' -iname "*-0*.subm" -delete' newline];
end
if Settings.Delete_outputs 
    cleanup_command = [cleanup_command Settings.wsl 'find ' windows2unix(Settings.WorkDir) ...
        ' -iname "*.std*" -delete' newline];
end
if Settings.Delete_MDPout
    cleanup_command = [cleanup_command Settings.wsl 'rm ' Settings.MDP_out_File newline]; %#ok<*UNRCH>
end
if Settings.Delete_MDlog
    cleanup_command = [cleanup_command Settings.wsl 'rm ' windows2unix(Log_File) newline];
end
if Settings.Delete_ConfOut
    cleanup_command = [cleanup_command Settings.wsl 'rm ' windows2unix(ConfOut_File) newline];
end
if Settings.Delete_TPR
    cleanup_command = [cleanup_command Settings.wsl 'rm ' windows2unix(Settings.Traj_Conf_File) newline];
end
if Settings.Delete_Backups
    cleanup_command = [cleanup_command Settings.wsl 'find ' windows2unix(Settings.WorkDir) ...
        ' -iname "#*#" -delete' newline]; %#ok<*AGROW>
end
if Settings.Delete_cpt
    cleanup_command = [cleanup_command Settings.wsl 'ls -tp | grep "cpt" | tail -n +2 | xargs rm; for i in *cpt; do mv "$i" "$(echo "$i" | sed -E ''s/-[0-9][0-9][0-9]\.cpt/.cpt/'')"; done' newline]; %#ok<*AGROW>
end

% Ensure the initial geometry file is available in gro format for post-processing
if ~strcmp(Settings.CoordType,'gro')
    SuperCellFileGro = fullfile(Settings.WorkDir,[Settings.JobName '.gro']);
    trjconv_command = [Settings.gmx_loc ' trjconv -s ' windows2unix(Settings.Traj_Conf_File) ...
        ' -f ' windows2unix(Settings.SuperCellFile) ' -o ' windows2unix(SuperCellFileGro)];

    trjconv_command = regexprep(trjconv_command,'gmx',['echo 0 ' Settings.pipe ' gmx'],'once');
    cleanup_command = [cleanup_command trjconv_command newline];
end

% Ensure cleanup command only runs after job is complete
cleanup_command = ['if [[ -f "' ConfOut_File '" ]]; then' newline ...
    regexprep(cleanup_command,{'\n(.)' '^(.)'},{'\n    $1' '    $1'}) ...
    'fi' newline];

%% Submitting / running job
if ~Settings.BatchMode % Running job locally
    
    if Settings.Verbose
        disp('Starting Job Locally.')
    end
    if ~Found_DataMatch || Settings.Liquid_Interface
        Output = MD_Preminimization(MinDir);
        if Output.Aborted
            return
        end
    elseif strcmpi(Settings.Structure,'liquid')
        Create_Liquid_Box(MinDir)
    end
    
    if Settings.PreEquilibration > 0
        Output = MD_PreEquilibrate(MinDir);
        if Output.Aborted
            return
        end
    end

    % Make sure the initial configuration file is also available in gro format
    if ~strcmp(Settings.CoordType,'gro')
        [errcode,outp] = system(trjconv_command,'-echo');

        % Catch error
        if errcode ~= 0
            if Settings.Verbose
                disp(outp);
            end
            error(['Error running TRJCONV. Problem command: ' newline trjconv_command]);
        end
    end

    if Settings.Submit_Jobs
        [err,outp] = system([mdrun_command CheckPoint_In]);
        mdlog_file = fullfile(Settings.WorkDir,[Settings.JobName '.mdlog']);
        fidBS = fopen(mdlog_file,'wt');
        fwrite(fidBS,regexprep(outp,'\r',''));
        fclose(fidBS);
        if err ~= 0 
            error(['Problem with mdrun command. See log file: ' mdlog_file]);
        end
        
        if Settings.RunPostProcessor
            Settings.PostProcessFlagFile = fullfile(Settings.WorkDir,'POSTPROCESS_COMPLETE');
            MD_Postprocessor(Settings);
        end
    end
    if ~isempty(cleanup_command)
        [~,~] = system(cleanup_command);
    end

elseif Settings.BatchMode && ~isempty(Settings.qsub) % Running job in batch mode if possible
    
    if Settings.RunPostProcessor
        Settings.PostProcessFlagFile = fullfile(Settings.WorkDir,'POSTPROCESS_COMPLETE');
        SettingsFile = fullfile(Settings.WorkDir,'PostSettings.mat');
        save(SettingsFile,'Settings');        
        postprocess_command = ['if [[ ! -f "' Settings.PostProcessFlagFile '" && -f "' ConfOut_File '" ]]; then' newline ...
            Settings.postprocess ...
            '    matlab -batch "MD_Postprocessor(''' SettingsFile ''')"' newline ...
            'fi'];
        
        cleanup_command = [postprocess_command newline ...
            cleanup_command];
    end    
    
    % Prepare batch script
    Batch_Text = Batch_Template; % Copy over a new batch template for each job replicate
    EXT1 = ['if [[ ! -f "' ConfOut_File '" ]]; then'];
    EXT2 = 'fi';
    Batch_Text = strrep(Batch_Text,'##EXT1##',EXT1);
    Batch_Text = strrep(Batch_Text,'##EXT2##',EXT2);
    if isempty(Premintxt)
        Batch_Text = strrep(Batch_Text,['##PREMIN##' newline],''); % add in pre-minimization for the first job only
    else
        Batch_Text = strrep(Batch_Text,['##PREMIN##' newline],[EXT1 newline '    ' Premintxt EXT2 newline]); % add in pre-minimization for the first job only
    end
    Batch_Text = strrep(Batch_Text,'##MDRUN##',[mdrun_command CheckPoint_In]);
    Batch_Text = strrep(Batch_Text,'##CLEANUP##',cleanup_command);
    Batch_Text = strrep(Batch_Text,'##TASKNAME##',[Settings.Salt '_' Settings.JobName]);
    Batch_Text = strrep(Batch_Text,'##ERROR##',windows2unix([Settings.WorkDir filesep Settings.JobName]));
    Batch_Text = strrep(Batch_Text,'##DIRECTORY##',windows2unix(Settings.WorkDir));

    % Shorten max time to allow for pre-minimization
    if ~Found_DataMatch || Settings.Liquid_Interface || Settings.PreEquilibration > 0
        Batch_Text = regexprep(Batch_Text,'-maxh ([0-9]|\.)',['-maxh ' num2str(Settings.JobSettings.Hours-0.17)]); % ~10 minutes
    elseif strcmpi(Settings.Structure,'liquid')
        Batch_Text = regexprep(Batch_Text,'-maxh ([0-9]|\.)',['-maxh ' num2str(Settings.JobSettings.Hours-0.017)]); % ~1 minute
    end

    % Open and save batch script
    batch_file = fullfile(Settings.WorkDir,[Settings.JobName '.subm']);
    fidBS = fopen(batch_file,'wt');
    fwrite(fidBS,regexprep(Batch_Text,{'\r' Settings.wsl},{'' ''}));
    fclose(fidBS);

    % Submit or run job
    if Settings.Verbose
        disp('Batch job input files produced.')
    end
    
    if Settings.Submit_Jobs
        pdir = pwd;
        cd(Settings.WorkDir);
        
        [~,output] = system([Settings.qsub ' ' batch_file]);
        if Settings.Verbose
            disp('Job (Link 1) submitted.')
        end

        PrevJobID = regexp(output,'[0-9]+','match','ONCE');
        cpt_output_prev = CheckPoint_File;

        % Make additional links
        for jdx = 2:Settings.JobSettings.N_Calc 

            Index = ['-' num2str(jdx,'%03.0f')];

            cpt_output_i = fullfile(Settings.WorkDir,[Settings.JobName Index '.cpt']);
            JobName_i = [Settings.JobName Index];

            % Edit mdrun
            mdrun_command_i = regexprep(mdrun_command,'-cpo [^\s]+ *',['-cpo ' cpt_output_i ' ']);
            mdrun_command_i = [mdrun_command_i ' -cpi ' cpt_output_prev];

            % Prevent jobs from re-starting if they finish, or fail without producing a checkpoint
            EXT1 = ['if [[ ! -f "' ConfOut_File '" && -f "' cpt_output_prev '" ]]; then'];
            EXT2 = 'fi';
            
            % Preminimization should be completed by the second job link
            Batch_Text_i = strrep(Batch_Template,['##PREMIN##' newline],'');
            Batch_Text_i = strrep(Batch_Text_i,'##ERROR##',fullfile(Settings.WorkDir,JobName_i));
            Batch_Text_i = strrep(Batch_Text_i,'##TASKNAME##',[Settings.Salt '_' JobName_i]);
            Batch_Text_i = strrep(Batch_Text_i,'##DIRECTORY##',Settings.WorkDir);
            Batch_Text_i = strrep(Batch_Text_i,'##EXT1##',EXT1);
            Batch_Text_i = strrep(Batch_Text_i,'##EXT2##',EXT2);
            Batch_Text_i = strrep(Batch_Text_i,'##MDRUN##',mdrun_command_i);
            Batch_Text_i = strrep(Batch_Text_i,'##CLEANUP##',cleanup_command);

            % Save batch script
            batch_file = fullfile(Settings.WorkDir,[JobName_i '.subm']);
            fidBS = fopen(batch_file,'wt');
            fwrite(fidBS,regexprep(Batch_Text_i,{'\r' Settings.wsl},{'' ''}));
            fclose(fidBS);

            if Settings.slurm
                [~,output] = system([Settings.qsub ' --depend=afterany:' PrevJobID ' ' batch_file]);
            else
                [~,output] = system([Settings.qsub ' -W depend=afterany:' PrevJobID ' ' batch_file]);
            end
            
            if Settings.Verbose
                disp(['Job (Link ' num2str(jdx) ') submitted.'])
            end

            % Update previous job stuff
            cpt_output_prev = cpt_output_i;
            PrevJobID = regexp(output,'[0-9]+','match','ONCE');
        end
        cd(pdir); % Return to previous location
    end
else
    if Settings.Verbose
        disp('Unable to run job in batch mode: no batch program available')
    end
end
