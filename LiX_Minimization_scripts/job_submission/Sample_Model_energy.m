% Sample_Model_energy
function Output = Sample_Model_energy(Settings,varargin)
%#ok<*UNRCH>
%#ok<*NASGU>
%#ok<*ASGLU>

% Parse optional inputs
p = inputParser;
p.FunctionName = 'Sample_Model_energy';

addOptional(p,'Delete_output_files',true,@(x)validateattributes(x,{'logical'},{'nonempty'}))
addOptional(p,'N_atoms',100,@(x)validateattributes(x,{'numeric'},{'nonempty'}))
addOptional(p,'Use_Conv_Cell',false,@(x)validateattributes(x,{'logical'},{'nonempty'}))
addOptional(p,'Extra_Properties',true,@(x)validateattributes(x,{'logical'},{'nonempty'}))

parse(p,varargin{:})


%% Main Settings
Use_Conv_Cell = p.Results.Use_Conv_Cell; % When true, minimizes in terms of the conventional unit cell (may be slower, especially for 5-5)
Extra_Properties = p.Results.Extra_Properties; % Gathers components of energy at the end in addition to the total energy.
Delete_output_files = p.Results.Delete_output_files; % When true: runs calculation in a temp folder and deletes the calculation folder when finished
N_atoms = p.Results.N_atoms; % Minimum number of atoms to include in super cell

%% ADDITIONAL SETTINGS
Settings.Longest_Cutoff = max([Settings.MinMDP.RList_Cutoff Settings.MinMDP.RCoulomb_Cutoff Settings.MinMDP.RVDW_Cutoff]);
Settings.Table_Length = Settings.Longest_Cutoff + 1.01; % nm. This should be at least equal rc+1 where rc is the largest cutoff

if Settings.SinglePrecision
    Settings.Table_StepSize = 0.002; % Step size of tabulated potentials in nm
else
    Settings.Table_StepSize = 0.0005; % 0.0005 Step size of tabulated potentials in nm
end

if ~isfield(Settings,'gmx')
    [~,Settings] = MD_Batch_Template(Settings);
end

if Settings.MinMDP.Parallel_Min % Since using matlab parallel, make sure gromacs runs in serial
    env.OMP_NUM_THREADS = getenv('OMP_NUM_THREADS');
    env.GMX_PME_NUM_THREADS = getenv('GMX_PME_NUM_THREADS');
    env.GMX_PME_NTHREADS = getenv('GMX_PME_NTHREADS');
    env.GMX_OPENMP_MAX_THREADS = getenv('GMX_OPENMP_MAX_THREADS');
    env.KMP_AFFINITY = getenv('KMP_AFFINITY');
    
    setenv('OMP_NUM_THREADS','1');
    setenv('GMX_PME_NUM_THREADS','1');
    setenv('GMX_PME_NTHREADS','1');
    setenv('GMX_OPENMP_MAX_THREADS','1');
    setenv('KMP_AFFINITY','disabled');
    Settings.mdrun_opts = ' -pin on -ntmpi 1 -ntomp 1';
    Settings.gmx = Settings.gmx_loc;
end


% Get system geometry
% Generate name for model with current scaling parameters
Model_Scaled = ModelName(Settings);

if Delete_output_files
    % Keep the model name short
    Model_Scaled = Settings.Theory;
end

% Load topology template location
Topology_Template = fullfile(Settings.home,'templates','Gromacs_Templates',...
'Topology.template');

% Load Topology template
Topology_Template = fileread(Topology_Template);

% Add in global parameters to Topology template
Topology_Template = strrep(Topology_Template,'##GENPAIRS##',Settings.Top.gen_pairs);
Topology_Template = strrep(Topology_Template,'##FUDGELJ##',num2str(Settings.Top.fudgeLJ));
Topology_Template = strrep(Topology_Template,'##FUDGEQQ##',num2str(Settings.Top.fudgeQQ));

% Load mdp template location
Settings.MinMDP.Template = fullfile(Settings.home,'templates','Gromacs_Templates',...
'MDP.template');

% Load MDP template
Settings.MinMDP.Template = fileread(Settings.MinMDP.Template);

% Add in global parameters to MDP template
Settings.MinMDP.Template = strrep(Settings.MinMDP.Template,'##NSTEPS##',pad(num2str(Settings.MinMDP.nsteps_point),18));
Settings.MinMDP.Template = strrep(Settings.MinMDP.Template,'##INTEGR##',pad(Settings.MinMDP.point_integrator,18));
Settings.MinMDP.Template = strrep(Settings.MinMDP.Template,'##TIMEST##',pad(num2str(Settings.MinMDP.dt),18));

% Boolean: check if parameters can be input into JC without a table;
Table_Req = IsGmxTableRequired(Settings);

% Get current label
Label = Settings.Geometry.Label;

% Get Metal and Halide info from Current Salt
[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);
if strcmp(Metal,'M')
    Metal_Info.atomic_mass = 1;
    Metal_Info.atomic_number = 1;
else
    Metal_Info = elements('Sym',Metal);
end
if strcmp(Halide,'X')
    Halide_Info.atomic_mass = 1;
    Halide_Info.atomic_number = 1;
else
    Halide_Info = elements('Sym',Halide);
end

% Copy Topology Template
Topology_Temp_Ions = Topology_Template;

% Insert element info into topology template
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##MET##',pad(Metal,2));
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##METZ##',pad(num2str(Metal_Info.atomic_number),3));
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##METMASS##',pad(num2str(Metal_Info.atomic_mass),7));
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##MCHRG##',pad(num2str(Settings.S.Q),2));

Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##HAL##',pad(Halide,2));
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##HALZ##',pad(num2str(Halide_Info.atomic_number),3));
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##HALMASS##',pad(num2str(Halide_Info.atomic_mass),7));
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##XCHRG##',pad(num2str(-Settings.S.Q),2));

% Copy MDP Template
Settings.MinMDP.Temp_Ions = Settings.MinMDP.Template;

% Insert salt components into MDP template
Settings.MinMDP.Temp_Ions = strrep(Settings.MinMDP.Temp_Ions,'##MET##',Metal);
Settings.MinMDP.Temp_Ions = strrep(Settings.MinMDP.Temp_Ions,'##HAL##',Halide);

% Copy Topology template
Topology_Temp_Struct = Topology_Temp_Ions;

% Get Structure Template filename
if Settings.Geometry.Conv
    conv = '';
else
    switch Settings.Structure
        case {'FiveFive' 'Rocksalt' 'Sphalerite'}
            conv = '_Primitive';
        otherwise
            
            conv = '';
    end
end
Coordinate_File = fullfile(Settings.home,'templates',...
    [upper(Settings.CoordType) '_Templates'],...
    [Settings.Structure conv '.' Settings.CoordType]);

% Load Template text from disk
Coordinate_text = fileread(Coordinate_File);

% Add Metal and Halide symbols
if strcmp(Settings.CoordType,'gro')
    Met = pad(Metal,2,'left');
    Hal = pad(Halide,2,'left');
elseif strcmp(Settings.CoordType,'g96')
    Met = pad(Metal,2,'right');
    Hal = pad(Halide,2,'right');
end
Coordinate_text = strrep(Coordinate_text,'##MET##',Met);
Coordinate_text = strrep(Coordinate_text,'##HAL##',Hal);

% Calculate size of supercell
N_Supercell = ceil((N_atoms/Settings.Geometry.N)^(1/3));

% Add number of unit cells to topology file
Topology_Temp_Struct = strrep(Topology_Temp_Struct,'##N##',num2str(N_Supercell));
Topology_Temp_Struct = strrep(Topology_Temp_Struct,'##GEOM##',Settings.Structure);

% Update File Base Name
FileBase = [Settings.Salt '_' Label '_' Model_Scaled '_Energy'];

% Update Directory
if Delete_output_files
    Settings.WorkDir = tempname;
else
    Settings.WorkDir = fullfile(Settings.project,...
        Settings.Project_Directory_Name,Settings.Salt,...
        Settings.Structure,Model_Scaled,OptTxt);
end

% Create directory if it does not exist
if ~exist(Settings.WorkDir,'dir')
    mkdir(Settings.WorkDir)
end

% Load topology template and MDP template
Topology_text = Topology_Temp_Struct;
Settings.MinMDP.Temp_Model = Settings.MinMDP.Temp_Ions;

% Update Topology and MDP files
Settings.MinMDP.Temp_Model = strrep(Settings.MinMDP.Temp_Model,'##COULOMB##',pad('PME',18));
Settings.MinMDP.Temp_Model = strrep(Settings.MinMDP.Temp_Model,'##FOURIER##',pad(num2str(Settings.MinMDP.Fourier_Spacing),18));
Settings.MinMDP.Temp_Model = strrep(Settings.MinMDP.Temp_Model,'##PMEORDER##',pad(num2str(Settings.MinMDP.PME_Order),18));
Settings.MinMDP.Temp_Model = strrep(Settings.MinMDP.Temp_Model,'##EWALDTOL##',pad(num2str(Settings.MinMDP.Ewald_rtol),18));

if Table_Req
    
    % Define the function type as 1 (needed for custom functions)
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##NBFUNC##','1');

    % Define the combination rules (Lorenz-berthelot)
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##COMBR##','1');

    % Define all the repulsive parameters as 1.0 (already included in potentials)
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METMETA##','1.0');
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##HALHALA##','1.0');
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METHALA##','1.0');
    
    % Generate tables
    TableName = [Settings.Salt '_' Model_Scaled '_Table'];
    [Settings.TableFile_MX,C6] = MakeTables(Settings,'MDP_Minimize',true,...
    	'TableName',TableName);
    
    % Dispersion coefficients
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METMETC##',num2str(C6.MM,'%.10e'));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##HALHALC##',num2str(C6.XX,'%.10e'));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METHALC##',num2str(C6.MX,'%.10e'));

    % Modify the MDP file
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##VDWTYPE##',pad('user',18));
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##VDWMOD##',pad('none',18));
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##CUTOFF##',pad('group',18));
    Settings.MDP_Template = regexprep(Settings.MDP_Template,'ewald-rtol-lj.+?\n','');
    Settings.MDP_Template = regexprep(Settings.MDP_Template,'lj-pme-comb-rule.+?\n','');
    Settings.MDP_Template = regexprep(Settings.MDP_Template,'verlet-buffer-tolerance.+?\n','');
elseif contains(Settings.Theory,'LJ')
    
    Settings.TableFile_MX = '';

    % Definte the function type as 1 (LJ)
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##NBFUNC##','1');

    % Define the combination rules (Lorenz-berthelot in sigma-epsilon form)
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##COMBR##','2');

    % Add parameters to topology text
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METMETC##',pad(num2str(Settings.S.S.MM,'%10.8e'),10));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##HALHALC##',pad(num2str(Settings.S.S.XX,'%10.8e'),10));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METHALC##',pad(num2str(Settings.S.S.MX,'%10.8e'),10));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METMETA##',num2str(Settings.S.E.MM,'%10.8e'));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##HALHALA##',num2str(Settings.S.E.XX,'%10.8e'));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METHALA##',num2str(Settings.S.E.MX,'%10.8e'));

    % Modify the MDP file
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##VDWTYPE##',pad(Settings.MinMDP.VDWType,18));
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##VDWMOD##',pad(Settings.MinMDP.vdw_modifier,18));
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##CUTOFF##',pad(Settings.MinMDP.CutOffScheme,18));
    Settings.MDP_Template = regexprep(Settings.MDP_Template,'energygrp-table.+?\n','');
    Settings.MDP_Template = regexprep(Settings.MDP_Template,'ewald-rtol-lj.+?\n','');
    Settings.MDP_Template = regexprep(Settings.MDP_Template,'lj-pme-comb-rule.+?\n','');
    
    % Add in Verlet Settings
    if strcmp(Settings.MinMDP.CutOffScheme,'Verlet')
        Settings.MDP_Template = strrep(Settings.MDP_Template,'##VerletBT##',pad(num2str(Settings.MinMDP.VerletBT),18));
    else
        Settings.MDP_Template = regexprep(Settings.MDP_Template,'verlet-buffer-tolerance.+?\n','');
    end
    
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

    Settings.TableFile_MX = '';

    % Definte the function type as 1 (LJ)
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##NBFUNC##','1');

    % Define the combination rules (Lorenz-berthelot in sigma-epsilon form)
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##COMBR##','2');

    % Get JC parameters
    [MX_JC_Param,MM_JC_Param,XX_JC_Param] = JC_Potential_Parameters(Settings);
    
    % Add parameters to topology text
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METMETC##',pad(num2str(MM_JC_Param.sigma,'%10.8e'),10));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##HALHALC##',pad(num2str(XX_JC_Param.sigma,'%10.8e'),10));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METHALC##',pad(num2str(MX_JC_Param.sigma,'%10.8e'),10));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METMETA##',num2str(MM_JC_Param.epsilon,'%10.8e'));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##HALHALA##',num2str(XX_JC_Param.epsilon,'%10.8e'));
    Settings.Topology_Template = strrep(Settings.Topology_Template,'##METHALA##',num2str(MX_JC_Param.epsilon,'%10.8e'));

    % Modify the MDP file
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##VDWTYPE##',pad(Settings.MinMDP.VDWType,18));
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##VDWMOD##',pad(Settings.MinMDP.vdw_modifier,18));
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##CUTOFF##',pad(Settings.MinMDP.CutOffScheme,18));
    Settings.MDP_Template = regexprep(Settings.MDP_Template,'energygrp-table.+?\n','');
    Settings.MDP_Template = regexprep(Settings.MDP_Template,'ewald-rtol-lj.+?\n','');
    Settings.MDP_Template = regexprep(Settings.MDP_Template,'lj-pme-comb-rule.+?\n','');
    
    % Add in Verlet Settings
    if strcmp(Settings.MinMDP.CutOffScheme,'Verlet')
        Settings.MDP_Template = strrep(Settings.MDP_Template,'##VerletBT##',pad(num2str(Settings.MinMDP.VerletBT),18));
    else
        Settings.MDP_Template = regexprep(Settings.MDP_Template,'verlet-buffer-tolerance.+?\n','');
    end
    
else
    error(['Warning: Unknown model type: "' Settings.Theory '.'])
end

if Settings.MDP.Disp_Correction
    Settings.MDP_Template = [Settings.MDP_Template newline newline ...
        '; Long-range dispersion correction' newline ...
        'DispCorr                 = EnerPres          ; apply long range dispersion corrections for Energy and pressure'];
end

if Settings.MinMDP.Verbose
    TotalTimer = tic;
    display_option = 'iter';
    pltfcn = 'psplotbestf';
    disp(['Calculating ' Settings.Salt ' ' Settings.Structure ' ' Model_Scaled ' Energy...'])
else
    display_option = 'none';
    pltfcn = [];
end

if ispc
    OptDirUnix = windows2unix(Settings.WorkDir);
    rm_command = ['wsl find ' OptDirUnix ' -iname "#*#" ^| xargs rm -f'];
else
    rm_command = ['find ' Settings.WorkDir ' -iname "#*#" | xargs rm -f'];
end
[~,~] = system(rm_command);

% Calculate with gromacs
OptimizationLoop(Settings,'Extra_Properties',Extra_Properties)
Output = GrabEnergy(Settings.WorkDir,FileBase,'Extra_Properties',Extra_Properties);
Output.Coulomb = Output.CoulSR + Output.CoulLR;

% Record output into object
if Settings.Geometry.Conv
    Output.a = Settings.Geometry.a;
    Output.b = Settings.Geometry.b;
    Output.c = Settings.Geometry.c;
else
    switch Settings.Structure
        case 'Rocksalt'
            Output.a = Settings.Geometry.a*sqrt(2);
            Output.b = Settings.Geometry.b*sqrt(2);
            Output.c = Settings.Geometry.c*sqrt(2);
        case 'Sphalerite'
            Output.a = Settings.Geometry.a*sqrt(2);
            Output.b = Settings.Geometry.b*sqrt(2);
            Output.c = Settings.Geometry.c*sqrt(2);
        case 'FiveFive'
            Output.a = Settings.Geometry.c;
            Output.b = (3/sqrt(3))*Settings.Geometry.a;
            Output.c = Settings.Geometry.b;
        otherwise
            Output.a = Settings.Geometry.a;
            Output.b = Settings.Geometry.b;
            Output.c = Settings.Geometry.c;
    end
end

if Settings.MinMDP.Verbose
    Settings.Geometry.a = Output.a;
    Settings.Geometry.b = Output.b;
    Settings.Geometry.c = Output.c;
    Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
    disp(['Completed: ' Settings.Salt ' ' Settings.Structure ' ' Model_Scaled ' Single Point Energy Calculation. Time ' Telap])
    disp(['Energy is ' num2str(Output.Pot,'%4.10f') ' kJ/mol'])
    DOF = {'a' 'b' 'c'};
    for Didx = 1:3
        disp(['Lattice Parameter ' DOF{Didx} ' = ' ...
            num2str(Settings.Geometry.(DOF{Didx}),'%2.8f') ' ' char(0197) '.']);
    end
    disp(['Fractional Coordinates for ' Metal ': ']);
    disp(num2str(Settings.Geometry.FC_Metal(:,:),'%2.8f '))
    disp(['Fractional Coordinates for ' Halide ': ']);
    disp(num2str(Settings.Geometry.FC_Halide(:,:),'%2.8f '))
end


% Delete the calculation directory
if Delete_output_files
    rmdir(Settings.WorkDir,'s')
end

% Return environmental variables back if changed, turn off parallel pool
if Settings.MinMDP.Parallel_Min
    setenv('OMP_NUM_THREADS',env.OMP_NUM_THREADS);
    setenv('GMX_PME_NUM_THREADS',env.GMX_PME_NUM_THREADS);
    setenv('GMX_PME_NTHREADS',env.GMX_PME_NTHREADS);
    setenv('GMX_OPENMP_MAX_THREADS',env.GMX_OPENMP_MAX_THREADS);
    setenv('KMP_AFFINITY',env.KMP_AFFINITY);
    delete(gcp('nocreate'));
end

end
