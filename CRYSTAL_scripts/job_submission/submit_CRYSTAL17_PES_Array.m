%% Generates crystal structure inputs for CRYSTAL17, plus batch job submission scripts
% Also capable of generating single ion or atom inputs for Alkali Metals or Halides
% submit_CRYSTAL17_PES_Array

%% LDA type functionals
% SVWN (SAME AS LSDA in GAUSSIAN)
 
%% GGA type functionals
% BLYP 
% PBEXC 
% PBESOLXC 
% SOGGAXC
% B97
 
%% Global Hybrid functionals
% B3PW
% B3LYP
% PBE0
% PBESOL0
% WC1LYP
% B97H
% PBE0-13
% PW1PW
 
%% Range-Separated Hybrid functionals
% HSE06
% HSEsol
% SC-BLYP
% HISS
% RSHXLDA
% LC-wPBE
% LC-wPBEsol
% LC-wBLYP
% wB97
% wB97X
% LC-BLYP
% CAM-B3LYP
 
%% meta-GGA functionals
% M06L
% M05
% M052X
% M06
% M062X
% M06HF
 
%% Double Hybrids (Not usable without CRYSCOR)
% B2PLYP
% B2GPPLYP
% mPW2PLYP

%% PURE EXCHANGE FUNCTIONALS
% LDA functionals
% LDA
% VBH
% GGA functionals
% BECKE
% mPW91
% PBE
% PBESOL
% PWGGA
% SOGGA
% WCGGA

%% PURE CORRELATION FUNCTIONALS
% LDA functionals
% PWLSD
% PZ
% VBH
% VWN
% LYP
% P86
% PBE
% PBESOL
% PWGGA
% WL

% Note that not all functionals available in CRYSTAL have been parametrized for D3(BJ). A list of
% available D3(BJ) dispersion corrected DFT methods follows: BLYP, PBE, B97, B3LYP, PBE0,
% PW1PW, M06, HSE06, HSEsol, LC-wPBE, PBESOL, CAM-B3LYP.

%% GLOST Settings (for parallel execution of jobs)
Glost_On = false; % Not compatible with restart jobs function (only works on cedar/graham?)
Glost_Name = 'ALL_NoDisp_FULL'; % Name your GLOST job here. Only used if Glost_On is true
Glost_Nodes = 1; % How many nodes to request

%% Vibration and QHA Thermodynamics Settings
QHA = false; % Quasi-Harominc Approximation. Do not combine with Calc_Vib. Used for calculation of volume dependant properties such as thermal expansion
QHA_Step = 6; % (percentrage) A volume range is defined from a -s% compression to a +2s% expansion (default 3%)
QHA_Points = 7; % Number of volumes to query using QHA: 4, 7, or 13 (default 4)
QHA_RANGE = [0.97 1.15 7]; % Set custom range here: [Min vol fraction   Max vol fraction   number of points] overrides QHA_Step and QHA_Points. Leave empty to not use

Calc_Vib = false; % Calculate vibration switch, automatically find lowest energy fully optimized configuration from data for input
T = '101 0 1000'; % Temperature (K): number of points, lower bound, upper bound
P = '101 0 10000'; % Pressure (MPa): number of points, lower bound, upper bound. Not used in QHA
Phonon_Dispersion = true; % Activates phonon dispersion keyword
Supercell = 2; % Used when Phonon dispersion is true or when QHA is on, creates supercell
PDOS = true; % Output the phonon density of states (ignored with QHA)
PDOS_Range = [1000 1001]; % Range of output density of states from 0 to first input. Second input = number of bins.
PDOS_projected = true; % Enables calculation of the Phonon density of states for each atom separately

%% Unit Cell Optimization Settings
Cell_optimization = true; % Switch to turn on optimization runs
opt_type = 'CELLONLY'; % 'FULLOPTG' = Optimize both unit cell vectors and atom positions (fixed symmetry), 'ATOMONLY' = Atoms only, 'CELLONLY' = cell only, 'ITATOCEL' = iteratively cell, atoms, cell, atoms ...
Internal_Coordinates = false; % Run geometry optimization in redundant internal coordinates when true (DOES NOT WORK WITH CELLONLY)
opt_conv = 8; % Convergence critera for optimization
max_cycle_opt = 500; % Maximum number of steps allowed for cell optimization
NUMGRCEL = false; % If true, compute numerical gradients (computationaly costly)
TOLDEG = '0.0003'; %'0.00003' Convergence criterion on the Root-Mean-Squared of the gradient between optimization steps
TOLDEX = '0.0012'; %'0.00012' Convergence criterion on the Root-Mean-Squared of the displacement between optimization steps
External_P = 0; % Adds external pressure in atmospheres. Set to 0 to ignore this option. Does not work with Vibrational analysis or QHA, even for pre-opt geometry
Updating_Scheme = 'BFGS'; % Default is 'BFGS', other options: 'BERNY' or 'POWELL'

%% Geometry Search Options
AutoFindCellParams = true; % When true, automatically start from fully optimized cell parameters if available.
AutoFindDataTypes = 1; % 0 = single point, 1 = cell optimized, 2 = full optimized, 3 = atom optimized only, 4 = full optimized in space group P1, 5 = full optimized at 1 atm pressure, 
Relax_Basis_Set = 'pob-TZVP'; % Leave blank to constrain the auto find cell parameters to the chosen basis set. Otherwise allows search for optimized cell parameters with this basis set
Relax_XC = ''; % Leave blank to constrain the auto find cell parameters to the chosen basis set, Otherwise allows search for optimized cell parameters with this XC functional.
Relax_Dispersion = 'Off'; % Set as 'Off' to constrain the auto find cell parameters to the chosen dispersion. Otherwise allows search for optimized cell parameters with this dispersion
Relax_gCP = 'Off'; %  Set as 'Off' to constrain the auto find cell parameters to the chosen gCP settings. Otherwise allows search for optimized cell parameters with this gCP setting
% If you wish to use a different basis set for metal vs halide then use the
% following format: 'M[Metal Basis set]_H[Halide Basis Set]'

% Option to skip jobs that already have data associated with them (set to
% false to rerun such jobs)
Skip_Completed = true;
Skip_Running = true; % Skip jobs that are detected to be running, only works on orcinus currently.
Stop_Running = false; % Ends jobs that are detected to be running and replaces with a new job. (Only works if Skip_Running is on)

%% Theory and Correction Options
% 'HF' 'HSE06' 'HSEsol' 'M06' 'PBE' 'PW1PW' 'B3LYP'
% 'HF' 'HSE06' 'HSEsol' 'LSDA' 'M06' 'M062X' 'M06L' 'PBE' 'PW1PW' 'B3LYP'
Theories = {'HF'};
Dispersion = ''; % Set to 'D3' (includes BJ damping) / 'D3TB' (3 body terms included) / 'D2' / '' (no dispersion);
Counterpoise = ''; % Set to 'gCP' to enable the geometrical counterpoise method for basis set superposition error estimation 

%% Structure and Salt Options
Salt_or_Ion_Types = {'LiI'}; % Either Lithium Halide salt OR a single Lithium/halide ion
Structures = {'Rocksalt' 'Wurtzite' 'FiveFive' 'NiAs' 'Sphalerite'}; % Available: {'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive' 'All' 'Pair' 'Atom' 'Ion'}
Initialize_As_Ions = true; % When True, initial basis set will be charge separated into positive and negative ions
P1_Symmetry = false; % When false, use salt geometry with the full symmetry of the space group. When true, use space group P1.

%% Batch Job Specification Settings
num_links = 1; % Link together batch jobs with restart function
Restart_From_Previous = false; % Only use when previous job did not finish, or you would like to restart a job with different parameters (does not work with array jobs)
Restart_WF = false; % Restarts the SCF using GUESSP when true. Requires a fort.20/fort.9/fort.79 or (filename).f20 file present in the job directory
submit = true; % Run the jobs, or just generate submission script?
submission_folder = '/home/scheiber/project/CRYSTAL_LiHalides'; % Main folder
hours = 120; % Max time for job (hours)
mins = 0; % Max time for job (minutes)
nCores = 1; % Number of cores to request for calculation (currently limited to 1)
nTasks_per_Node = 1; % Cores per node to request
mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', use '3G' for Cedar)
template_file = 'Crystal17.template'; % Name of template file

%% SCF AND CONVERGENCE SETTINGS
conv = 11; % Convergence criteria
PM_net = 15; %12 Number of K-points along one dimension in Pack-Monkhorst Net (meaniningless for ions)
G_net = 30; %24 Number of K-points along one dimension in Gilat Net (meaniningless for ions)
integration = [15 15 15 20 30]; %[8 8 8 8 16]; %[15 15 15 20 30]; % 5 Integration parameters
Save_WF = true; % Save WF every 2 iterations of the SCF process when true (saved as fort.79 in temp folder)
max_SCF = 500; % Maximum number of SCF cycles
Latvec = 12000; %4500 % Maximum size of the cluster of classified lattice vectors (default 3500, meaningless for ions)
ILAsize = 12000; % ILA is an array containing a list of contributions to be computed in the evaluation of the Coulomb series. (default 6000)
LIMBEK = 900; %500 % Size of local arrays for integration weights [default 400] (for DFT only)
GridSize = 'XXLGRID'; % DFT integration grid size. Default is 'XLGRID'
BIPOSIZE = 9000000; % Size of buffer for Coulomb integrals bipolar expansion Default 4000000
EXCHSIZE = 9000000; % Size of buffer for Exchange integrals bipolar expansion Default 4000000

%% Convergence Acceleration techniques
Anderson = false; % Default = false; Do not use with Broyden or DIIS.
Broyden = [0.0001 50 2]; %[0.0001 50 2]; % [W0 IMIX ISTART]. Suggested: [0.0001 50 2]. Leave empty to not use. Do not use with Anderson or DIIS. Overwrites Fock_mixing once active.
DIIS = 'NODIIS'; % Default = 'DIIS'. One of: 'DIIS' , 'DIISALLK' , 'SLOSHING' , 'NODIIS'.
HISTDIIS = []; % Limit of DIIS history. Blank for no limit.
THREDIIS = []; % Postpones the activation of DIIS by performing Fock mixing until the difference in energy is below this value.
THRKDIIS = []; % All the Fock matrices having a weight lower than this are discarded from DIIS history
Fock_mixing = 50; % 50% By default, DIIS is used.
Level_shift = []; %[3 1]; % Leave empty to turn off. Warning, do not combine with DIIS

%% BASIS SET SETTINGS
% Available: 'pob-TZVP' / 'pob-TZVP-rev2' (Li, Na, K, F, Cl, Br ONLY) /
% 'pob-TZVPP' (Li, Na, F, Cl ONLY) / '7-311G' / 'STO-3G' / 'def2-TZVPD' /
% 'cc-pV5Z' (Li, Na, F, Cl, BR ONLY) / 'aug-cc-pV5Z' (Li, Na, F, Cl, BR ONLY) /
% 'aug-cc-pV5ZD' (Li, Na, F, Cl, BR ONLY) / 'def2-TZVP' (Li, Na, F, Cl, BR ONLY) /
% 'def2-QZVP' (Li, Na, F, Cl, BR, I ONLY) / 'def2-QZVPD' (Li, Na, F, Cl, BR, I ONLY) /
% 'def2-SVPD' (Li, Na, F, Cl, BR ONLY) /
% 'ECP28MDF-VTZ' (I ONLY) / 'ECP28MDF-AVTZ' (I ONLY) /
% 'ECP28MDF-VQZ' (I ONLY) / 'ECP28MDF-AVQZ' (I ONLY) /
% 'ECP28MDF-V5Z' (I ONLY) / 'ECP28MDF-AV5Z' (I ONLY)
Metal_Basis_set = {'pob-TZVP'};
Halide_Basis_set = {'def2-QZVP'};

% Switch to activate exponent scaling (for pob_TZVP / pob_TZVPP / pob-TZVP-rev2 basis sets only)
Scale_exponents = false;
Scale_M1 = 2.092403731668759; %2.092403731668759;
Scale_M2 = 1; 
Scale_H1 = 1;
Scale_H2 = 1;
Scale_H3 = 1;

% Prune exponents (for pob_TZVP / pob_TZVPP / pob-TZVP-rev2 basis sets only)
Prune_M1 = false; % Most Dif: Li      Mid Dif:         Least Dif:
Prune_M2 = false; % Most Dif:         Mid Dif:         Least Dif: Li
Prune_H1 = false; % Most Dif: Cl, Br  Mid Dif: F       Least Dif: I
Prune_H2 = false; % Most Dif: F, I    Mid Dif: Cl, Br  Least Dif: 
Prune_H3 = false; % Most Dif:         Mid Dif: I       Least Dif: F, Cl, Br

% Add exponents (not compatible with '7-311G' / 'STO-3G' / 'def2-TZVPD' basis sets), '' for nothing.
Add_M = '';
Add_H = ''; %['0 0 1 0 1.0' newline '  0.1500000000      1.000000000000'];
if ~isempty(Add_M)
    Aug_Label_M = 'D'; % Label for augmented metal basis set
else
    Aug_Label_M = '';
end
if ~isempty(Add_H)
    Aug_Label_H = 'D'; % Label for augmented halide basis set
else
    Aug_Label_H = '';
end

%% BSSE with Counterpoise Settings
BSSE_id = ''; % 'Metal', 'Halide', or '' (no BSSE correction);
BSSE_Nstars = 2; % Number of stars of ghost atoms included in calculation USE MORE WITH BetaBeO
BSSE_Rmax = 45; % maximum distance for ghost atoms in atomic units
findBSSEparameters = false; % Instead of chosen cell parameters, find them from data

%% Add point charges in a cluster around central unit cell (BSSE ONLY)
N_Cluster = 0; % an N_Cluster x N_Cluster x N_Cluster placement of charges will be produced. Use N_Cluster = 0 to turn off.

%% CRYSTAL DEFAULT GEOMETRIES
% For crystals choose a single number or range. For optimization, initial
% conditions are found from lowest energy known structure
a_BetaBeO = 7.3607;
c_over_a_BetaBeO = 4.2663/7.3607; % Theoretical c/a = 1/sqrt(3)
BBO_MetCoord = '0.3360 0.6640 0.0000'; %'0.336 -0.336  0.000';
BBO_HalCoord = '0.3100 0.3100 0.0000'; %'0.310  0.310  0.000';

a_CsCl = 3.24902899;

a_FiveFive = 5.1300;         
b_over_a_FiveFive = 3/2; % Theoretical b/a = 3/2
c_over_a_FiveFive = cosd(30); % Theoretical c/a = cosd(30)
FF_MetCoord = '0.25 0.166666666666667 0.5';
FF_HalCoord = '0.25 0.333333333333333 0.0';

a_NiAs = 4.06113217;
c_over_a_NiAs = 6.76014149/4.06113217; %1.7316 LiF LSDA, 1.7275 LiF PBE

a_Rocksalt = 3.7814; %5.4903;

a_Sphalerite = 6.0456;

a_Wurtzite = 3.1499;
c_over_a_Wurtzite = sqrt(8/3); % Perfect wurtzite c/a = sqrt(8/3);
Wurtz_metcoord = '0.33333333333 0.666666666667 0.000';
Wurtz_halcoord = '0.33333333333 0.666666666667 0.375';

BL_pair = 1.0:0.1:6.0; % for pair calculations: the bond length

%% Begin code
if ispc % for testing
    addpath('C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\submission')
    addpath('C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\analysis')
    home = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase'; % PC
    server = 'Unbearabull';
    Glost_On = false;
elseif isunix
    home = '/home/scheiber/ThesisCodeBase'; % Cedar/Graham/orcinus
    server = getenv('HOSTNAME');
    if strcmpi(server(1:3),'sea') % GLOST not active on orcinus
        Glost_On = false;
    elseif strcmpi(server(1:3),'ced')
        if strcmp(mempernode,'0')
            mempernode = 'MaxMemPerCPU';
        end
        Account = 'rrg-patey-ad';
        Glost_Cores_Per_Node = 48;
    else
        if strcmp(mempernode,'0')
            mempernode = 'MaxMemPerCPU';
        end
        Account = 'def-patey';
        Glost_Cores_Per_Node = 32;
    end
else
    home = input('Please input home bin directory.\n','s');
end

% Some error catches
if findBSSEparameters && isempty(BSSE_id)
    warning('findBSSEparameters should only be active when using BSSE: Disabling.')
    findBSSEparameters = false;
end
if Cell_optimization && findBSSEparameters
    error('Cannot do Cell Optimization and find BSSE Parameters simultaneously.')
elseif Cell_optimization && ~isempty(BSSE_id)
    error('Cannot do Cell Optimization and BSSE Calculation simultaneously.')
elseif findBSSEparameters && (Calc_Vib || QHA)
    error('Cannot do Frequency Analysis and find BSSE Parameters simultaneously.')
elseif QHA && Calc_Vib
    error('Choose only one of QHA or Calc_Vib.')
elseif (QHA || Calc_Vib) && (AutoFindDataTypes ~= 2)
    error('Always run QHA or vibrational analysis jobs from fully minimized structures.')
end

% Crystal types
if strcmp('All',Structures)
    Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive' 'Pair'};
else
    if ~iscell(Structures)
        Structures = {Structures};
    end
end

% Salt/ion types
if ~iscell(Salt_or_Ion_Types)
    Salt_or_Ion_Types = {Salt_or_Ion_Types};
end

% Add newlines to extra exponent text
if ~isempty(Add_M)
    Add_M = [newline Add_M];
    NumM = sum(Add_M == 10)/2;
else
    NumM = 0;
end
if ~isempty(Add_H)
    Add_H = [newline Add_H];
    NumH = sum(Add_H == 10)/2;
else
    NumH = 0;
end

% Move to outer directory
if ispc
    Maindir = pwd;
elseif isunix
    Maindir = submission_folder;
    cd(Maindir)
end

% Restart job functions
if Restart_From_Previous
    Restart_txt = '1';
else
    Restart_txt = '0';
end
if Restart_WF
    Restart_WF_txt = '1';
else
    Restart_WF_txt = '0';
end

if P1_Symmetry
    SymmetryMod = '_SG1';
    OutSymTxt = ' Space Group 1';
else
    SymmetryMod = '';
    OutSymTxt = '';
end

% External pressure settings
if External_P == 0 || Calc_Vib || QHA
    DirPTxt = '';
    OutPTxt = ' (P = 0 Atm)';
    ExtP_AU = 0;
    EXTPRESS = '';
else
    DirPTxt = ['_P' sprintf('%1.2f',External_P)];
    OutPTxt = [' (P = ' sprintf('%1.2f',External_P) ' Atm)'];
    
    AU_Per_GPa = 1/(29421.912); % (Ha/Bohr^3) / GPa
    GPa_Per_Atm = 0.000101325; % GPa / Atm
    AU_Per_Atm = AU_Per_GPa*GPa_Per_Atm; % 1 atm per Ha/Bohr^3
    ExtP_AU = External_P*AU_Per_Atm; % Convertion to Ha/Bohr^3 for input
    
    EXTPRESS = ['EXTPRESS' newline num2str(ExtP_AU,'%1.8E') newline];
end

%% Loop through Salt/ion types
for id=1:length(Salt_or_Ion_Types)
    Current_Salt = Salt_or_Ion_Types{id};
    
    %% Loop through crystal types
    for idx = 1:length(Structures)
        Structure = Structures{idx};
        
        % Change directory to crystal type
        if strcmp('Ion',Structure)
            BSSE = [];
            % Check to make sure it is acceptable ion type
            [I1,I2] = Separate_Metal_Halide(Current_Salt);
            if isempty(I1) && isempty(I2)
                warning([Current_Salt ' is not a known single-ion type,' ...
                   ' excluding from submission.']);
               continue
            end

            Ion_dir = [Current_Salt '_Ion'];
            if ~exist(Ion_dir,'dir')
                mkdir(Ion_dir);
            end
            cd(Ion_dir)
        elseif strcmp('Atom',Structure)
            BSSE = [];
            % Check to make sure it is acceptable atom type
            [I1,I2] = Separate_Metal_Halide(Current_Salt);
            if isempty(I1) && isempty(I2)
                warning([Current_Salt ' is not a known single-atom type,' ...
                   ' excluding from submission.']);
               continue
            end

            Atom_dir = [Current_Salt '_Atom'];
            if ~exist(Atom_dir,'dir')
                mkdir(Atom_dir);
            end
            cd(Atom_dir)
        else
            % Check to make sure it is acceptable salt type
            [Metal,Halide] = Separate_Metal_Halide(Current_Salt);
            if isempty(Metal) && isempty(Halide)
                warning([Current_Salt ' is not a known salt type,' ...
                   ' excluding from submission.']);
               continue
            elseif ( isempty(Metal) || isempty(Halide) ) && ~strcmp(Structure,'Pair')
                warning([Current_Salt ' can only be run as a pair potential,' ...
                   ' excluding ' Current_Salt ' ' Structure ' from submission.']);
               continue
            % For Pair potentials of Li-Li or X-X
            elseif isempty(Metal)
                Metal = Halide;
                Double = 'Halide';
            elseif isempty(Halide)
                Halide = Metal;
                Double = 'Metal';
            end
            
            if isempty(BSSE_id)
                BSSE = [];
                Salt_Crystal_dir = [Current_Salt filesep ...
                    Structure];
                % Make directory
                [~, ~, ~] = mkdir(Salt_Crystal_dir);
                cd(Salt_Crystal_dir)
                
            elseif strcmp('Metal',BSSE_id)
                if strcmp(Structure,'Pair')
                    BSSE = [];
                else
                    BSSE = [1 BSSE_Nstars BSSE_Rmax];
                end
                
                % Move to BSSE directory for specified metal
                BSSE_dir = ['BSSE_Estimation' filesep Metal];
                if ~exist(BSSE_dir,'dir')
                    mkdir(BSSE_dir);
                end
                cd(BSSE_dir)
                
            elseif strcmp('Halide',BSSE_id)
                if strcmp(Structure,'Pair')
                    BSSE = [];
                elseif contained_in_cell(Structure,{'Rocksalt' 'Sphalerite' 'CsCl'})
                    BSSE = [2 BSSE_Nstars BSSE_Rmax];
                elseif contained_in_cell(Structure,{'Wurtzite' 'NiAs'})
                    BSSE = [3 BSSE_Nstars BSSE_Rmax];
                elseif contained_in_cell(Structure,{'BetaBeO' 'FiveFive'})
                    BSSE = [5 BSSE_Nstars BSSE_Rmax];
                end
                
                % Move to BSSE directory for specified Halide
                BSSE_dir = ['BSSE_Estimation' filesep Halide];
                if ~exist(BSSE_dir,'dir')
                    mkdir(BSSE_dir);
                end
                cd(BSSE_dir)
            else
                error('Unknown input for "BSSE" variable.')
            end
        end

        %% Structure specific unit cell stuff   
        if strcmp('Wurtzite',Structure)
            a = a_Wurtzite;
            if P1_Symmetry
                MFC = textscan(Wurtz_metcoord,'%f %f %f',1);
                HFC = textscan(Wurtz_halcoord,'%f %f %f',1);
                Crystal_Info = P1_Setting(Structure,[MFC{:}],[HFC{:}]);
                Space_group = '1';
                Cell_param = ['##APAR## ##BPAR## ##CPAR## ' num2str(Crystal_Info.alpha) ' ' ...
                    num2str(Crystal_Info.beta) ' ' num2str(Crystal_Info.gamma) newline];
                Coordinates = [num2str(Crystal_Info.N) newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(1,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(2,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(1,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(2,:),12),' +',' ')];
            else
                Space_group = '186';
                Cell_param = ['##APAR## ##CPAR##' newline];
                Coordinates = ['2' newline ...
                    '##Metal## ' Wurtz_metcoord newline ...
                    '##Halide## ' Wurtz_halcoord]; 
            end
            
            labend = '_W';
        elseif strcmp('Rocksalt',Structure)
            a = a_Rocksalt;
            if P1_Symmetry
                MFC = [0.0 0.0 0.0];
                HFC = [0.5 0.0 0.0];
                Crystal_Info = P1_Setting(Structure,MFC,HFC);
                Space_group = '1';
                Cell_param = ['##APAR## ##BPAR## ##CPAR## ' num2str(Crystal_Info.alpha) ' ' ...
                    num2str(Crystal_Info.beta) ' ' num2str(Crystal_Info.gamma) newline];
                Coordinates = [num2str(Crystal_Info.N) newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(1,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(2,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(1,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(2,:),12),' +',' ')];
            else
                Space_group = '225';
                Cell_param = ['##APAR##' newline];
                Coordinates = ['2' newline ...
                    '##Metal## 0.0 0.0 0.0' newline ...
                    '##Halide## 0.5 0.0 0.0'];
            end
            labend = '_R';
        elseif strcmp('Sphalerite',Structure)
            a = a_Sphalerite;
            if P1_Symmetry
                MFC = [0.0 0.0 0.0];
                HFC = [0.25 0.25 0.25];
                Crystal_Info = P1_Setting(Structure,MFC,HFC);
                Space_group = '1';
                Cell_param = ['##APAR## ##BPAR## ##CPAR## ' num2str(Crystal_Info.alpha) ' ' ...
                    num2str(Crystal_Info.beta) ' ' num2str(Crystal_Info.gamma) newline];

                Coordinates = [num2str(Crystal_Info.N) newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(1,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(2,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(3,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(4,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(1,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(2,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(3,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(4,:),12),' +',' ')];
            else
                Space_group = '216';
                Cell_param = ['##APAR##' newline];
                Coordinates = ['2' newline ...
                    '##Metal## 0.0 0.0 0.0' newline ...
                    '##Halide## 0.25 0.25 0.25'];
            end
            labend = '_S';
        elseif strcmp('CsCl',Structure)
            a = a_CsCl;
            if P1_Symmetry
                MFC = [0.0 0.0 0.0];
                HFC = [0.5 0.5 0.5];
                Crystal_Info = P1_Setting(Structure,MFC,HFC);
                Space_group = '1';
                Cell_param = ['##APAR## ##BPAR## ##CPAR## ' num2str(Crystal_Info.alpha) ' ' ...
                    num2str(Crystal_Info.beta) ' ' num2str(Crystal_Info.gamma) newline];
                Coordinates = [num2str(Crystal_Info.N) newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(1,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(1,:),12),' +',' ')];
            else
                Space_group = '221';
                Cell_param = ['##APAR##' newline];
                Coordinates = ['2' newline ...
                    '##Metal## 0.0 0.0 0.0' newline ...
                    '##Halide## 0.5 0.5 0.5'];
            end
            labend = '_C';
        elseif strcmp('NiAs',Structure)
            a = a_NiAs;
            if P1_Symmetry
                MFC = [0.0 0.0 0.0];
                HFC = [1/3 2/3 0.25];
                Crystal_Info = P1_Setting(Structure,MFC,HFC);
                Space_group = '1';
                Cell_param = ['##APAR## ##BPAR## ##CPAR## ' num2str(Crystal_Info.alpha) ' ' ...
                    num2str(Crystal_Info.beta) ' ' num2str(Crystal_Info.gamma) newline];
                Coordinates = [num2str(Crystal_Info.N) newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(1,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(2,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(1,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(2,:),12),' +',' ')];
            else
                Space_group = '194';
                Cell_param = ['##APAR## ##CPAR##' newline];
                Coordinates = ['2' newline ...
                    '##Metal## 0.0 0.0 0.0' newline ...
                    '##Halide## 0.33333333333 0.66666666667 0.25'];
            end
            labend = '_N';
        elseif strcmp('BetaBeO',Structure)
            a = a_BetaBeO;
            if P1_Symmetry
                MFC = textscan(BBO_MetCoord,'%f %f %f',1);
                HFC = textscan(BBO_HalCoord,'%f %f %f',1);
                MFC = [MFC{:}];
                HFC = [HFC{:}];
                Crystal_Info = P1_Setting(Structure,MFC,HFC);
                Space_group = '1';
                Cell_param = ['##APAR## ##BPAR## ##CPAR## ' num2str(Crystal_Info.alpha) ' ' ...
                    num2str(Crystal_Info.beta) ' ' num2str(Crystal_Info.gamma) newline];
                Coordinates = [num2str(Crystal_Info.N) newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(1,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(2,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(3,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(4,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(1,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(2,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(3,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(4,:),12),' +',' ')];
            else
                Space_group = '136';
                Cell_param = ['##APAR## ##CPAR##' newline];
                Coordinates = ['2' newline ...
                    '##Metal## ' BBO_MetCoord newline ...
                    '##Halide## ' BBO_HalCoord];
            end
            labend = '_B';
        elseif strcmp('FiveFive',Structure)
            a = a_FiveFive;
            if P1_Symmetry
                MFC = textscan(FF_MetCoord,'%f %f %f',1);
                HFC = textscan(FF_HalCoord,'%f %f %f',1);
                MFC = [MFC{:}];
                HFC = [HFC{:}];
                Crystal_Info = P1_Setting(Structure,MFC,HFC);
                Space_group = '1';
                Cell_param = ['##APAR## ##BPAR## ##CPAR## ' num2str(Crystal_Info.alpha) ' ' ...
                    num2str(Crystal_Info.beta) ' ' num2str(Crystal_Info.gamma) newline];
                Coordinates = [num2str(Crystal_Info.N) newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(1,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(2,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(3,:),12),' +',' ') newline ...
                    '##Metal## ' regexprep(num2str(Crystal_Info.FC_Metal(4,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(1,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(2,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(3,:),12),' +',' ') newline ...
                    '##Halide## ' regexprep(num2str(Crystal_Info.FC_Halide(4,:),12),' +',' ')];
            else
                Space_group = '58';
                Cell_param = ['##APAR## ##BPAR## ##CPAR##' newline];
                Coordinates = ['2' newline ...
                    '##Metal## ' FF_MetCoord newline ...
                    '##Halide## ' FF_HalCoord];
            end
            labend = '_F';
        elseif strcmp('Pair',Structure)
            a = BL_pair;
            Space_group = '1';
            labend = '_P';
            Periodicity = 'MOLECULE';
            Cell_param = '';
            Coordinates = ['2' newline ...
                '##Metal## 0.0 0.0 0.0' newline ...
                '##Halide## ##APAR## 0.0 0.0'];
            kpoints = '';
            Latvectext = '';
            
            % For Pair BSSE calculations, turn the metal or halide into a ghost
            if strcmp(BSSE_id,'Halide')
                Coordinates = strrep(Coordinates,'##Metal##','0');
            elseif strcmp(BSSE_id,'Metal')
                Coordinates = strrep(Coordinates,'##Halide##','0');
            end

            if Calc_Vib
                submit_Crystal = 'submit_CRYSTAL17.pl';
                params = '';
                display_text = [' Vibration Analysis' OutSymTxt];
                optTag = '_Vib';
            elseif QHA
                submit_Crystal = 'submit_CRYSTAL17.pl';
                params = '';
                display_text = [' Vibration Analysis' OutSymTxt];
                optTag = '_QHA';
            elseif Cell_optimization
                submit_Crystal = 'submit_CRYSTAL17.pl';
                params = '';               
                if strcmpi(opt_type,'FULLOPTG') || strcmpi(opt_type,'ITATOCEL')
                    display_text = [' Full Opt' OutSymTxt OutPTxt];
                elseif strcmpi(opt_type,'CELLONLY')
                    display_text = [' Cell Opt' OutSymTxt OutPTxt];
                elseif strcmpi(opt_type,'ATOMONLY')
                    display_text = [' Atom Opt' OutSymTxt OutPTxt];
                else
                    error(['Unknown Optimization Type: ' opt_type])
                end
                optTag = '_CellOpt';
            elseif AutoFindCellParams
                submit_Crystal = 'submit_CRYSTAL17.pl';
                params = '';
                optTag = '';
            else
                submit_Crystal = 'submit_CRYSTAL17_array.pl';
                params = [' ' regexprep(num2str(a.*100),' +',',')];
                display_text = [' ' num2str(min(a)) ' - ' num2str(max(a)) OutSymTxt];
                optTag = '';
            end
            % For M-M or X-X calculations
            if strcmp(Metal,Halide)
                if strcmp(Double,'Metal')
                    Coordinates = strrep(Coordinates,'##Halide##','##Metal##');
                elseif strcmp(Double,'Halide')
                    Coordinates = strrep(Coordinates,'##Metal##','##Halide##');
                end
            end
            
        elseif strcmp('Ion',Structure)
            display_text = ' Ion';
        elseif strcmp('Atom',Structure)
            display_text = ' Atom';
        else
            error(['Error: Unknown Crystal Structure Type: ' Structure]);
        end

        if contained_in_cell(Structure,{'Wurtzite' 'Rocksalt' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'})
            Periodicity = ['CRYSTAL' newline '0 0 0'];
            kpoints = ['SHRINK' newline num2str(PM_net) ' ' num2str(G_net) newline];
            if Calc_Vib && Cell_optimization
                submit_Crystal = 'submit_CRYSTAL17.pl';
                params = '';
                display_text = [' Vibration Analysis' OutSymTxt OutPTxt];
                optTag = '_Vib';
                
                Opt_Input = Cell_Optimization_Text('FULLOPTG',opt_conv,...
                    max_cycle_opt,NUMGRCEL,TOLDEG,TOLDEX,EXTPRESS,...
                    Internal_Coordinates,Updating_Scheme);
                Vib_Input = Vibrational_Analysis_Text(QHA,QHA_Step,QHA_Points,QHA_RANGE,T,P,Phonon_Dispersion,Supercell,PDOS,PDOS_Range,PDOS_projected,Opt_Input);
            elseif QHA && Cell_optimization
                submit_Crystal = 'submit_CRYSTAL17.pl';
                params = '';
                display_text = [' QHA Vib. Analysis' OutSymTxt OutPTxt];
                optTag = '_QHA';              
                
                Opt_Input = Cell_Optimization_Text('FULLOPTG',opt_conv,max_cycle_opt,...
                    NUMGRCEL,TOLDEG,TOLDEX,EXTPRESS,Internal_Coordinates,Updating_Scheme);
                Vib_Input = Vibrational_Analysis_Text(QHA,QHA_Step,QHA_Points,...
                    QHA_RANGE,T,P,Phonon_Dispersion,Supercell,PDOS,PDOS_Range,...
                    PDOS_projected,Opt_Input);
            elseif Calc_Vib
                submit_Crystal = 'submit_CRYSTAL17.pl';
                params = '';
                display_text = [' Vibration Analysis' OutSymTxt OutPTxt];
                optTag = '_Vib';
                
                Vib_Input = Vibrational_Analysis_Text(QHA,QHA_Step,QHA_Points,QHA_RANGE,T,P,Phonon_Dispersion,Supercell,PDOS,PDOS_Range,PDOS_projected,'');
            elseif QHA
                submit_Crystal = 'submit_CRYSTAL17.pl';
                params = '';
                display_text = [' QHA Vib. Analysis' OutSymTxt OutPTxt];
                optTag = '_QHA';
                
                Vib_Input = Vibrational_Analysis_Text(QHA,QHA_Step,QHA_Points,QHA_RANGE,T,P,Phonon_Dispersion,Supercell,PDOS,PDOS_Range,PDOS_projected,'');
            elseif Cell_optimization
                submit_Crystal = 'submit_CRYSTAL17.pl';
                params = '';
                if strcmpi(opt_type,'FULLOPTG') || strcmpi(opt_type,'ITATOCEL')
                    display_text = [' Full Opt' OutSymTxt OutPTxt];
                elseif strcmpi(opt_type,'CELLONLY')
                    display_text = [' Cell Opt' OutSymTxt OutPTxt];
                elseif strcmpi(opt_type,'ATOMONLY')
                    display_text = [' Atom Opt' OutSymTxt OutPTxt];
                else
                    error(['Incorrect Optimization Type: ' opt_type])
                end
                optTag = '_CellOpt';
            elseif AutoFindCellParams
                submit_Crystal = 'submit_CRYSTAL17.pl';
                params = '';
                optTag = '';
                Glost_On = false;
            else
                submit_Crystal = 'submit_CRYSTAL17_array.pl';
                params = [' ' regexprep(num2str(a.*100),' +',',')];
                optTag = '';
                Glost_On = false;
            end
            Latvectext = ['LATVEC' newline num2str(Latvec) newline];
        elseif contained_in_cell(Structure,{'Atom' 'Ion'})
            Periodicity = 'MOLECULE';
            Space_group = '1';
            Cell_param = '';
            Coordinates = ['1' newline ...
                '##Metal####Halide## 0.0 0.0 0.0'];
            labend = '';
            kpoints = '';
            submit_Crystal = 'submit_CRYSTAL17.pl';
            params = '';
            Latvectext = '';
            optTag = '';
        end

        label = [Current_Salt labend];

        %% Edit basis set
        % Remove scaling if turned off
        if ~Scale_exponents
            Scale_M1 = 1;
            Scale_M2 = 1;
            Scale_H1 = 1;
            Scale_H2 = 1;
            Scale_H3 = 1;
        end

        % Basis set superposition error
        if isempty(BSSE)
            BSSE_text = '';
        else
            % Beta BeO cannot be used with automatic ATOMBSSE routine
            if strcmpi(Structure,'BetaBeO')
                Periodicity = 'MOLECULE';
                Space_group = '1';
                Cell_param = '';
                kpoints = '';
                Latvectext = '';
                BSSE_text = '';
            else
                BSSE_text = ['ATOMBSSE' newline num2str(BSSE(1)) ' ' ...
                    num2str(BSSE(2)) ' ' num2str(BSSE(3)) newline];
            end
        end

        % Misc inputs
        workDir = pwd;
        nTasks = num2str(nCores);
        hours_calc = num2str(hours);
        mins_calc = num2str(mins);
        nCores_per_Node = num2str(nTasks_per_Node);
        integration_params = regexprep(num2str(integration),' +',' ');
        max_cycle = num2str(max_SCF);
        convergence = num2str(conv);
        ILA = num2str(ILAsize);
        Limbek = num2str(LIMBEK);
        links = num2str(num_links);
        biposize = num2str(BIPOSIZE);
        exchsize = num2str(EXCHSIZE);
        
        if Save_WF
            SAVEWF = ['SAVEWF' newline];
        else
            SAVEWF = '';
        end

        % DIIS acceleration
        if isempty(HISTDIIS)
            HISTDIIS_text = '';
        else
            HISTDIIS_text = [newline 'HISTDIIS' newline num2str(HISTDIIS)];
        end

        if isempty(THREDIIS)
            THREDIIS_text = '';
        else
            THREDIIS_text = [newline 'THREDIIS' newline num2str(THREDIIS)];
        end

        if isempty(THRKDIIS)
            THRKDIIS_text = '';
        else
            THRKDIIS_text = [newline 'THRKDIIS' newline num2str(THRKDIIS)];
        end
        DIIS_Text = [DIIS HISTDIIS_text THREDIIS_text THRKDIIS_text];

        % Anderson method
        if Anderson
            Anderson_text = ['ANDERSON' newline];
        else
            Anderson_text = ''; %#ok<*UNRCH>
        end

        % Broyden method
        if isempty(Broyden)
            Broyden_text = '';
        else
            Broyden_text = ['BROYDEN' newline regexprep(num2str(Broyden),' +',' ') newline];
        end

        % Fock mixing
        if isempty(Fock_mixing)
            FMIXING = '';
        else
            FMIXING = ['FMIXING' newline num2str(Fock_mixing) newline];
        end

        % Level shifting
        if isempty(Level_shift)
            LEVSHIFT = '';
        else
            LEVSHIFT = ['LEVSHIFT' newline regexprep(num2str(Level_shift),' +',' ') newline];
        end

        for i = 1:length(Theories)
            
            Theory = Theories{i};  

            % Use Theory Selection Subroutine to get Theory text
            [DFT,disptext,dispswitch,parameters,version] = ...
                Select_Theory_for_CRYSTAL17(Theory,Dispersion,Limbek,...
                Structure,GridSize);

            % Dispersion options
            if strcmp(Dispersion,'D3TB') && dispswitch
                dispersion_text = ['DFTD3' newline version newline ...
                    'ABC' newline parameters newline 'END' newline];
            elseif strcmp(Dispersion,'D2') && dispswitch
                dispersion_text = ['DFTD3' newline 'VERSION' newline '2' newline ...
                    'END' newline];
            elseif ~isempty(parameters) && dispswitch
                dispersion_text = ['DFTD3' newline version newline parameters ...
                    newline 'END' newline];
            else
                dispersion_text = '';
            end
            
            % Geometrical counterpoise correction (gCP) options
            if strcmpi(Counterpoise,'gCP')
                gCP_text = ['GCP' newline 'METHOD' newline 'dft/pobtz' newline 'PRINTEMISS' newline 'END' newline];
                gCP_label = '_gCP';
            else
                gCP_text = '';
                gCP_label = '';
            end
            
            Theory_dir = [label '_' regexprep(Theory,'\((.?)\)','_$1') disptext gCP_label];
            dir1 = [workDir filesep Theory_dir];
            % create/change dir
            if ~exist(dir1,'dir')
                mkdir(dir1);
            end
            cd(dir1)

            Coordinates_Copy = Coordinates; % Save Coordinates variable
            for j = 1:length(Metal_Basis_set)
                Current_M_basis = [Metal_Basis_set{j} Aug_Label_M];
                % Input coordinates
                [match,~] = regexp(Current_Salt,...
                    {'Li' 'Na' 'K' 'Rb' 'F' 'Cl' 'Br' 'I'},'match','tokens');
                matches = match(~cellfun('isempty',match));
                
                for jj = 1:length(Halide_Basis_set)
                    Coordinates = Coordinates_Copy; % Reload fresh coordinates variable
            
                    Current_H_basis = [Halide_Basis_set{jj} Aug_Label_H];
                    
                    if strcmp(Current_M_basis,Current_H_basis) % If both basis sets are the same, combine them.
                        Current_basis = Current_M_basis;
                        Aug_Label = Aug_Label_M;
                        OneBasis = true;
                    elseif strcmp(Structure,'Ion') || strcmp(Structure,'Atom') % If only one basis set type used
                        OneBasis = true;
                        if contained_in_cell(matches{1}{1},{'Li' 'Na' 'K' 'Rb'})
                            Current_basis = Current_M_basis; 
                            Aug_Label = Aug_Label_M;
                        elseif contained_in_cell(matches{1}{1},{'F' 'Cl' 'Br' 'I'})
                            Current_basis = Current_H_basis;
                            Aug_Label = Aug_Label_H;
                        else
                            error(['Unknown Match: ' matches])
                        end
                    else % If using a mixed basis set
                        Current_basis = ['M' Current_M_basis '_H' Current_H_basis];
                        OneBasis = false;
                    end

                    for k=1:length(matches)
                        if strcmp('Li',matches{k}{1})
                            Coordinates = strrep(Coordinates,'##Metal##','3');
                        elseif strcmp('Na',matches{k}{1})
                            Coordinates = strrep(Coordinates,'##Metal##','11');
                        elseif strcmp('K',matches{k}{1})
                            Coordinates = strrep(Coordinates,'##Metal##','19');
                        elseif strcmp('Rb',matches{k}{1})
                            Coordinates = strrep(Coordinates,'##Metal##','237');
                        elseif strcmp('F',matches{k}{1})
                            Coordinates = strrep(Coordinates,'##Halide##','9');
                        elseif strcmp('Cl',matches{k}{1})
                            Coordinates = strrep(Coordinates,'##Halide##','17');
                        elseif strcmp('Br',matches{k}{1})
                            if strcmpi(Current_H_basis,'7-311G')
                                Coordinates = strrep(Coordinates,'##Halide##','235');
                            else
                                Coordinates = strrep(Coordinates,'##Halide##','35');
                            end
                        elseif strcmp('I',matches{k}{1})
                            if strcmpi(Current_H_basis,'STO-3G')
                                Coordinates = strrep(Coordinates,'##Halide##','53');
                            else
                                Coordinates = strrep(Coordinates,'##Halide##','253');
                            end
                        else
                            error(['Error: Unknown Element Type: ' matches{k}{1}]);
                        end
                    end
                    Coordinates = strrep(Coordinates,'##Metal##','');
                    Coordinates = strrep(Coordinates,'##Halide##','');

                    run_min = false;
                    if findBSSEparameters
                        [a,b,c,CoordsArray] = FindCrystalParams(Current_Salt,...
                            Structure,Theory,Current_basis,home,...
                            BSSE_id);
                        params = [' ' regexprep(num2str(a.*100),' +',',')];
                        display_text = [' ' num2str(min(a),'%5.4f') ' - ' num2str(max(a),'%5.4f') OutSymTxt];
                    elseif Calc_Vib || QHA
                        if isempty(Relax_Basis_Set)
                            Reference_Basis = Current_basis;
                        else
                            Reference_Basis = Relax_Basis_Set;
                        end
                        if isempty(Relax_XC)
                            Reference_Theory = Theory;
                        else
                            Reference_Theory = Relax_XC;
                        end
                        if strcmpi(Relax_Dispersion,'Off')
                            Reference_Dispersion = disptext;
                        elseif isempty(Relax_Dispersion)
                            Reference_Dispersion = '';
                        else
                            Reference_Dispersion = ['_' Relax_Dispersion];
                        end
                        if strcmpi(Relax_gCP,'Off')
                            Reference_gCP = gCP_label;
                        elseif isempty(Relax_gCP)
                            Reference_gCP = '';
                        else
                            Reference_gCP = ['_' Relax_gCP];
                        end

                        [a_out,b_out,c_out,Coordinates_out,run_min] = FindMinEParams(...
                            Current_Salt,Structure,Reference_Theory,Theory,...
                            Reference_Basis,Current_basis,home,Reference_Dispersion,disptext,...
                            Reference_gCP,gCP_label,false,AutoFindDataTypes,false);

                        % No data found, run minimization prior to freqency calc
                        if run_min
                            Opt_Input = Cell_Optimization_Text('FULLOPTG',...
                                opt_conv,max_cycle_opt,NUMGRCEL,TOLDEG,...
                                TOLDEX,EXTPRESS,Internal_Coordinates,Updating_Scheme);
                            Vib_Input_temp = Vibrational_Analysis_Text(QHA,...
                                QHA_Step,QHA_Points,QHA_RANGE,T,P,...
                                Phonon_Dispersion,Supercell,PDOS,PDOS_Range,...
                                PDOS_projected,Opt_Input);

                            % Data found
                        else                       
                            a = a_out;
                            b = b_out;
                            c = c_out;
                            Coordinates = Coordinates_out;
                        end
                    elseif Cell_optimization
                        if AutoFindCellParams
                            if isempty(Relax_Basis_Set)
                                Reference_Basis = Current_basis;
                            else
                                Reference_Basis = Relax_Basis_Set;
                            end
                            if isempty(Relax_XC)
                                Reference_Theory = Theory;
                            else
                                Reference_Theory = Relax_XC;
                            end
                            if strcmpi(Relax_Dispersion,'Off')
                                Reference_Dispersion = disptext;
                            elseif isempty(Relax_Dispersion)
                                Reference_Dispersion = '';
                            else
                                Reference_Dispersion = ['_' Relax_Dispersion];
                            end
                            if strcmpi(Relax_gCP,'Off')
                                Reference_gCP = gCP_label;
                            elseif isempty(Relax_gCP)
                                Reference_gCP = '';
                            else
                                Reference_gCP = ['_' Relax_gCP];
                            end
                            
                            [a_out,b_out,c_out,Coordinates_out,run_min] = FindMinEParams(...
                                Current_Salt,Structure,Reference_Theory,Theory,...
                                Reference_Basis,Current_basis,home,Reference_Dispersion,disptext,...
                                Reference_gCP,gCP_label,false,AutoFindDataTypes,P1_Symmetry);
                            
                            if ~run_min                       
                                a = a_out;
                                b = b_out;
                                c = c_out;
                                Coordinates = Coordinates_out;
                            end
                        else
                            run_min = true;
                        end
                    elseif AutoFindCellParams
                        if isempty(Relax_Basis_Set)
                            Reference_Basis = Current_basis;
                        else
                            Reference_Basis = Relax_Basis_Set;
                        end
                        if isempty(Relax_XC)
                            Reference_Theory = Theory;
                        else
                            Reference_Theory = Relax_XC;
                        end
                        if strcmpi(Relax_Dispersion,'Off')
                            Reference_Dispersion = disptext;
                        elseif isempty(Relax_Dispersion)
                            Reference_Dispersion = '';
                        else
                            Reference_Dispersion = ['_' Relax_Dispersion];
                        end
                        if strcmpi(Relax_gCP,'Off')
                            Reference_gCP = gCP_label;
                        elseif isempty(Relax_gCP)
                            Reference_gCP = '';
                        else
                            Reference_gCP = ['_' Relax_gCP];
                        end
                        
                        [a_out,b_out,c_out,Coordinates_out,run_min] = FindMinEParams(...
                            Current_Salt,Structure,Reference_Theory,Theory,...
                            Reference_Basis,Current_basis,home,Reference_Dispersion,disptext,...
                            Reference_gCP,gCP_label,false,AutoFindDataTypes,P1_Symmetry);
                        if ~run_min % Run min is true if no unit cell parameters found
                            a = a_out;
                            b = b_out;
                            c = c_out;
                            Coordinates = Coordinates_out;
                        end
                    end
                    
                    if ~Calc_Vib && ~QHA && ~Cell_optimization && ...
                            ~findBSSEparameters && ~strcmp('Ion',Structure) && ...
                            ~strcmp('Atom',Structure)
                        if length(a) > 1
                            display_text = [' ' num2str(min(a),'%5.2f') ' - ' num2str(max(a),'%5.2f') OutSymTxt];
                        else
                            display_text = [' ' num2str(a,'%5.4f') OutSymTxt];
                        end
                    end
                    
                    if ~isempty(BSSE_id) && (strcmpi(Structure,'BetaBeO') || N_Cluster > 0)
                        % Save unit cell data into mat file
                        FC = textscan(Coordinates,'%*f %f %f %f','HeaderLines',1);
                        FC_Metal = {Metal FC{1}(1) FC{2}(1) FC{3}(1)};
                        FC_Halide = {Halide FC{1}(2) FC{2}(2) FC{3}(2)};
                        %Coordinates = Gen_Cluster_BetaBeO(Current_Salt,Coordinates,BSSE_id,...
                        %    BSSE_Nstars,BSSE_Rmax,a,b,c);
                        Coordinates = Gen_Cluster_Salt(Current_Salt,Structure,Coordinates,BSSE_id,...
                            BSSE_Nstars,BSSE_Rmax,a,b,c,N_Cluster);
                        
                        
                        
                    end
                
                    Basis_dir = ['BASS_' Current_basis];
                    dir2 = [workDir filesep Theory_dir filesep Basis_dir];
                    % create/change dir
                    if ~exist(dir2,'dir')
                        mkdir(dir2);
                    end
                    cd(dir2)

                    % Open template for input file
                    fid = fopen([home filesep 'templates' filesep template_file],'rt');
                    X = fread(fid);
                    fclose(fid);
                    input_text = char(X.');

                    %% Replace strings
                    % Vibration switch 
                    if Calc_Vib || QHA
                        if run_min
                            input_text = strrep(input_text,['##VIBRATION##' newline],Vib_Input_temp);
                        else
                            input_text = strrep(input_text,['##VIBRATION##' newline],Vib_Input);
                        end
                    else
                        input_text = strrep(input_text,['##VIBRATION##' newline],'');
                    end
                
                    % Geometry
                    input_text = strrep(input_text,'##PERIODICITY##',Periodicity);
                    input_text = strrep(input_text,'##SPACEGROUP##',Space_group);
                    input_text = strrep(input_text,['##CELLPARAM##' newline],Cell_param);
                    if ~findBSSEparameters
                        input_text = strrep(input_text,'##GEOM##',Coordinates);
                    end
                    input_text = strrep(input_text,['##LATVEC##' newline],Latvectext);
                    input_text = strrep(input_text,['##BSSE##' newline],BSSE_text);

                    % DFT Settings
                    input_text = strrep(input_text,['##DFT##' newline],DFT);
                    input_text = strrep(input_text,['##DFTD3##' newline],[dispersion_text gCP_text]);

                    % Hamiltonian and SCF parameters
                    input_text = strrep(input_text,['##KPOINTS##' newline],kpoints);
                    input_text = strrep(input_text,'##INTTOL##',integration_params);
                    input_text = strrep(input_text,['##SAVEWF##' newline],SAVEWF);
                    input_text = strrep(input_text,'##MAXCYCL##',max_cycle);
                    input_text = strrep(input_text,'##CONV##',convergence);
                    input_text = strrep(input_text,'##ILAS##',ILA);
                    input_text = strrep(input_text,'##DIIS##',DIIS_Text);
                    input_text = strrep(input_text,['##ANDERSON##' newline],Anderson_text);
                    input_text = strrep(input_text,['##BROYDEN##' newline],Broyden_text);
                    input_text = strrep(input_text,['##FMIXING##' newline],FMIXING);
                    input_text = strrep(input_text,['##LEVSHIFT##' newline],LEVSHIFT);
                    input_text = strrep(input_text,'##BIPOSIZ##',biposize);
                    input_text = strrep(input_text,'##EXCHSIZ##',exchsize);

                    % load the basis sets
                    if OneBasis % If both metal and halide are the same basis set
                        switch lower(Current_basis)
                            case lower(['pob-TZVP' Aug_Label])
                                [Li,Na,K,Rb,F,Cl,Br,I] = Prune_Scale_pob_TZVP(home,Prune_M1,...
                                    Prune_M2,Prune_H1,Prune_H2,Prune_H3,Scale_M1,...
                                    Scale_M2,Scale_H1,Scale_H2,Scale_H3,...
                                    Add_M,Add_H,NumM,NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['pob-TZVPP' Aug_Label])
                                [Li,Na,K,F,Cl] = Prune_Scale_pob_TZVPP(home,Prune_M1,...
                                    Prune_M2,Prune_H1,Prune_H2,Prune_H3,Scale_M1,...
                                    Scale_M2,Scale_H1,Scale_H2,Scale_H3,...
                                    Add_M,Add_H,NumM,NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['pob-TZVP-rev2' Aug_Label])
                                [Li,Na,K,F,Cl,Br] = Prune_Scale_pob_TZVP_rev2(home,Prune_M1,...
                                    Prune_M2,Prune_H1,Prune_H2,Prune_H3,Scale_M1,...
                                    Scale_M2,Scale_H1,Scale_H2,Scale_H3,...
                                    Add_M,Add_H,NumM,NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-AVTZ' Aug_Label])
                                I = Prune_Scale_ECP28MDF_AVTZ(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-VTZ' Aug_Label])
                                I = Prune_Scale_ECP28MDF_VTZ(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-VQZ' Aug_Label])
                                I = Prune_Scale_ECP28MDF_VQZ(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-AVQZ' Aug_Label])
                                I = Prune_Scale_ECP28MDF_AVQZ(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-V5Z' Aug_Label])
                                I = Prune_Scale_ECP28MDF_V5Z(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-AV5Z' Aug_Label])
                                I = Prune_Scale_ECP28MDF_AV5Z(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['def2-TZVPD' Aug_Label])
                                [Li,Na,K,Rb,F,Cl,Br,I] = Prune_Scale_def2_TZVPD(home,Structure,BSSE_id,...
                                Initialize_As_Ions);
                            case lower('7-311G')
                                [Li,Na,K,Rb,F,Cl,Br,I] = Prune_Scale_7_311G(home,Structure,...
                                Initialize_As_Ions);
                            case lower('STO-3G')
                                [Li,Na,K,Rb,F,Cl,Br,I] = Prune_Scale_STO_3G(home,Structure,...
                                Initialize_As_Ions);
                            case lower({'def2-SVPD' 'def2-TZVP' 'def2-QZVP' 'def2-QZVPD' 'cc-pV5Z' 'aug-cc-pV5Z' 'aug-cc-pV5ZD'})
                                [Li,Na,F,Cl,Br,I] = Prune_Scale_Basis(strrep(Current_basis,'-','_'),home,Structure,BSSE_id,...
                                Initialize_As_Ions);
                            otherwise
                                error(['Unknown Basis Set: ' Current_basis])
                        end
                    else % Two basis sets used
                        % First Load the metal basis set
                        switch lower(Current_M_basis)
                            case lower(['pob-TZVP' Aug_Label_M])
                                [Li,Na,K,Rb,~,~,~,~] = Prune_Scale_pob_TZVP(home,Prune_M1,...
                                    Prune_M2,Prune_H1,Prune_H2,Prune_H3,Scale_M1,...
                                    Scale_M2,Scale_H1,Scale_H2,Scale_H3,...
                                    Add_M,Add_H,NumM,NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['pob-TZVPP' Aug_Label_M])
                                [Li,Na,K,~,~] = Prune_Scale_pob_TZVPP(home,Prune_M1,...
                                    Prune_M2,Prune_H1,Prune_H2,Prune_H3,Scale_M1,...
                                    Scale_M2,Scale_H1,Scale_H2,Scale_H3,...
                                    Add_M,Add_H,NumM,NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['pob-TZVP-rev2' Aug_Label_M])
                                [Li,Na,K,~,~,~] = Prune_Scale_pob_TZVP_rev2(home,Prune_M1,...
                                    Prune_M2,Prune_H1,Prune_H2,Prune_H3,Scale_M1,...
                                    Scale_M2,Scale_H1,Scale_H2,Scale_H3,...
                                    Add_M,Add_H,NumM,NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['def2-TZVPD' Aug_Label_M])
                                [Li,Na,K,Rb,~,~,~,~] = Prune_Scale_def2_TZVPD(home,Structure,BSSE_id,...
                                Initialize_As_Ions);
                            case lower('7-311G')
                                [Li,Na,K,Rb,~,~,~,~] = Prune_Scale_7_311G(home,Structure,...
                                Initialize_As_Ions);
                            case lower('STO-3G')
                                [Li,Na,K,Rb,~,~,~,~] = Prune_Scale_STO_3G(home,Structure,...
                                Initialize_As_Ions);
                            case lower({'def2-SVPD' 'def2-TZVP' 'def2-QZVP' 'def2-QZVPD' 'cc-pV5Z' 'aug-cc-pV5Z' 'aug-cc-pV5ZD'})
                                [Li,Na,~,~,~,~] = Prune_Scale_Basis(strrep(Current_M_basis,'-','_'),home,Structure,BSSE_id,...
                                Initialize_As_Ions);
                            otherwise
                            error(['Unknown Basis Metal Set: ' Current_M_basis])
                        end
                    
                        % Load the halide basis set
                        switch lower(Current_H_basis)
                            case lower(['pob-TZVP' Aug_Label_H])
                                [~,~,~,~,F,Cl,Br,I] = Prune_Scale_pob_TZVP(home,Prune_M1,...
                                    Prune_M2,Prune_H1,Prune_H2,Prune_H3,Scale_M1,...
                                    Scale_M2,Scale_H1,Scale_H2,Scale_H3,...
                                    Add_M,Add_H,NumM,NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['pob-TZVPP' Aug_Label_H])
                                [~,~,~,F,Cl] = Prune_Scale_pob_TZVPP(home,Prune_M1,...
                                    Prune_M2,Prune_H1,Prune_H2,Prune_H3,Scale_M1,...
                                    Scale_M2,Scale_H1,Scale_H2,Scale_H3,...
                                    Add_M,Add_H,NumM,NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['pob-TZVP-rev2' Aug_Label_H])
                                [~,~,~,F,Cl,Br] = Prune_Scale_pob_TZVP_rev2(home,Prune_M1,...
                                    Prune_M2,Prune_H1,Prune_H2,Prune_H3,Scale_M1,...
                                    Scale_M2,Scale_H1,Scale_H2,Scale_H3,...
                                    Add_M,Add_H,NumM,NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-AVTZ' Aug_Label_H])
                                I = Prune_Scale_ECP28MDF_AVTZ(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-VTZ' Aug_Label_H])
                                I = Prune_Scale_ECP28MDF_VTZ(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-VQZ' Aug_Label_H])
                                I = Prune_Scale_ECP28MDF_VQZ(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-AVQZ' Aug_Label_H])
                                I = Prune_Scale_ECP28MDF_AVQZ(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-V5Z' Aug_Label_H])
                                I = Prune_Scale_ECP28MDF_V5Z(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['ECP28MDF-AV5Z' Aug_Label_H])
                                I = Prune_Scale_ECP28MDF_AV5Z(home,Add_H,...
                                    NumH,Structure,BSSE_id,...
                                    Initialize_As_Ions);
                            case lower(['def2-TZVPD' Aug_Label_H])
                                [~,~,~,~,F,Cl,Br,I] = Prune_Scale_def2_TZVPD(home,Structure,BSSE_id,...
                                Initialize_As_Ions);
                            case lower('7-311G')
                                [~,~,~,~,F,Cl,Br,I] = Prune_Scale_7_311G(home,Structure,...
                                Initialize_As_Ions);
                            case lower('STO-3G')
                                [~,~,~,~,F,Cl,Br,I] = Prune_Scale_STO_3G(home,Structure,...
                                Initialize_As_Ions);
                            case lower({'def2-SVPD' 'def2-TZVP' 'def2-QZVP' 'def2-QZVPD' 'cc-pV5Z' 'aug-cc-pV5Z' 'aug-cc-pV5ZD'})
                                [~,~,F,Cl,Br,I] = Prune_Scale_Basis(strrep(Current_H_basis,'-','_'),home,Structure,BSSE_id,...
                                Initialize_As_Ions);
                            otherwise
                                error(['Unknown Basis Set: ' Current_basis])
                        end
                    end

                    % Preprocessing basis set
                    basis_set = '';
                    for k=1:length(matches)
                        if strcmp('Li',matches{k}{1}) %#ok<*AGROW>
                            basis_set = [basis_set newline Li];
                            Met_Bas = Li;
                        elseif strcmp('Na',matches{k}{1})
                            basis_set = [basis_set newline Na];
                            Met_Bas = Na;
                        elseif strcmp('K',matches{k}{1})
                            basis_set = [basis_set newline K];
                            Met_Bas = K;
                        elseif strcmp('Rb',matches{k}{1})
                            basis_set = [basis_set newline Rb];
                            Met_Bas = Rb;
                        elseif strcmp('F',matches{k}{1})
                            basis_set = [basis_set newline F];
                            Hal_Bas = F;
                        elseif strcmp('Cl',matches{k}{1})
                            basis_set = [basis_set newline Cl];
                            Hal_Bas = Cl;
                        elseif strcmp('Br',matches{k}{1})
                            basis_set = [basis_set newline Br];
                            Hal_Bas = Br;
                        elseif strcmp('I',matches{k}{1})
                            basis_set = [basis_set newline I];
                            Hal_Bas = I;
                        else
                            error(['Error obtaining basis set: ' Current_basis ' for ' Current_Salt])
                        end
                    end
                    basis_set = regexprep(basis_set,'^\n','');
                    clearvars Li Na K Rb F Cl Br I

                    % Prep basis set for BSSE with BetaBeO
                    if ~isempty(BSSE_id) && strcmpi(Structure,'BetaBeO')
                        basis_set = Ghost_Basis_BetaBeO(Current_Salt,Met_Bas,Hal_Bas,BSSE_id);
                    end
                    
                    % Basis Set input
                    input_text = strrep(input_text,'##BASIS##',basis_set);
                    
                    %% For Single Ion/Atom or pair calculations
                    if strcmp(Periodicity,'MOLECULE') && ~Cell_optimization && ~strcmpi(Structure,'BetaBeO')

                        % Add bond length if salt pair
                        if strcmp(Structure,'Pair')
                            for k=1:length(a)

                                % Copy Input text to new variable so as not to overwrite it
                                input_text_a = input_text;
                                a_param = num2str(a(k),'%15.14f');

                                % Generate Title Text
                                Title = [Current_Salt ' ' Structure ' ' Theory disptext gCP_label ...
                                    ' with ' Current_basis ' Basis Set and Bond Length of ' num2str(a(k)) ' Angstroms.'];

                                param_dir = ['APAR_' num2str(a(k),'%05.2f') SymmetryMod];

                                % create/change dir
                                dir3 = [workDir filesep Theory_dir filesep Basis_dir filesep param_dir];
                                if ~exist(dir3,'dir')
                                    mkdir(dir3);
                                end
                                cd(dir3)

                                % Replace strings Unique to Unit Cell with Specific Length
                                if findBSSEparameters
                                    input_text_a = strrep(input_text_a,'##GEOM##',Coordinates);
                                end
                                input_text_a = strrep(input_text_a,'##TITLE##',Title);
                                input_text_a = strrep(input_text_a,'##APAR##',a_param);
                                input_text_a = strrep(input_text_a,[newline '##CELLOPT##'],'');

                                if Skip_Running && strcmpi(server(1:3),'sea')
                                    JobRunning = CheckRunning(dir3,Stop_Running);
                                    if JobRunning
                                        continue
                                    end
                                end
                                
                                % Save file
                                fid2 = fopen( [Theory_dir '_' Current_basis '-' num2str(a(k),'%05.2f') '.d12'],'wt');
                                fwrite(fid2,input_text_a);
                                fclose(fid2);

                                cd(dir2) % Return to outer directory
                            end
                        else
                            % Generate Title Text
                            Title = [Current_Salt ' ' Structure ' ' Theory disptext gCP_label ...
                                ' with ' Current_basis ' Basis Set.'];

                            % Replace Title String
                            input_text = strrep(input_text,'##TITLE##',Title);
                            input_text = strrep(input_text,[newline '##CELLOPT##'],'');

                            if Skip_Running && strcmpi(server(1:3),'sea')
                                JobRunning = CheckRunning(dir2,Stop_Running);
                                if JobRunning
                                    continue
                                end
                            end
                            
                            % Save file
                            fid2 = fopen( [Theory_dir '_' Current_basis '.d12'],'wt');
                            fwrite(fid2,input_text);
                            fclose(fid2);

                            a = nan; % to prevent error
                        end
                    elseif Calc_Vib || QHA
                        if Skip_Completed
                            SkipSwitch = FindCompletedJobs(Current_Salt,Structure,Theory,Current_basis,...
                                Dispersion,Counterpoise,P1_Symmetry,External_P,Cell_optimization,...
                                opt_type,opt_conv,QHA,Calc_Vib,Supercell,a,home);
                            if SkipSwitch
                                continue
                            end
                        end
                        
                        % Running minimization
                        if run_min
                            % Generate Formatted a (possibly b and c) cell parameters
                            if length(a) == 1
                                a_param = num2str(a,'%15.14f');
                            else
                                error('Please ensure the length of initial cell inputs is unity for optimization runs.');
                            end

                            if strcmp(Structure,'Wurtzite')
                                c_param = num2str(a*c_over_a_Wurtzite,'%15.14f');
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            elseif strcmp(Structure,'NiAs')
                                c_param = num2str(a*c_over_a_NiAs,'%15.14f');
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            elseif strcmp(Structure,'BetaBeO')
                                c_param = num2str(a*c_over_a_BetaBeO,'%15.14f');
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            elseif strcmp(Structure,'FiveFive')
                                b_param = num2str(a*b_over_a_FiveFive,'%15.14f');
                                c_param = num2str(a*c_over_a_FiveFive,'%15.14f');
                                input_text = strrep(input_text,'##BPAR##',b_param);
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            end
                        else
                            a_param = num2str(a,'%15.14f');
                            b_param = num2str(b,'%15.14f');
                            c_param = num2str(c,'%15.14f');

                            input_text = strrep(input_text,'##BPAR##',b_param);
                            input_text = strrep(input_text,'##CPAR##',c_param);
                        end

                        input_text = strrep(input_text,'##APAR##',a_param);

                        % Generate Title Text
                        if QHA
                            Title = ['Quasi-Harmonic Approx. Vib. Analysis of ' Current_Salt ' ' Structure ' ' Theory disptext gCP_label ...
                            ' with ' Current_basis ' Basis Set and Cell Vector |a| = ' a_param];
                            if Supercell > 1
                                param_dir = ['QHA_SS' num2str(Supercell) SymmetryMod];
                            else
                                param_dir = ['QHA' SymmetryMod];
                            end
                        else
                            Title = ['Vibrational Analysis of ' Current_Salt ' ' Structure ' ' Theory disptext gCP_label ...
                            ' with ' Current_basis ' Basis Set and Cell Vector |a| = ' a_param];

                            if Phonon_Dispersion && Supercell > 1
                                param_dir = ['VIBANALYSIS_SS' num2str(Supercell) SymmetryMod];
                            else
                                param_dir = ['VIBANALYSIS' SymmetryMod];
                            end
                        end

                        % create/change dir
                        dir3 = [workDir filesep Theory_dir filesep Basis_dir filesep param_dir];
                        if ~exist(dir3,'dir')
                            mkdir(dir3);
                        end
                        
                        cd(dir3)

                        % Replace strings Unique to Unit Cell with Specific Length
                        input_text = strrep(input_text,'##TITLE##',Title);
                        input_text = strrep(input_text,['##CELLOPT##' newline],'');

                        if Skip_Running && strcmpi(server(1:3),'sea')
                            JobRunning = CheckRunning(dir3,Stop_Running);
                            if JobRunning
                                continue
                            end
                        end
                        
                        % Save file
                        if QHA
                            fid2 = fopen( [Theory_dir '_' Current_basis '_QHA.d12'],'wt');
                        else
                            fid2 = fopen( [Theory_dir '_' Current_basis '_Vib.d12'],'wt');
                        end
                        fwrite(fid2,input_text);
                        fclose(fid2);

                    %% For Unit Cell optimization
                    elseif Cell_optimization
                        if Skip_Completed
                            SkipSwitch = FindCompletedJobs(Current_Salt,Structure,Theory,Current_basis,...
                                Dispersion,Counterpoise,P1_Symmetry,External_P,Cell_optimization,...
                                opt_type,opt_conv,QHA,Calc_Vib,Supercell,a,home);
                            if SkipSwitch
                                continue
                            end
                        end
                        % Generate cell optimization text
                        CellOpt_text = Cell_Optimization_Text(opt_type,opt_conv,...
                            max_cycle_opt,NUMGRCEL,TOLDEG,TOLDEX,EXTPRESS,...
                            Internal_Coordinates,Updating_Scheme);

                        % Generate Formatted a (possibly b and c) cell parameters
                        if length(a) == 1
                            a_param = num2str(a,'%15.14f');
                        else
                            error('Please ensure the length of initial cell inputs is unity for optimization runs.');
                        end

                        if run_min
                            if strcmp(Structure,'Wurtzite')
                                c_param = num2str(a*c_over_a_Wurtzite,'%15.14f');
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            elseif strcmp(Structure,'NiAs')
                                c_param = num2str(a*c_over_a_NiAs,'%15.14f');
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            elseif strcmp(Structure,'BetaBeO')
                                c_param = num2str(a*c_over_a_BetaBeO,'%15.14f');
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            elseif strcmp(Structure,'FiveFive')
                                b_param = num2str(a*b_over_a_FiveFive,'%15.14f');
                                c_param = num2str(a*c_over_a_FiveFive,'%15.14f');

                                input_text = strrep(input_text,'##BPAR##',b_param);
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            end
                        else
                            a_param = num2str(a,'%15.14f');
                            b_param = num2str(b,'%15.14f');
                            c_param = num2str(c,'%15.14f');

                            input_text = strrep(input_text,'##BPAR##',b_param);
                            input_text = strrep(input_text,'##CPAR##',c_param);
                        end

                        % Generate Title Text
                        Title = ['Unit Cell Optimization of ' Current_Salt ' ' Structure ' ' Theory disptext gCP_label ...
                            ' with ' Current_basis ' Basis Set and Initial Cell Vector |a| = ' num2str(a)];

                        if strcmp(opt_type,'CELLONLY')
                            param_dir = ['CELLOPT' DirPTxt SymmetryMod];
                        elseif strcmp(opt_type,'FULLOPTG') || strcmp(opt_type,'ITATOCEL')
                            param_dir = ['FULLOPT' DirPTxt SymmetryMod];
                        elseif strcmp(opt_type,'ATOMONLY')
                            param_dir = ['ATOMOPT' DirPTxt SymmetryMod];
                        else
                            error(['Unknown Optimization Type: ' opt_type])
                        end               

                        % create/change dir
                        [~,~,~] = mkdir(param_dir);
                        dir3 = [workDir filesep Theory_dir filesep Basis_dir filesep param_dir];
                        cd(dir3)

                        % Replace strings Unique to Unit Cell with Specific Length
                        input_text = strrep(input_text,'##TITLE##',Title);
                        input_text = strrep(input_text,'##APAR##',a_param);
                        input_text = strrep(input_text,'##CELLOPT##',CellOpt_text);

                        if Skip_Running && strcmpi(server(1:3),'sea')
                            JobRunning = CheckRunning(dir3,Stop_Running);
                            if JobRunning
                                continue
                            end
                        end
                        
                        % Save file
                        fid2 = fopen( [Theory_dir '_' Current_basis '_CellOpt.d12'],'wt');
                        fwrite(fid2,input_text);
                        fclose(fid2);
                        
                    elseif AutoFindCellParams
                        if Skip_Completed
                            SkipSwitch = FindCompletedJobs(Current_Salt,Structure,Theory,Current_basis,...
                                Dispersion,Counterpoise,P1_Symmetry,External_P,Cell_optimization,...
                                opt_type,opt_conv,QHA,Calc_Vib,Supercell,a,home);
                            if SkipSwitch
                                continue
                            end
                        end
                        % Generate Formatted a (possibly b and c) cell parameters
                        if length(a) > 1
                            error('Please ensure the length of initial cell inputs is unity for auto-find-parameter runs.');
                        end

                        if run_min % If no parameters were found, use defaults
                            disp('No parameters found in search, using defaults.')
                            a_param = num2str(a,'%15.14f');
                            if strcmp(Structure,'Wurtzite')
                                c_param = num2str(a*c_over_a_Wurtzite,'%15.14f');
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            elseif strcmp(Structure,'NiAs')
                                c_param = num2str(a*c_over_a_NiAs,'%15.14f');
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            elseif strcmp(Structure,'BetaBeO')
                                c_param = num2str(a*c_over_a_BetaBeO,'%15.14f');
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            elseif strcmp(Structure,'FiveFive')
                                b_param = num2str(a*b_over_a_FiveFive,'%15.14f');
                                c_param = num2str(a*c_over_a_FiveFive,'%15.14f');

                                input_text = strrep(input_text,'##BPAR##',b_param);
                                input_text = strrep(input_text,'##CPAR##',c_param);
                            end
                        else
                            a_param = num2str(a,'%15.14f');
                            b_param = num2str(b,'%15.14f');
                            c_param = num2str(c,'%15.14f');

                            input_text = strrep(input_text,'##BPAR##',b_param);
                            input_text = strrep(input_text,'##CPAR##',c_param);
                        end
                        optTag = ['-' num2str(a,'%05.2f')];

                        % Generate Title Text
                        Title = [Current_Salt ' ' Structure ' ' Theory disptext gCP_label ...
                            ' with ' Current_basis ' Basis Set and Cell Vector |a| = ' num2str(a)];

                        param_dir = ['APAR_' num2str(a,'%05.2f') SymmetryMod];

                        % create/change dir
                        dir3 = [workDir filesep Theory_dir filesep Basis_dir filesep param_dir];
                        if ~exist(dir3,'dir')
                            mkdir(dir3);
                        end
                        cd(dir3)
                        
                        % If BetaBeO BSSE, save structure data to mat file
                        if ~isempty(BSSE_id) && strcmpi(Structure,'BetaBeO')
                            save([dir3 filesep Current_Salt labend '_' Theory '_' ...
                                Current_basis '-' num2str(a,'%05.2f') '.mat'],...
                                'a','b','c','FC_Metal','FC_Halide')
                        end
                        FC_Metal = {Metal FC{1}(1) FC{2}(1) FC{3}(1)};
                        FC_Halide = {Halide FC{1}(2) FC{2}(2) FC{3}(2)};
                        

                        % Replace strings Unique to Unit Cell with Specific Length
                        input_text = strrep(input_text,'##TITLE##',Title);
                        input_text = strrep(input_text,'##APAR##',a_param);
                        input_text = strrep(input_text,[newline '##CELLOPT##'],'');

                        if Skip_Running && strcmpi(server(1:3),'sea')
                            JobRunning = CheckRunning(dir3,Stop_Running);
                            if JobRunning
                                continue
                            end
                        end
                        
                        % Save file
                        fid2 = fopen( [Theory_dir '_' Current_basis '-' num2str(a,'%05.2f') '.d12'],'wt');
                        fwrite(fid2,input_text);
                        fclose(fid2);
                        
                    %% For Periodic Systems
                    else
                        % Skip jobs for which data already exists
                        if Skip_Completed
                            SkipSwitch = FindCompletedJobs(Current_Salt,Structure,Theory,Current_basis,...
                                Dispersion,Counterpoise,P1_Symmetry,External_P,Cell_optimization,...
                                opt_type,opt_conv,QHA,Calc_Vib,Supercell,a,home);
                            a = a(SkipSwitch);
                            
                            % Update display text
                            if length(a) > 1
                                display_text = [' ' num2str(min(a),'%2.2f') ' - ' num2str(max(a),'%2.2f')];
                            else
                                display_text = [' ' num2str(a,'%2.2f')];
                            end
                        end
                        
                        for k=1:length(a)

                            % Copy Input text to new variable so as not to overwrite it
                            input_text_a = input_text;
                            a_param = num2str(a(k),'%15.14f');

                            if findBSSEparameters
                                CoordsCell = CoordsArray{k};
                                fcnwrapper = @(x) num2str(x,'%6.4f');
                                CoordsCell(:,2:4) = cellfun(fcnwrapper,CoordsCell(:,2:4),'UniformOutput',0);

                                Coordinates = ['2' newline ...
                                    atomlabel2num(CoordsCell{1,1}) ' ' CoordsCell{1,2} ' ' CoordsCell{1,3} ' ' CoordsCell{1,3} newline ...
                                    atomlabel2num(CoordsCell{2,1}) ' ' CoordsCell{2,2} ' ' CoordsCell{2,3} ' ' CoordsCell{2,3}];

                                input_text_a = strrep(input_text_a,'##GEOM##',Coordinates);

                                b_param = num2str(b(k),'%15.14f');
                                c_param = num2str(c(k),'%15.14f');
                                if strcmp(Structure,'Wurtzite')
                                    input_text_a = strrep(input_text_a,'##CPAR##',c_param);
                                elseif strcmp(Structure,'NiAs')
                                    input_text_a = strrep(input_text_a,'##CPAR##',c_param);
                                elseif strcmp(Structure,'BetaBeO')
                                    input_text_a = strrep(input_text_a,'##CPAR##',c_param);
                                elseif strcmp(Structure,'FiveFive')
                                    input_text_a = strrep(input_text_a,'##BPAR##',b_param);
                                    input_text_a = strrep(input_text_a,'##CPAR##',c_param);
                                end
                            else
                                % Generate Formatted cell parameters
                                if strcmp(Structure,'Wurtzite')
                                    c_param = num2str(a(k)*c_over_a_Wurtzite,'%15.14f');
                                    input_text_a = strrep(input_text_a,'##CPAR##',c_param);
                                elseif strcmp(Structure,'NiAs')
                                    c_param = num2str(a(k)*c_over_a_NiAs,'%15.14f');
                                    input_text_a = strrep(input_text_a,'##CPAR##',c_param);
                                elseif strcmp(Structure,'BetaBeO')
                                    c_param = num2str(a(k)*c_over_a_BetaBeO,'%15.14f');
                                    input_text_a = strrep(input_text_a,'##CPAR##',c_param);
                                elseif strcmp(Structure,'FiveFive')
                                    b_param = num2str(a(k)*b_over_a_FiveFive,'%15.14f');
                                    c_param = num2str(a(k)*c_over_a_FiveFive,'%15.14f');
                                    input_text_a = strrep(input_text_a,'##BPAR##',b_param);
                                    input_text_a = strrep(input_text_a,'##CPAR##',c_param);
                                end
                            end

                            % Generate Title Text
                            Title = [Current_Salt ' ' Structure ' ' Theory disptext gCP_label ...
                                ' with ' Current_basis ' Basis Set and Cell Vector |a| = ' num2str(a(k))];

                            param_dir = ['APAR_' num2str(a(k),'%05.2f') SymmetryMod];

                            % create/change dir
                            dir3 = [workDir filesep Theory_dir filesep Basis_dir filesep param_dir];
                            if ~exist(dir3,'dir')
                                mkdir(dir3);
                            end
                            cd(dir3)

                            % Replace strings Unique to Unit Cell with Specific Length
                            input_text_a = strrep(input_text_a,'##TITLE##',Title);
                            input_text_a = strrep(input_text_a,'##APAR##',a_param);
                            input_text_a = strrep(input_text_a,[newline '##CELLOPT##'],'');
                            
                            % Check if the job is already in progress
                            if Skip_Running && strcmpi(server(1:3),'sea')
                                JobRunning = CheckRunning(dir3,Stop_Running);
                                if JobRunning
                                    continue
                                end
                            end

                            % Save file
                            fid2 = fopen( [Theory_dir '_' Current_basis '-' num2str(a(k),'%05.2f') '.d12'],'wt');
                            fwrite(fid2,input_text_a);
                            fclose(fid2);

                            cd(dir2) % Return to outer directory
                        end
                    end
                
                    % Submit job
                    if ~isempty(a)
                        if Glost_On

                            % Current job info
                            current_job = [Theory_dir '_' Current_basis optTag];
                            inputfile = fullfile(dir3,[current_job '.d12']);                 
                            crystal_tempdir = fullfile('/scratch','$USER','$SLURM_JOB_ID',current_job);
                            temp_inputfile = fullfile(crystal_tempdir,[current_job '.d12']);

                            % Generate job line
                            Glost_newline = ['export PATH=$PATH && export CRY17_SCRDIR=' crystal_tempdir ' && mkdir -p ' crystal_tempdir ...
                                ' && cd ' crystal_tempdir ' && cp ' inputfile ' ' temp_inputfile ...
                                ' && runcry17 ' current_job ...
                                ' && cp ' current_job '.out ' dir3 filesep current_job '.out' ...
                                ' && cd ' dir3 ' && rm -r ' crystal_tempdir filesep '*' ...
                                ' && rmdir ' crystal_tempdir ' && echo "Job finished at" date' ...
                                ' > ' current_job '.glostlog'];

                            % Add new command to command list
                            if ~exist('Glost_txt','var')
                                Glost_txt = Glost_newline;
                            else
                                Glost_txt = [Glost_txt newline Glost_newline];
                            end

                        else
                            if submit
                                submit_txt = ' 1';
                            else
                                submit_txt = ' 0';
                            end
                            
                            command = [submit_Crystal ' -1 ' nTasks ' ' nCores_per_Node ...
                            ' -1 ' mempernode ' ' hours_calc ' ' mins_calc ' ' Theory_dir '_' Current_basis optTag...
                            ' ' Theory_dir '_' Current_basis optTag ' ' Theory_dir '_' Current_basis optTag ...
                            ' ' links ' ' Restart_txt ' ' Restart_WF_txt params submit_txt];
                            disp([Theory_dir ' ' Current_basis display_text])
                            if ispc
                                disp(command); % for debugging on pc
                            else
                                system(command);
                            end
                        end
                    end
                    cd(dir1)
                end
            end
            cd(workDir)
        end
        cd(Maindir)
    end
end

%% Run Final Glost command
if Glost_On
    
    % Move to GLOST Jobs folder
    Glost_Directory = [Maindir filesep 'GLOST_INPUTS'];
    if ~exist(Glost_Directory,'dir')
        mkdir(Glost_Directory)
    end
    cd(Glost_Directory)
    
    % Save jobslist file
    FID = fopen([Glost_Name '.joblist'],'w');
    fprintf(FID,'%s',Glost_txt);
    fclose(FID);
    
    % Make submission script
    Glost_submission = [ ...
        '#!/bin/bash' newline ...
        '#SBATCH --time=' num2str(hours) ':0:00' newline ...
        '#SBATCH --nodes=' num2str(Glost_Nodes) newline ...
        '#SBATCH --tasks-per-node=' num2str(Glost_Cores_Per_Node) newline ...
        '#SBATCH --mem-per-cpu=' mempernode newline ...
        '#SBATCH --account=' Account newline ...
        '#SBATCH --job-name=' Glost_Name newline ...
        '#SBATCH --error=' Glost_Name '.stde' newline ...
        '#SBATCH --export=ALL' newline ...
        '#SBATCH --exclusive' newline newline ...
        '# set EXE environment' newline ...
        'PATH=$PATH' newline ...
        'export PATH' newline newline ...
        '# Load GLOST module along with the modules required' newline ...
        'module load nixpkgs/16.09  intel/2016.4  openmpi/2.0.2 glost/0.3.1' newline newline ...
        'echo "Starting run at: `date`"' newline newline ...
        '# Move to GLOST start directory' newline ...
        'cd ' Glost_Directory newline newline ...        
        '# Run GLOST with the job list' newline ...
        'mpiexec glost_launch ' Glost_Name '.joblist' newline ....
        'echo "Job finished at"' newline ...
        'date' newline ...
        '################### Job Ended ###################' newline ...
        'exit 0'];
    
    % Save submission script
    FID = fopen([Glost_Name '.subm'],'w');
    fprintf(FID,'%s',Glost_submission);
    fclose(FID);

    % Submit job
    command = ['sbatch ' Glost_Name '.subm'];
    disp(['Submitting New GLOST Job Set: ' Glost_Name '...']);
    if submit && ~ispc
        system(command);
        disp('Job Sumitted.');
    else
        disp(command); % for debugging
    end
end