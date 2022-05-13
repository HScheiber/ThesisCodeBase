function PlotFail = Compare_PES_Structures(app)
% Turn off warning for adding specific sheet
warning('off','MATLAB:xlswrite:AddSheet');

%% Load Plot Figures and pull up load bar
LoadbarObj = waitbar(0,'Loading Settings...','windowstyle','modal');
Mainfigh = app.figureh;
CrosshairButton = app.CrosshairButton;

%% Load Plot Settings
PrintTable = app.settings.PrintTable;
PreviewTable = app.settings.PreviewTable;
FormattedTable = app.settings.FormattedTables;
Table_Filename = app.settings.Table_Filename;
Color_Scheme = app.settings.Color_Scheme;
Basis_Sets = app.settings.Basis_set;
Salts = app.settings.Salts;
Structures = app.settings.Structures;
fs = app.settings.fs;
lw = app.settings.lw;
X_Axis = app.settings.x_axis;
plot_Time = app.settings.plot_Time;
plot_BSSE = app.settings.plot_BSSE;
BSSE_Correct = app.settings.BSSE_Correct;
X_align = app.settings.x_align;
Y_align = app.settings.E_align;
plot_style = app.settings.plot_style;
Y_Axis = app.settings.Y_Axis_Units;
plotTheories = app.settings.plotDFT;
Scaling_Type = app.settings.ScalingType;
Scaling_Params = str2double(app.settings.Scaling);
Data_Types = app.settings.Data_Types;
Emp_Damp = app.settings.Emp_Damp;
Auto_Plot = app.settings.Auto_Plot;
plot_Experimental = app.settings.plot_Experimental;
Interpolation_Active = app.settings.InterpolationSwitch;
Interpolation_Method = app.settings.InterpolationMethod;
InterpStepSize = app.settings.InterpolationSize;

%% Clear old axes, shift focus to new axis
waitbar(0.25,LoadbarObj,'Loading Settings...');
arrayfun(@delete,findall(Mainfigh,'type','Axes'))
arrayfun(@delete,findall(Mainfigh,'type','Text'))
arrayfun(@delete,findall(Mainfigh,'type','Legend'))
arrayfun(@delete,findall(Mainfigh,'type','UIControl'))

% Find colour type
if contained_in_cell(Color_Scheme,{'Accent' 'Dark2' 'Paired' 'Pastel1' 'Pastel2' 'Set1' 'Set2' 'Set3'})
    Ctype = 'qual';
elseif contained_in_cell(Color_Scheme,{'BrBG' 'PiYG' 'PRGn' 'PuOr' 'RdBu' 'RdGy' 'RdYlBu' 'RdYlGn'})
    Ctype = 'div';
else
    Ctype = 'seq';
end

if PreviewTable
    PrintTable = true;
end

% Override x_axis and y_axis control in case of printing table
if PrintTable
    X_Axis = "a";
    plot_BSSE = false;
    if strcmp(Y_Axis,'P')
        Y_Axis = 'L';
    end
end

PlotFail = false; % Initialize as false

% Get scaling type tag
switch lower(Scaling_Type)
    case 'disp'
        Scaling_Tag = 'D';
    case 'rep'
        Scaling_Tag = 'R';
    case 'epsilon'
        Scaling_Tag = 'E';
    case 'sigma'
        Scaling_Tag = 'S';
    case 'mmd'
        Scaling_Tag = 'MMD';
    case 'xxd'
        Scaling_Tag = 'XXD';
    case 'mxd'
        Scaling_Tag = 'MXD';
end

% Generate list of scaling tags for empirical dispersion
NScale = length(Scaling_Params);
Scaling_Tags = cell(NScale,1);
for i = 1:NScale
    if Scaling_Params(i) == 1
        Scaling_Tags{i} = '';
    else
    	Scaling_Tags{i} = ['_' Scaling_Tag regexprep(num2str(Scaling_Params(i),...
            '%7.5f'),{'-' '\.'},{'N' 'P'})];
    end
end

% Generate list of damping tokens for empirical damping
NDamp = length(Emp_Damp);
Dampnum = find(contains({'No Damping' 'BJ Damping' ...
    'Tang Damping' 'MMDRE Damping' 'MANoC Damping' 'EHFSK Damping' 'WY Damping'},Emp_Damp))-1;
Damp_numc = arrayfun(@num2str, Dampnum, 'UniformOutput', 0);
[dmp{1:NDamp}] = deal('d');
Damp_Tokens = strrep(join({dmp{:}; Damp_numc{:}}',''),'d0','');

%% Permanent settings
waitbar(0.50,LoadbarObj,'Loading Settings...');
home = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase';
NStr = length(app.settings.Structures);
Structures = cell(2,NStr);
Structures(1,:) = app.settings.Structures;
Basis_Sets = strrep(Basis_Sets,'CRYSTAL Old Basis','7-311G');

% Get Labels
for i=1:NStr
    if strcmp(Structures{1,i},'Rocksalt')
        Structures{2,i} = 'R';
    elseif strcmp(Structures{1,i},'Wurtzite')
        Structures{2,i} = 'W';
    elseif strcmp(Structures{1,i},'Sphalerite')
        Structures{2,i} = 'S';
    elseif strcmp(Structures{1,i},'CsCl')
        Structures{2,i} = 'C';
    elseif strcmp(Structures{1,i},'NiAs')
        Structures{2,i} = 'N';
    elseif strcmp(Structures{1,i},'BetaBeO')
        Structures{2,i} = 'B';
    elseif strcmp(Structures{1,i},'FiveFive')
        Structures{2,i} = 'F';
    elseif strcmp(Structures{1,i},'Pair')
        Structures{2,i} = 'P';
    end
end

Data_type_code = [];
for i=1:length(Data_Types)
    if strcmp(Data_Types{i},'Rigid Structure')
        Data_type_code(end+1) = 0;
    elseif strcmp(Data_Types{i},'Cell Opt')
        Data_type_code(end+1) = 1;
    elseif strcmp(Data_Types{i},'Full Opt')
        Data_type_code(end+1) = 2;
    elseif strcmp(Data_Types{i},'Full Opt (No Sym)')
        Data_type_code(end+1) = 4;
    elseif strcmp(Data_Types{i},'Full Opt (P = 1 atm)')
        Data_type_code(end+1) = 5;
    elseif strcmp(Data_Types{i},'Atom Optimized')
        Data_type_code(end+1) = 3;
    end
end
if isempty(Data_type_code)
    error('No Data Type chosen!')
end

%% Unit Conversions
waitbar(0.75,LoadbarObj,'Loading Settings...');
Ha_kJ = 4.35974e-21; % KiloJoules per Hartree
NA = 6.0221409e23; % Molecules per mole
Ha_kJ_NA = Ha_kJ*NA; % KiloJoule per Hartree per mole
m3_per_A3 = 1e-30; % Cubic meter per cubic angstrom
cm3_per_A3 = 1e-24; % Cubic centimeter per cubic angstrom
GJ_per_kJ = 1e-6; % GigaJ per kJ

% If there is two plots
plot_extra = plot_Time || plot_BSSE;

% Lattice energy or cohesive energy
if strcmp(Y_Axis,'L')
    energy_type = 'Lattice';
    energy_title = 'Lattice Energy';
elseif strcmp(Y_Axis,'C')
    energy_type = 'Cohesive';
    energy_title = 'Cohesive Energy';
elseif strcmp(Y_Axis,'P')
    energy_type = 'Lattice';
    energy_title = 'Pressure';
end

quantum_data_dir = [home filesep 'data' filesep 'CRYSTAL'];

empirical_data_dir = fullfile(home,'data','GROMACS');
experimental_data_dir = fullfile(home,'data','CRYSTAL');
MD_data_dir = fullfile(home,'data','PURE_SALT_MD');
waitbar(1,LoadbarObj,'Loading Settings...');

%% Load Data
waitbar(0,LoadbarObj,'Loading Quantum Data...');
Nsalt = length(Salts);
Nstrct = size(Structures,2);
NBasis = length(Basis_Sets);
Salt_Data = struct();
if plot_BSSE
    Metal_BSSE_Data = struct();
    Halide_BSSE_Data = struct();
end

%% Load quantum data
for i=1:Nstrct
    waitbar(i/Nstrct,LoadbarObj,'Loading Quantum Data...');
    Structure = Structures{1,i};
    Salt_Data.(Structure) = struct();

    for j=1:Nsalt
        Salt = Salts{j};
    
        try   
            if BSSE_Correct      
                X = load([quantum_data_dir filesep Salt '_' Structure '_BSSE_Corrected_' energy_type '_Energies.mat'],'BSSE_Corrected_Data');
                Salt_Data.(Structure) = catstruct(Salt_Data.(Structure),X.BSSE_Corrected_Data);
                clearvars X 
            else
                X = load([quantum_data_dir filesep Salt '_' Structure '_' energy_type '_Energies.mat'],'Data');
                Salt_Data.(Structure) = catstruct(Salt_Data.(Structure),X.Data);
                clearvars X
            end
        catch
            %error(['Error: No Quantum Data for ' Salt ' ' energy_type ' PES available.'])
        end
    end
end

%% Remove empty fields
for i=1:Nstrct
Structure = Structures{1,i};
    Filedstemp = fields(Salt_Data.(Structure));
    for j = 1:numel(Filedstemp)
        if isempty(Salt_Data.(Structure).(Filedstemp{j}))
            Salt_Data.(Structure) = rmfield(Salt_Data.(Structure),Filedstemp{j});
        end
    end
end
waitbar(0,LoadbarObj,'Loading Empirical Data...');

%% Load all empirical data from storage
Emp_Salt_Data = struct();
Contains_Emp = regexp(plotTheories,'TF|JC');
if sum([Contains_Emp{:}]) > 0
    if strcmp(energy_type,'Lattice')
        for j=1:Nsalt
            waitbar(j/Nsalt,LoadbarObj,'Loading Empirical Data...');
            %Emp_Salt_Data.(Salt) = struct();
            Salt = Salts{j};
            for i=1:Nstrct
                Structure = Structures{1,i};

                try
                    X = load(fullfile(empirical_data_dir,[Salt '_' Structure '_Lattice_Energies.mat']),'Data');
                    if ~isempty(X.Data)
                        Emp_Salt_Data.(Salt).(Structure) = X.Data.(Salt).(Structure);
                        clearvars X
                    end
                catch
                    warning(['Unable to load file: ' fullfile(empirical_data_dir,[Salt '_' Structure '_Lattice_Energies.mat. Skipping.'])]);
                    continue
                end
            end
        end
    end
end

%% Remove empty fields in BSSE data
waitbar(0,LoadbarObj,'Selecting Data...');
if plot_BSSE
    if isempty(Metal_BSSE_Data) || isempty(Halide_BSSE_Data)
        warning('No BSSE data, or incomplete BSSE Data for selected data.')
        plot_BSSE = false;
    else
        Filedstemp = fields(Metal_BSSE_Data);
        for i = 1:numel(Filedstemp)
            if isempty(Metal_BSSE_Data.(Filedstemp{i}))
                Metal_BSSE_Data = rmfield(Metal_BSSE_Data,Filedstemp{i});
            end
        end

        Filedstemp = fields(Halide_BSSE_Data);
        for i = 1:numel(Filedstemp)
            if isempty(Halide_BSSE_Data.(Filedstemp{i}))
                Halide_BSSE_Data = rmfield(Halide_BSSE_Data,Filedstemp{i});
            end
        end
    end
end

%% Find unique list of all theory types used in quantum calculations and empirical calculations
waitbar(0.25,LoadbarObj,'Selecting Data...');
Qfields_Str = {};
Efields_Str = {};
for i=1:Nstrct
    Structure = Structures{1,i};
    Symbol = Structures{2,i};
    
    for j=1:Nsalt
        Salt = Salts{j};
        % Quantum
        DataFields = fieldnames(Salt_Data.(Structure));
        
        ind = ~cellfun(@isempty,regexp(DataFields,Salt));

        Qfields_Str = union(Qfields_Str,...
            unique(strrep(DataFields(ind),...
            [Salt '_' Symbol '_'],'')));
        
        % Empirical
        if isfield(Emp_Salt_Data,Salt)
            if isfield(Emp_Salt_Data.(Salt),Structure)
                EmpDataFields = fieldnames(Emp_Salt_Data.(Salt).(Structure));

                Efields_Str = unique(union(Efields_Str,EmpDataFields));
            end
        end
    end
end

%% Separate quantum from empirical models
Empindx = contains(plotTheories,{'JC' 'TF'});
EmpTheories = plotTheories(Empindx);
QTheories = plotTheories(~Empindx);
NEmpTher = length(EmpTheories);

% Replace different TF models with their tags
EmpTheories = strrep(EmpTheories,'_D3','1');
EmpTheories = strrep(EmpTheories,'_Update','2');
EmpTheories = strrep(EmpTheories,'_D4','3');

Empirical_Theories = cell(1,NEmpTher*NScale*NDamp);
x=0;
for i = 1:NEmpTher 
    for j = 1:NDamp
        for k = 1:NScale
            x=x+1;
            Empirical_Theories{x} = [EmpTheories{i} Damp_Tokens{j} Scaling_Tags{k}];
        end
    end
end

    
%% Modify list to only contain those that match selected theories
waitbar(0.60,LoadbarObj,'Selecting Data...');
QuantumTheories = intersect(Qfields_Str,QTheories);
EmpiricalTheories = intersect(Efields_Str,Empirical_Theories);

% Number of theories used
N_Quantum = length(QuantumTheories);
N_Empirical = length(EmpiricalTheories);
N_emp = N_Empirical*Nstrct*Nsalt;
N_QM = NBasis*N_Quantum*Nstrct*Nsalt;
N_Col = N_QM + N_emp;

if plot_Experimental
    N_Col = N_Col + Nsalt;
end

if N_Col == 0 % No data selected
    msgbox('Warning: No Data Selected.')
    PlotFail = true;
    return
elseif N_Col < 3
    N_Col = 3;
end

%% Pre-assign sub_type array and other variables
waitbar(0.8,LoadbarObj,'Selecting Data...');
Legend_Info = cell(N_Col,1);
Time_list = nan(N_Col,1);
if ~PrintTable
    if plot_extra 
        if strcmp(plot_style,'Scatter') || strcmp(plot_style,'Line')
            axes1 = axes(Mainfigh,'Position',[0.056470588235294 0.574718668268892 0.783340336134454 0.382728140241747]);
        elseif strcmp(plot_style,'Bar')
            axes1 = axes(Mainfigh,'Position',[0.0494791666666667 0.574718668268892 0.95 0.382728140241747]);
        end
    else
        axes1 = axes(Mainfigh,'Position',[0.0625 0.11318795430945 0.730042016806722 0.81619937694704]);
    end
    hold(axes1,'on');
end

% Initialize counters
waitbar(1,LoadbarObj,'Selecting Data...');
q=0; % Keeps track of each item to plot
p = cell(N_Col,1);
if X_align || PrintTable
    struct_data = cell(N_Col,16);
end
if Y_align || PrintTable
    Ymin = cell(N_Col,1);
end

if (strcmp(plot_style,'Bar') || plot_Time) && ~PrintTable
    % Initialize bar plot stuff
    Bar_Y = nan(N_Col,1);
    Bar_Label = cell(N_Col,1);
    Bar_Info = cell(N_Col,4);
end

%% Plot Quantum Data from CRYSTAL17
waitbar(0,LoadbarObj,'Plotting Quantum Data...');
for g = 1:Nsalt % Loop through salt types
    Salt = Salts{g};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    
    for i = 1:N_Quantum % Loop through levels of theory
        DFTmethod = QuantumTheories{i};
        
        for h = 1:Nstrct % Loop through structure types
            Structure = Structures{1,h};
            Symbol = Structures{2,h};
            Method = [Salt '_' Symbol '_' DFTmethod];
            
%             if strcmp(Structure,'Rocksalt')
%                 indh = h;
%             end

            % Check if this DFT method is available in data for specified
            % structure type
            if isfield(Salt_Data.(Structure), Method)

                B = size(Salt_Data.(Structure).(Method),1); % Number of available basis sets

                % Loop through chosen basis sets
                for j=1:NBasis

                    % Current basis set
                    Basis_Set = Basis_Sets{j};

                    % Loop through available basis sets
                    for k=1:B

                        % Current Basis set
                        Basis = Salt_Data.(Structure).(Method){k,1};

                        % If current basis set matches desired, plot it
                        if strcmpi(Basis,Basis_Set)
                            
                            % Get structural info
                            a = Salt_Data.(Structure).(Method){k,2}{2};
                            b = Salt_Data.(Structure).(Method){k,2}{3};
                            c = Salt_Data.(Structure).(Method){k,2}{4};
                            BL = Salt_Data.(Structure).(Method){k,2}{5};
                            Volume = Salt_Data.(Structure).(Method){k,2}{6};
                            FC = Salt_Data.(Structure).(Method){k,2}{7}; % Frac coords
                            Energy = Salt_Data.(Structure).(Method){k,2}{8}.*Ha_kJ_NA; % kJ/mol
                            DT = Salt_Data.(Structure).(Method){k,2}{10}; % data type
                            if strcmpi(Structure,'Pair')
                                alpha = zeros(1,length(DT));
                                beta = zeros(1,length(DT));
                                gamma = zeros(1,length(DT));
                            else
                                alpha = Salt_Data.(Structure).(Method){k,2}{11}; % alpha angle
                                beta = Salt_Data.(Structure).(Method){k,2}{12}; % beta angle
                                gamma = Salt_Data.(Structure).(Method){k,2}{13}; % gamma angle
                            end

                            % Remove unselected data
                            DatInd = ismember(DT,Data_type_code);
                            if ~DatInd % If no data meets criteria, move on
                                continue
                            end
                            a = a(DatInd);
                            b = b(DatInd);
                            c = c(DatInd);
                            BL = BL(DatInd);
                            Volume = Volume(DatInd); % Cubic Ang per Unit Cell
                            FC = FC(DatInd);
                            Energy = Energy(DatInd); % kJ/mol
                            DT = DT(DatInd);
                            alpha = alpha(DatInd); % alpha angle
                            beta = beta(DatInd); % beta angle
                            gamma = gamma(DatInd); % gamma angle
                            
                            if ismember(2,Data_type_code) && ~isempty(Salt_Data.(Structure).(Method){k,6})
                                Thermdat = Salt_Data.(Structure).(Method){k,6}{end,21};
                                N_SuperCell = Salt_Data.(Structure).(Method){k,6}{end,1};
  
                                % Therm Data at lowest pressure
                                T = Thermdat(:,1);
                                a_f = Thermdat(:,2);
                                b_f = Thermdat(:,4);
                                c_f = Thermdat(:,6);
                                
                                a_Therm = interp1(T,a_f,298.15,'spline')/N_SuperCell;
                                b_Therm = interp1(T,b_f,298.15,'spline')/N_SuperCell;
                                c_Therm = interp1(T,c_f,298.15,'spline')/N_SuperCell;
                                if strcmp(Structure,'Rocksalt') || strcmp(Structure,'Sphalerite')
                                    a_Therm = a_Therm*sqrt(2);
                                    b_Therm = b_Therm*sqrt(2);
                                    c_Therm = c_Therm*sqrt(2);
                                end
                            else
                                a_Therm = nan;
                                b_Therm = nan;
                                c_Therm = nan;
                            end
                            
                            % unit conversion for pressure calculations
                            if strcmp(Y_Axis,'P') || X_Axis == "P"
                                Vol_m3 = (Volume.*m3_per_A3./Atom_Pairs_Per_Unit_Cell(Structure)).*NA; % m^3 per mol
                                Energy_GJ = Energy.*GJ_per_kJ; % GJ/mol
                            end
                            
                            if Interpolation_Active
                                if strcmp(Structure,'Pair') && X_Axis ~= "bl"
                                    continue
                                end
                                if X_Axis == "a"
                                    % Generate finer grid
                                    a_interp = a(1):InterpStepSize:a(end);
                                    if strcmp(Y_Axis,'P')
                                        % Interpolate the volume and energy to match
                                        Vol_interp = interp1(a,Vol_m3,a_interp,Interpolation_Method);
                                        Energy_interp = interp1(a,Energy_GJ,a_interp,Interpolation_Method);
                                        
                                        % Take numerical derivative
                                        P = (-diff(Energy_interp)./diff(Vol_interp)); % In GigaPascales
                                        
                                        X_uncor = a_interp(1:end-1);
                                        Y_uncor = P;
                                    else
                                        Energy_interp = interp1(a,Energy,a_interp,Interpolation_Method);
                                        X_uncor = a_interp;
                                        Y_uncor = Energy_interp;
                                    end
                                elseif X_Axis == "bl"
                                    % Generate finer grid
                                    BL_interp = BL(1):InterpStepSize:BL(end);
                                    if strcmp(Y_Axis,'P')                                        
                                        % Interpolate the volume and energy to match
                                        Vol_interp = interp1(BL,Vol_m3,BL_interp,Interpolation_Method);
                                        Energy_interp = interp1(BL,Energy_GJ,BL_interp,Interpolation_Method);
                                        
                                        % Take numerical derivative
                                        P = (-diff(Energy_interp)./diff(Vol_interp)); % In GigaPascales
                                        
                                        X_uncor = BL_interp(1:end-1);
                                        Y_uncor = P;
                                    else
                                        Energy_interp = interp1(BL,Energy,BL_interp,Interpolation_Method);
                                        X_uncor = BL_interp;
                                        Y_uncor = Energy_interp;
                                    end
                                elseif X_Axis == "V"
                                    % Unit conversion
                                    Vol_cm3_mol = (Volume.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(Structure)).*NA; % cm^3 per mol
                                    
                                    % Generate finer grid
                                    Vol_cm3_mol_interp = Vol_cm3_mol(1):InterpStepSize:Vol_cm3_mol(end);
                                    if strcmp(Y_Axis,'P')                                        
                                        % Interpolate the volume and energy to match
                                        Vol_interp = interp1(Vol_cm3_mol,Vol_m3,Vol_cm3_mol_interp,Interpolation_Method);
                                        Energy_interp = interp1(Vol_cm3_mol,Energy_GJ,Vol_cm3_mol_interp,Interpolation_Method);
                                        
                                        % Take numerical derivative
                                        P = (-diff(Energy_interp)./diff(Vol_interp)); % In GigaPascales
                                        
                                        X_uncor = Vol_cm3_mol_interp(1:end-1);
                                        Y_uncor = P;
                                    else
                                        Energy_interp = interp1(Vol_cm3_mol,Energy,Vol_cm3_mol_interp,Interpolation_Method);
                                        X_uncor = Vol_cm3_mol_interp;
                                        Y_uncor = Energy_interp;
                                    end
                                elseif X_Axis == "P"
                                    % Generate finer grid
                                    Vol_interp = Vol_m3(1):InterpStepSize/1e5:Vol_m3(end);
                                    
                                    Energy_interp = interp1(Vol_m3,Energy_GJ,Vol_interp,Interpolation_Method);
                                    P = (-diff(Energy_interp)./diff(Vol_interp)); % In GigaPascales
                                    
                                    if strcmp(Y_Axis,'P')                                                                               
                                        X_uncor = P;
                                        Y_uncor = P;
                                    else
                                        Energy_interp = Energy_interp./GJ_per_kJ; % Convert back to kJ/mol;
                                        X_uncor = P;
                                        Y_uncor = Energy_interp(1:end-1);
                                    end
                                elseif X_Axis == "d"
                                    % Mass
                                    MassMetal = elements('Symbol',Metal,'atomic_mass'); % g/mol
                                    MassHalide = elements('Symbol',Halide,'atomic_mass'); % g/mol
                                    Molar_Mass = MassMetal + MassHalide; % g/mol
                                    % Volume
                                    Vol_cm3 = (Volume.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(Structure)).*NA; % cm^3 / mol
                                    
                                    Density = Molar_Mass./Vol_cm3; % g / cm^3
                                    
                                    % Generate finer grid
                                    Density_interp = Density(1):-InterpStepSize:Density(end);
                                    if strcmp(Y_Axis,'P')                                        
                                        % Interpolate the volume and energy to match
                                        Vol_interp = interp1(Density,Vol_m3,Density_interp,Interpolation_Method);
                                        Energy_interp = interp1(Density,Energy_GJ,Density_interp,Interpolation_Method);
                                        
                                        % Take numerical derivative
                                        P = (-diff(Energy_interp)./diff(Vol_interp)); % In GigaPascales
                                        
                                        X_uncor = Density_interp(1:end-1);
                                        Y_uncor = P;
                                    else
                                        Energy_interp = interp1(Density,Energy,Density_interp,Interpolation_Method);
                                        X_uncor = Density_interp;
                                        Y_uncor = Energy_interp;
                                    end
                                end
                            % No interpolation
                            else
                                if X_Axis == "a"
                                    if strcmp(Y_Axis,'P')
                                        % Take numerical derivative
                                        P = (-diff(Energy_GJ)./diff(Vol_m3)); % In GigaPascales
                                        
                                        X_uncor = a(1:end-1);
                                        Y_uncor = P;
                                    else
                                        X_uncor = a;
                                        Y_uncor = Energy;
                                    end
                                elseif X_Axis == "bl"  
                                    if strcmp(Y_Axis,'P')
                                        % Take numerical derivative
                                        P = (-diff(Energy_GJ)./diff(Vol_m3)); % In GigaPascales
                                        
                                        X_uncor = BL(1:end-1);
                                        Y_uncor = P;
                                    else
                                        X_uncor = BL;
                                        Y_uncor = Energy;
                                    end
                                elseif X_Axis == "V"
                                    % Unit conversion
                                    Vol_cm3_mol = (Volume.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(Structure)).*NA; % cm^3 per mol
                                    if strcmp(Y_Axis,'P')
                                        % Take numerical derivative
                                        P = (-diff(Energy_GJ)./diff(Vol_m3)); % In GigaPascales
                                        
                                        X_uncor = Vol_cm3_mol(1:end-1);
                                        Y_uncor = P;
                                    else
                                        X_uncor = Vol_cm3_mol;
                                        Y_uncor = Energy;
                                    end
                                elseif X_Axis == "P"
                                    % Take numerical derivative
                                    P = (-diff(Energy_GJ)./diff(Vol_m3)); % In GigaPascales
                                    if strcmp(Y_Axis,'P')
                                        X_uncor = P;
                                        Y_uncor = P;
                                    else
                                        X_uncor = P;
                                        Y_uncor = Energy(1:end-1);
                                    end
                                elseif X_Axis == "d"
                                    % Mass
                                    MassMetal = elements('Symbol',Metal,'atomic_mass'); % g/mol
                                    MassHalide = elements('Symbol',Halide,'atomic_mass'); % g/mol
                                    Molar_Mass = MassMetal + MassHalide; % g/mol
                                    
                                    % Volume
                                    Vol_cm3 = (Volume.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(Structure)).*NA; % cm^3 / mol
                                    
                                    Density = Molar_Mass./Vol_cm3; % g / cm^3
                                    if strcmp(Y_Axis,'P')
                                        % Take numerical derivative
                                        P = (-diff(Energy_GJ)./diff(Vol_m3)); % In GigaPascales
                                        
                                        X_uncor = Density(1:end-1);
                                        Y_uncor = P;
                                    else
                                        X_uncor = Density;
                                        Y_uncor = Energy;
                                    end
                                end
                            end
                            
                            q=q+1;
                            Time_list(q) = Salt_Data.(Structure).(Method){k,3};
                            % Choose font size for basis set
                            BASL = length(Basis_Set);
                            if BASL <= 10
                                bassfont = '10';
                            elseif BASL <= 15
                                bassfont = '9';
                            elseif BASL <= 20
                                bassfont = '8';
                            elseif BASL <= 25
                                bassfont = '7.5';
                            elseif BASL <= 30
                                bassfont = '7';
                            else
                                bassfont = '6';
                            end
                            Basstext = ['\fontsize{' bassfont '}{0}\selectfont{' ...
                                regexprep(regexprep(Basis_Set,'^M',''),'_H','/') '}'];
                            
                            Legend_Info{q} = [Salt '(' Symbol ') ' regexprep(regexprep(DFTmethod,'_','-'),...
                                '(.+)\1+','$1') newline strrep(Basstext,'7-311G','CRYSTAL Basis')];
                            
                            % Realign data before plotting
                            if X_align || PrintTable
                                struct_data{q,1} = DFTmethod;
                                struct_data{q,2} = Structure;
                                [~,Indx] = min(Y_uncor);
                                struct_data{q,3} = X_uncor(Indx);
                                X = X_uncor - struct_data{q,3};
                                struct_data{q,4} = Salt;
                                struct_data{q,5} = b(Indx);
                                struct_data{q,6} = c(Indx);
                                struct_data{q,7} = Volume(Indx);
                                struct_data{q,8} = BL(Indx);
                                struct_data{q,9} = FC{Indx};
                                struct_data{q,10} = alpha(Indx);
                                struct_data{q,11} = beta(Indx);
                                struct_data{q,12} = gamma(Indx);
                                struct_data{q,13} = DT(Indx);
                                struct_data{q,14} = a_Therm;
                                struct_data{q,15} = b_Therm;
                                struct_data{q,16} = c_Therm;
                            else
                                X = X_uncor;
                            end

                            if Y_align || PrintTable
                                Ymin{q} = min(Y_uncor);
                                Y = Y_uncor - Ymin{q};
                            else
                                Y = Y_uncor;
                            end

                            % Plot Data
                            if ~PrintTable
                                hold on
                                if strcmp(plot_style,'Scatter')
                                    p{q} = scatter(axes1,X,Y,40,'k',...
                                        'Marker','o','Visible','on');
                                    p{q}.UserData = {Salt DFTmethod Structure Basis_Set};
                                elseif strcmp(plot_style,'Line')
                                    p{q} = plot(axes1,X,Y,'Color','k',...
                                        'LineWidth',lw,'LineStyle','-',...
                                        'Visible','on');
                                    p{q}.UserData = {Salt DFTmethod Structure Basis_Set};
                                end
                                if strcmp(plot_style,'Bar') || plot_Time
                                    Bar_Y(q) = -min(Y_uncor);
                                    Bar_Label{q} = ['\begin{tabular}{c} ' Salt '(' Symbol ')' ...
                                        ' \\ ' regexprep(DFTmethod,'_','-') ' \\ ' Basis_Set ' \end{tabular}'];
                                    Bar_Info(q,:) = {Salt DFTmethod Structure Basis_Set};
                                end
                            end
                        end 
                    end
                end
            end
        end
    end
    waitbar(g/Nsalt,LoadbarObj,'Plotting Quantum Data...');
end

if ismember(2,Data_type_code)
    % Load thermal expansion data
    ThermDat = load([MD_data_dir filesep 'T298K_Thermal_Expansion_Data.mat'],'Data');
    ThermDat = ThermDat.Data;
else
    % Does not exist for other systems
    ThermDat = struct;
end


%% Plot Empirical Potentials from GROMACS
waitbar(g/Nsalt,LoadbarObj,'Plotting Empirical Data...');
for g=1:Nsalt % Loop through salt types
    ThermSwitch_Salt = true;
    Salt = Salts{g};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    
    % Check to see if salt data is available
    if ~isfield(Emp_Salt_Data,Salt)
        continue
    elseif ~isfield(ThermDat,Salt)
        ThermSwitch_Salt = false;
    end
    
    for h = 1:Nstrct % Loop throught structures

        Structure = Structures{1,h};
        Symbol = Structures{2,h};
        
        switch Symbol
            case 'R'
                alpha = 90;
                beta = 90;
                gamma = 90;
            case 'W'
                alpha = 90;
                beta = 90;
                gamma = 120;
            case 'S'
                alpha = 90;
                beta = 90;
                gamma = 90;
            case 'C'
                alpha = 90;
                beta = 90;
                gamma = 90;
            case 'N'
                alpha = 90;
                beta = 90;
                gamma = 120;
            case 'B'
                alpha = 90;
                beta = 90;
                gamma = 90;
            case 'F'
                alpha = 90;
                beta = 90;
                gamma = 90;
            case 'P'
                alpha = nan;
                beta = nan;
                gamma = nan;
        end
        
        % Check to see if structure data is available
        if ~isfield(Emp_Salt_Data.(Salt),Structure)
            continue
        end
        if ThermSwitch_Salt
            if isfield(ThermDat.(Salt),Structure)
                ThermSwitch_Struc = true;
            else
                ThermSwitch_Struc = false;
            end
        else
            ThermSwitch_Struc = false;
        end
        
        for i=1:N_Empirical % loop through models
            Model = EmpiricalTheories{i};
            Model_legtxt = regexprep(Model,'_(D|E)([0-9])P([0-9]{5})',' \($1 \$\\times\$ $2\.$3\)');
            
            % Check to see if model data is available
            if ~isfield(Emp_Salt_Data.(Salt).(Structure),Model)
                continue
            elseif isempty(Emp_Salt_Data.(Salt).(Structure).(Model))
                continue
            end
            
            if ThermSwitch_Struc
                if isfield(ThermDat.(Salt).(Structure),Model)
                    ThermSwitch_Mod = true;
                else
                    ThermSwitch_Mod = false;
                end
            else
                ThermSwitch_Mod = false;
            end

            Current_data = Emp_Salt_Data.(Salt).(Structure).(Model);           
            
            % Get structural info
            a = [Current_data{:,1}];
            b = [Current_data{:,2}];
            c = [Current_data{:,3}];
            BL = [Current_data{:,4}];
            Volume = [Current_data{:,5}]; % Cubic Angstrom per unit cell
            FC = Current_data(:,6)'; % Frac coords
            DT = [Current_data{:,9}]; % data type
            Energy = [Current_data{:,7}]; % kJ/mol
            Angles = {Current_data{:,8}}; %#ok<CCAT1>

            % Remove unselected data
            DatInd = ismember(DT,Data_type_code);
            if ~DatInd % If no data meets criteria, move on
                continue
            end
            a = a(DatInd);
            b = b(DatInd);
            c = c(DatInd);
            BL = BL(DatInd);
            Volume = Volume(DatInd); % Cubic Ang per Unit Cell
            FC = FC(DatInd);
            Energy = Energy(DatInd); % kJ/mol
            DT = DT(DatInd);
            Angles = Angles(DatInd);
            

            if ismember(2,Data_type_code) && ThermSwitch_Mod
                a_Therm = ThermDat.(Salt).(Structure).(Model)(1);
                b_Therm = ThermDat.(Salt).(Structure).(Model)(2);
                c_Therm = ThermDat.(Salt).(Structure).(Model)(3);
            else
                a_Therm = nan;
                b_Therm = nan;
                c_Therm = nan;
            end
            
            % unit conversion for pressure calculations
            if strcmp(Y_Axis,'P') || X_Axis == "P"
                Vol_m3 = (Volume.*m3_per_A3./Atom_Pairs_Per_Unit_Cell(Structure)).*NA; % m^3 per mol
                Energy_GJ = Energy.*GJ_per_kJ; % GJ/mol
            end
            
            if Interpolation_Active
                if X_Axis == "a"
                    % Generate finer grid
                    a_interp = a(1):InterpStepSize:a(end);
                    if strcmp(Y_Axis,'P')
                        % Interpolate the volume and energy to match
                        Vol_interp = interp1(a,Vol_m3,a_interp,Interpolation_Method);
                        Energy_interp = interp1(a,Energy_GJ,a_interp,Interpolation_Method);

                        % Take numerical derivative
                        P = (-diff(Energy_interp)./diff(Vol_interp)); % In GigaPascales

                        X_uncor = a_interp(1:end-1);
                        Y_uncor = P;
                    else
                        Energy_interp = interp1(a,Energy,a_interp,Interpolation_Method);
                        X_uncor = a_interp;
                        Y_uncor = Energy_interp;
                    end
                elseif X_Axis == "bl"
                    % Generate finer grid
                    BL_interp = BL(1):InterpStepSize:BL(end);
                    if strcmp(Y_Axis,'P')                                        
                        % Interpolate the volume and energy to match
                        Vol_interp = interp1(BL,Vol_m3,BL_interp,Interpolation_Method);
                        Energy_interp = interp1(BL,Energy_GJ,BL_interp,Interpolation_Method);

                        % Take numerical derivative
                        P = (-diff(Energy_interp)./diff(Vol_interp)); % In GigaPascales

                        X_uncor = BL_interp(1:end-1);
                        Y_uncor = P;
                    else
                        Energy_interp = interp1(BL,Energy,BL_interp,Interpolation_Method);
                        X_uncor = BL_interp;
                        Y_uncor = Energy_interp;
                    end
                elseif X_Axis == "V"
                    % Unit conversion
                    Vol_cm3_mol = (Volume.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(Structure)).*NA; % cm^3 per mol

                    % Generate finer grid
                    Vol_cm3_mol_interp = Vol_cm3_mol(1):InterpStepSize:Vol_cm3_mol(end);
                    if strcmp(Y_Axis,'P')                                        
                        % Interpolate the volume and energy to match
                        Vol_interp = interp1(Vol_cm3_mol,Vol_m3,Vol_cm3_mol_interp,Interpolation_Method);
                        Energy_interp = interp1(Vol_cm3_mol,Energy_GJ,Vol_cm3_mol_interp,Interpolation_Method);

                        % Take numerical derivative
                        P = (-diff(Energy_interp)./diff(Vol_interp)); % In GigaPascales

                        X_uncor = Vol_cm3_mol_interp(1:end-1);
                        Y_uncor = P;
                    else
                        Energy_interp = interp1(Vol_cm3_mol,Energy,Vol_cm3_mol_interp,Interpolation_Method);
                        X_uncor = Vol_cm3_mol_interp;
                        Y_uncor = Energy_interp;
                    end
                elseif X_Axis == "P"
                    % Generate finer grid
                    Vol_interp = Vol_m3(1):InterpStepSize/1e5:Vol_m3(end);

                    Energy_interp = interp1(Vol_m3,Energy_GJ,Vol_interp,Interpolation_Method);
                    P = (-diff(Energy_interp)./diff(Vol_interp)); % In GigaPascales

                    if strcmp(Y_Axis,'P')                                                                               
                        X_uncor = P;
                        Y_uncor = P;
                    else
                        Energy_interp = Energy_interp./GJ_per_kJ; % Convert back to kJ/mol;
                        X_uncor = P;
                        Y_uncor = Energy_interp(1:end-1);
                    end
                elseif X_Axis == "d"
                    % Mass
                    MassMetal = elements('Symbol',Metal,'atomic_mass'); % g/mol
                    MassHalide = elements('Symbol',Halide,'atomic_mass'); % g/mol
                    Molar_Mass = MassMetal + MassHalide; % g/mol
                    
                    % Volume
                    Vol_cm3 = (Volume.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(Structure)).*NA; % cm^3 / mol

                    Density = Molar_Mass./Vol_cm3; % g / cm^3

                    % Generate finer grid
                    Density_interp = Density(1):-InterpStepSize:Density(end);
                    if strcmp(Y_Axis,'P')                                        
                        % Interpolate the volume and energy to match
                        Vol_interp = interp1(Density,Vol_m3,Density_interp,Interpolation_Method);
                        Energy_interp = interp1(Density,Energy_GJ,Density_interp,Interpolation_Method);

                        % Take numerical derivative
                        P = (-diff(Energy_interp)./diff(Vol_interp)); % In GigaPascales

                        X_uncor = Density_interp(1:end-1);
                        Y_uncor = P;
                    else
                        Energy_interp = interp1(Density,Energy,Density_interp,Interpolation_Method);
                        X_uncor = Density_interp;
                        Y_uncor = Energy_interp;
                    end
                end
            % No interpolation
            else
                if X_Axis == "a"
                    if strcmp(Y_Axis,'P')
                        % Take numerical derivative
                        P = (-diff(Energy_GJ)./diff(Vol_m3)); % In GigaPascales

                        X_uncor = a(1:end-1);
                        Y_uncor = P;
                    else
                        X_uncor = a;
                        Y_uncor = Energy;
                    end
                elseif X_Axis == "bl"  
                    if strcmp(Y_Axis,'P')
                        % Take numerical derivative
                        P = (-diff(Energy_GJ)./diff(Vol_m3)); % In GigaPascales

                        X_uncor = BL(1:end-1);
                        Y_uncor = P;
                    else
                        X_uncor = BL;
                        Y_uncor = Energy;
                    end
                elseif X_Axis == "V"
                    % Unit conversion
                    Vol_cm3_mol = (Volume.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(Structure)).*NA; % cm^3 per mol
                    if strcmp(Y_Axis,'P')
                        % Take numerical derivative
                        P = (-diff(Energy_GJ)./diff(Vol_m3)); % In GigaPascales

                        X_uncor = Vol_cm3_mol(1:end-1);
                        Y_uncor = P;
                    else
                        X_uncor = Vol_cm3_mol;
                        Y_uncor = Energy;
                    end
                elseif X_Axis == "P"
                    % Take numerical derivative
                    P = (-diff(Energy_GJ)./diff(Vol_m3)); % In GigaPascales
                    if strcmp(Y_Axis,'P')
                        X_uncor = P;
                        Y_uncor = P;
                    else
                        X_uncor = P;
                        Y_uncor = Energy(1:end-1);
                    end
                elseif X_Axis == "d"
                    % Mass
                    MassMetal = elements('Symbol',Metal,'atomic_mass'); % g/mol
                    MassHalide = elements('Symbol',Halide,'atomic_mass'); % g/mol
                    Molar_Mass = MassMetal + MassHalide; % g/mol

                    % Volume
                    Vol_cm3 = (Volume.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(Structure)).*NA; % cm^3 / mol

                    Density = Molar_Mass./Vol_cm3; % g / cm^3
                    if strcmp(Y_Axis,'P')
                        % Take numerical derivative
                        P = (-diff(Energy_GJ)./diff(Vol_m3)); % In GigaPascales

                        X_uncor = Density(1:end-1);
                        Y_uncor = P;
                    else
                        X_uncor = Density;
                        Y_uncor = Energy;
                    end
                end
            end

            q=q+1;
            Time_list(q) = sum([Current_data{:,8}]);
            Legend_Info{q} = [Salt ' ' Model_legtxt ' ' Structure];
            
            % Realign data before plotting
            if X_align || PrintTable
                struct_data{q,1} = Model;
                struct_data{q,2} = Structure;
                [~,Indx] = min(Y_uncor);
                struct_data{q,3} = X_uncor(Indx);
                X = X_uncor - struct_data{q,3};
                struct_data{q,4} = Salt;
                struct_data{q,5} = b(Indx);
                struct_data{q,6} = c(Indx);
                struct_data{q,7} = Volume(Indx);
                struct_data{q,8} = BL(Indx);
                struct_data{q,9} = FC{Indx};
                
                if isnan(sum(Angles{Indx}))
                    struct_data{q,10} = alpha;
                    struct_data{q,11} = beta;
                    struct_data{q,12} = gamma;
                else
                    struct_data{q,10} = Angles{Indx}(1);
                    struct_data{q,11} = Angles{Indx}(2);
                    struct_data{q,12} = Angles{Indx}(3);
                end
                struct_data{q,13} = DT(Indx);
                struct_data{q,14} = a_Therm;
                struct_data{q,15} = b_Therm;
                struct_data{q,16} = c_Therm;
            else
                X = X_uncor;
            end

            if Y_align || PrintTable
                Ymin{q} = min(Y_uncor);
                Y = Y_uncor - Ymin{q};
            else
                Y = Y_uncor;
            end

            % Plot Data
            if ~PrintTable
                hold on
                if strcmp(plot_style,'Scatter')
                    p{q} = scatter(axes1,X,Y,40,'k',...
                        'Marker','o','Visible','on');
                    p{q}.UserData = {Salt Model Structure ''};
                elseif strcmp(plot_style,'Line')
                    p{q} = plot(axes1,X,Y,'Color','k',...
                        'LineWidth',lw,'LineStyle','-',...
                        'Visible','on');
                    p{q}.UserData = {Salt Model Structure ''};
                end
                if strcmp(plot_style,'Bar') || plot_Time
                    Bar_Y(q) = -min(Y_uncor);
                    Bar_Label{q} = ['\begin{tabular}{c} ' Salt '(' Symbol ')'...
                        ' \\ ' regexprep(Model,'_','-') ' \end{tabular}'];
                    Bar_Info(q,:) = {Salt Model Structure ''};
                end
            end
        end
    end
    waitbar(g/Nsalt,LoadbarObj,'Plotting Empirical Data...');
end

%% Plot Experimental points
if plot_Experimental
    waitbar(0,LoadbarObj,'Plotting Experimental Data...');
    
    % Load the data 
    X = load([experimental_data_dir filesep 'Experimental_Lattice_Energies.mat'],'Data2');
    Experimental_Data = X.Data2;
    clearvars X
    
    % Number for Born-Haber cycle data
    for g=1:Nsalt
        Salt = Salts{g};
        [Metal,Halide] = Separate_Metal_Halide(Salt);
        x = 0; % Keep track of first of each salt
        X = nan(1,2);
        Y = nan(1,2);
        
        switch Salt
            case 'LiCl'
                FCW_c = 0.379; % Experimental wurtzite x/c
            case 'LiBr'
                FCW_c = 0.379; % Experimental wurtzite x/c
            otherwise
                FCW_c = 3/8; % ideal wurtzite x/c
        end

        % Get data if exists
        try
            Current_a = Experimental_Data.(Salt).a;
            Current_a_Wurtz = Experimental_Data.(Salt).Wurtzite.a;
            Current_c_Wurtz = Experimental_Data.(Salt).Wurtzite.c;
            a_Err = Experimental_Data.(Salt).da;
            a_Err_Wurtz = Experimental_Data.(Salt).Wurtzite.da;
            c_Err_Wurtz = Experimental_Data.(Salt).Wurtzite.dc;

            Current_E = Experimental_Data.(Salt).E;
            Current_E_Wurtz = Experimental_Data.(Salt).Wurtzite.E;
            E_Err = Experimental_Data.(Salt).dE;
            E_Err_Wurtz = Experimental_Data.(Salt).Wurtzite.dE;
        catch
            continue
        end

        % incriment counter
        x = x + 1;
        if x == 1
            q = q + 2;
        end

        % Don't know the time, leave it at 0
        Time_list(q-1) = 0;
        Time_list(q) = 0;

        % Manipulate data before plotting
        Y_uncor = [Current_E Current_E_Wurtz];
        Y_Err = [E_Err E_Err_Wurtz];
        b = [Current_a Current_a_Wurtz];
        c = [Current_a Current_c_Wurtz];
        BL = [Current_a.*0.5 Current_a_Wurtz.*sqrt(8/3)];
        Volume = [(Current_a.^3) sind(60).*(Current_a_Wurtz.^2).*Current_c_Wurtz]; % A^3 per unit cell
        FC{1} = {upper(Metal) 0 0 0 ; upper(Halide) -1/2 -1/2 -1/2};
        FC{2} = {upper(Metal) 1/3 2/3 FCW_c ; upper(Halide) 1/3 2/3 0};
        
        if X_Axis == "a"
            X_uncor = [Current_a Current_a_Wurtz];
            X_Err = [a_Err a_Err_Wurtz];
        elseif X_Axis == "bl"
            X_uncor = BL;
            X_Err = [a_Err.*0.5 a_Err_Wurtz.*sqrt(8/3)];
        elseif X_Axis == "V"
            % Unit conversion
            Vol_cm3_mol = [(Volume(1)*cm3_per_A3/4).*NA (Volume(2)*cm3_per_A3/2)*NA]; % cm^3 per mol
            
            X_uncor = Vol_cm3_mol; % cm^3 per mol

            X_Rel_Err = sqrt(3*((a_Err/Current_a)^2));
            X_Rel_Err_Wurtz = sqrt(2*((a_Err_Wurtz/Current_a_Wurtz)^2) + (c_Err_Wurtz/Current_c_Wurtz)^2);
            X_Err = [X_Rel_Err*X_uncor(1) X_Rel_Err_Wurtz*X_uncor(2)];
        elseif X_Axis == "P"
            X_uncor = 0.000101325; % 1 atm in GPa
            X_Err = 0;
        elseif X_Axis == "d"
            MassMetal = elements('Symbol',Metal,'atomic_mass'); % g/mol
            MassHalide = elements('Symbol',Halide,'atomic_mass'); % g/mol
            Molar_Mass = MassMetal + MassHalide; % g/mol

            % Volume
            Vol_cm3_mol = [(Volume(1)*cm3_per_A3/4).*NA (Volume(2)*cm3_per_A3/2)*NA]; % cm^3 / mol

            Density = Molar_Mass./Vol_cm3_mol; % g / cm^3
            X_uncor = Density;

            X_Rel_Err = sqrt(3*((a_Err/Current_a)^2));
            X_Rel_Err_Wurtz = sqrt(2*((a_Err_Wurtz/Current_a_Wurtz)^2) + (c_Err_Wurtz/Current_c_Wurtz)^2);
            X_Err = [X_Rel_Err*X_uncor(1) X_Rel_Err_Wurtz*X_uncor(2)];
        end

        % Sub-type for legend
        if x == 1
            Legend_Info{q-1} = [Salt '(R) Experiment'];
            Legend_Info{q} = [Salt '(W) Experiment'];
        end

        % Realign data before plotting
        if X_align || PrintTable
            struct_data{q-1,1} = 'Experiment';
            struct_data{q,1} = 'Experiment';
            struct_data{q-1,2} = 'Rocksalt';
            struct_data{q,2} = 'Wurtzite';
            struct_data{q-1,3} = X_uncor(1);
            struct_data{q,3} = X_uncor(2);
            X = [0 0];
            struct_data{q-1,4} = Salt;
            struct_data{q,4} = Salt;
            struct_data{q-1,5} = b(1);
            struct_data{q,5} = b(2);
            struct_data{q-1,6} = c(1);
            struct_data{q,6} = c(2);
            struct_data{q-1,7} = Volume(1);
            struct_data{q,7} = Volume(2);
            struct_data{q-1,8} = BL(1);
            struct_data{q,8} = BL(2);
            struct_data{q-1,9} = FC{1};
            struct_data{q,9} = FC{2};
            struct_data{q-1,10} = 90; % alpha
            struct_data{q,10} = 90; % alpha
            struct_data{q-1,11} = 90; % beta
            struct_data{q,11} = 90; % beta
            struct_data{q-1,12} = 90; % gamma
            struct_data{q,12} = 120; % gamma
            struct_data{q-1,13} = 5; % fully optimized (1 atm) data type code
            struct_data{q,13} = 5; % fully optimized (1 atm) data type code
            struct_data{q-1,14} = nan; % a_therm
            struct_data{q-1,15} = nan; % b_therm
            struct_data{q-1,16} = nan; % c_therm
            struct_data{q,14} = nan; % a_therm
            struct_data{q,15} = nan; % b_therm
            struct_data{q,16} = nan; % c_therm
        else
            X(x) = X_uncor(1);
            X(x+1) = X_uncor(2);
        end

        if Y_align || PrintTable
            Ymin{q-1} = Y_uncor(1);
            Ymin{q} = Y_uncor(2);
            Y = [0 0];
        else
            Y(x) = Y_uncor(1);
            Y(x+1) = Y_uncor(2);
        end

        if ~PrintTable
            hold on
            if strcmp(plot_style,'Scatter') || strcmp(plot_style,'Line')
                p{q-1} = errorbar(X(1),Y(1),Y_Err(1),Y_Err(1),X_Err(1),X_Err(1),'o','MarkerSize',6,...
                    'MarkerEdgeColor','none','MarkerFaceColor','k','Visible','on','Parent',axes1,...
                    'Color','k','LineWidth',1.5);
                p{q} = errorbar(X(2),Y(2),Y_Err(2),Y_Err(2),X_Err(2),X_Err(2),'o','MarkerSize',6,...
                    'MarkerEdgeColor','none','MarkerFaceColor','g','Visible','on','Parent',axes1,...
                    'Color','k','LineWidth',1.5);
                p{q-1}.UserData = {Salt 'Experiment' 'Rocksalt' ''};
                p{q}.UserData = {Salt 'Experiment' 'Wurtzite' ''};
            end
            if strcmp(plot_style,'Bar') || plot_Time
                Bar_Y(q-1) = -Y_uncor(1);
                Bar_Y(q) = -Y_uncor(2);
                Bar_Label{q-1} = ['\begin{tabular}{c} ' Salt '(Rocksalt)'...
                    ' \\  Experiment\end{tabular}'];
                Bar_Label{q} = ['\begin{tabular}{c} ' Salt '(Wurtzite)'...
                                    ' \\  Experiment\end{tabular}'];
                Bar_Info(q-1,:) = {Salt 'ExperimentRocksalt'};
                Bar_Info(q,:) = {Salt 'ExperimentWurtzite'};
            end
        end
        waitbar(g/Nsalt,LoadbarObj,'Plotting Experimental Data...');
    end
end

waitbar(0,LoadbarObj,'Setting Plot Style...');

% Bar plot
if strcmp(plot_style,'Bar') && ~PrintTable
    if isnan(Bar_Y) % No data?
        msgbox('Warning: No Data Currently Selected.')
        PlotFail = true;
        return
    end
    [Bar_Y_Sorted,Bar_SD_Sorted,BarCData,Bar_Label_Sorted] = SetPlotColoursBar(Bar_Y,Bar_Label,Bar_Info,Ctype,Color_Scheme);
    p = bar(axes1,Bar_Y_Sorted,'FaceColor','flat','Visible','on','BarWidth',1,...
        'CData',BarCData);
    hold on
    errorbar(p.XData(Bar_SD_Sorted > 0),Bar_Y_Sorted(Bar_SD_Sorted > 0),...
        Bar_SD_Sorted(Bar_SD_Sorted > 0),'r','LineStyle','none',...
        'LineWidth',lw,'CapSize',30)
    
    barmin = min(Bar_Y);
    barmax = max(Bar_Y);
    p = {p};
end

%% Print out a table of minimum energies if requested
if PrintTable
    waitbar(0,LoadbarObj,'Generating Tables...');
    % Get rid of empty entries
    struct_data = struct_data(~all(cellfun('isempty', struct_data(:,1:4)), 2),:);
    
    % No data?
    if isempty(struct_data)
        msgbox('Warning: No Data Currently Selected.');
        PlotFail = true;
        return
    end
    
    if FormattedTable
        %LoadbarObj = PrintLatexTables(LoadbarObj,struct_data,Ymin,Table_Filename,Data_Types,Basis_Sets);
        LoadbarObj = PrintCSV(LoadbarObj,struct_data,Ymin,Table_Filename,Basis_Sets);
        close(LoadbarObj);
        return
    end
    
    N_Col = size(struct_data,1);
    
    Species = struct_data(:,4);
    % Rename species to order correctly
    Species = strrep(Species,'LiF','AAA');
    Species = strrep(Species,'LiCl','BBB');
    Species = strrep(Species,'LiBr','CCC');
    Species = strrep(Species,'LiI','DDD');
    Species = strrep(Species,'NaCl','EEE');
    
    Structure_types = struct_data(:,2);
    Theory_types = regexprep(struct_data(:,1),'_','-');
    Theory_types = regexprep(Theory_types,'(JC|TF)','ZZZ$1'); %% Make sure empirical are sorted last
    Theory_types = strrep(Theory_types,'Experiment','AAA'); %% Make sure experiment are sorted first
    Energy_min = cell2mat(Ymin(:));
    a_min = cell2mat(struct_data(:,3));
    b_min = cell2mat(struct_data(:,5));
    c_min = cell2mat(struct_data(:,6));
    Volume_min = cell2mat(struct_data(:,7));
    BL_min = cell2mat(struct_data(:,8));
    alpha_min = cell2mat(struct_data(:,10));
    beta_min = cell2mat(struct_data(:,11));
    gamma_min = cell2mat(struct_data(:,12));
    DT_min = struct_data(:,13);
    
    % Fractional coords
    FC_min = cell(N_Col,1);
    Opt_Type = cell(N_Col,1);
    for i=1:N_Col
        FC = struct_data{i,9}';
        for x=1:size(FC,2)
            FC(2:4,x) = num2cell(mod([FC{2:4,x}],1)');
        end
        % Sort by metal first then halide
        FC(1,:) = regexprep(FC(1,:),'L(I|i)','AAA');
        FC(1,:) = regexprep(FC(1,:),'N(A|a)','AAB');
        FC(1,:) = regexprep(FC(1,:),'K','AAC');
        FC(1,:) = regexprep(FC(1,:),'R(B|b)','AAD');
        FC(1,:) = regexprep(FC(1,:),'C(S|s)','AAE');
        [~,SortIdx] = sort(FC(1,:));
        FC = FC(:,SortIdx);
        FC(1,:) = strrep(FC(1,:),'AAA','LI');
        FC(1,:) = strrep(FC(1,:),'AAB','Na');
        FC(1,:) = strrep(FC(1,:),'AAC','K');
        FC(1,:) = strrep(FC(1,:),'AAD','RB');
        FC(1,:) = strrep(FC(1,:),'AAE','CS');
        
        FC = strtrim(sprintf('(%5.4f, %5.4f, %5.4f)\n',FC{2:end,:}));
        FC_min{i} = FC;
        
        switch DT_min{i}
            case 0
                Opt_Type{i} = 'Energy Curve';
            case 1
                Opt_Type{i} = 'Cell';
            case 2
                Opt_Type{i} = 'Full';
            case 3
                Opt_Type{i} = 'Atom';
            case 4
                Opt_Type{i} = 'Full (P1)';
            case 5
                Opt_Type{i} = 'Full (1 atm)';
        end
        
    end
    waitbar(0.25,LoadbarObj,'Generating Table...');
    
    % Merge into table
    Data_Table = table(Species,Structure_types,Theory_types,Opt_Type,...
        Energy_min,a_min,b_min,c_min,alpha_min,beta_min,gamma_min,...
        Volume_min,BL_min,FC_min,...
        'VariableNames',{'Species' 'Structure' 'Theory' 'Optimization' ...
        'E_min_kJmol' 'a_Ang' 'b_Ang' 'c_Ang' 'alpha_Degree' 'beta_Degree' 'gamma_Degree' ...
        'V_Ang3' 'BL_Ang' 'Frac_Coords'});
    
    % Sort Table
    Data_Table = sortrows(Data_Table,{'Species','Structure','Theory' 'Optimization'},{'ascend','ascend','ascend' 'ascend'});
    waitbar(0.5,LoadbarObj,'Generating Table...');
    
    % Rename Species to correct names
    Data_Table{:,1} = strrep(Data_Table{:,1},'AAA','LiF');
    Data_Table{:,1} = strrep(Data_Table{:,1},'BBB','LiCl');
    Data_Table{:,1} = strrep(Data_Table{:,1},'CCC','LiBr');
    Data_Table{:,1} = strrep(Data_Table{:,1},'DDD','LiI');
    Data_Table{:,1} = strrep(Data_Table{:,1},'EEE','NaCl');
    Data_Table{:,3} = strrep(Data_Table{:,3},'AAA','Experiment');
    Data_Table{:,3} = strrep(Data_Table{:,3},'ZZZ','');
    
    N_Col = size(Data_Table,1);
    [~,~,ic] = unique(Data_Table(:,2));
    NColours = max(ic);
    if NColours < 3
        Colourset = cbrewer(Ctype,Color_Scheme,3);
    else
        Colourset = cbrewer(Ctype,Color_Scheme,NColours);
    end
    Rowcolours = Colourset(ic,:);
    waitbar(0.75,LoadbarObj,'Generating Table...');
    
    View_Structure = cell(N_Col,1);
    Click_text = '<span style="text-decoration: underline;"><font color="blue"><b>Click To View</b></font></span>';
    View_Structure(:) = {Click_text};
    
    % Separate Cohesive and Lattice energy by sheet
    if strcmp(Y_Axis,'L')
        Sheet_name = 'Lattice Energy';
    elseif strcmp(Y_Axis,'C')
        Sheet_name = 'Cohesive Energy';
    elseif strcmp(Y_Axis,'P')
        Sheet_name = 'Pressure';
    end
    
    if PreviewTable % Open in figure
        % Add click here buttons to table
        Data_Table = [Data_Table(:,1:14) cell2table(View_Structure)];
        
        VarNames = Data_Table.Properties.VariableNames;
        VarNames = regexprep(VarNames,'(^.)','<HTML><center><b>$1');
        VarNames = regexprep(VarNames,'(.$)','$1</b></center></HTML>');
        VarNames = strrep(VarNames,'kJmol','(kJ/mol)');
        VarNames = strrep(VarNames,'Ang3','(&#8491;<sup>3</sup>)');
        VarNames = strrep(VarNames,'Ang','(&#8491;)');
        VarNames = strrep(VarNames,'Degree','(&#176;)');
        VarNames = strrep(VarNames,'alpha','&#945;');
        VarNames = strrep(VarNames,'beta','&#946;');
        VarNames = strrep(VarNames,'gamma','&#947;');
        VarNames = strrep(VarNames,'Delta_','&Delta;');
        VarNames = strrep(VarNames,'_min','<sub>min</sub>');
        VarNames = strrep(VarNames,'No_Disp_','</b>No Dispersion<b><br/>');
        VarNames = strrep(VarNames,'Disp_','</b>Dispersion<b><br/>');
        VarNames = strrep(VarNames,'_',' ');

        Data_Table{:,14} = strrep(Data_Table{:,14},newline,' ');
        Data_Table{:,14} = strrep(Data_Table{:,14},'+',' ');
        
        % Convert to strings, remove NaN, set background colours
        Data_Cell = table2cell(Data_Table);
        Data_Cell = cellfun(@(M) M(~isnan(M)), Data_Cell, 'Uniform', 0);
        for i=1:size(Data_Cell,1) % Row
            for j=1:size(Data_Cell,2) % Column
                if isnumeric(Data_Cell{i,j})
                    if j == 5
                        Data_Cell{i,j} = num2str(Data_Cell{i,j},'%8.2f');
                    else
                        Data_Cell{i,j} = num2str(Data_Cell{i,j},'%8.4f');
                    end
                end
                if isempty(Data_Cell{i,j})
                    Data_Cell{i,j} = '<HTML><table border=0 width=500 bgcolor="rgb(255,255,255)"><TR><TD> - </TD></TR> </table></HTML>';
                else
                    Data_Cell{i,j} = sprintf(['<HTML><table border=0 width=500 bgcolor="rgb(%f,%f,%f)"><TR><TD>'...
                        Data_Cell{i,j} '</TD></TR> </table></HTML>'], Rowcolours(i,:).*255);
                end
            end
        end
        waitbar(1,LoadbarObj,'Generating Table...');
        
        % Create full screen figure
        TabFigure = figure('Units','Normalized','Outerposition',[0 0 1 1],...
            'WindowState','maximize','Name',Table_Filename,'ToolBar','none');
        
        % Create UITable
        UITable = uitable(TabFigure,'Data',Data_Cell,'Units','Normalized','ColumnName',VarNames,...
            'Position',[0 0 1 1],'CellSelectionCallback', {@myCellSelectionCB,Data_Table,home});
        
        % Display the uitable and get its underlying Java object handle
        warning('off','MATLAB:ui:javaframe:PropertyToBeRemoved')
        jscrollpane = findjobj(UITable);
        jtable = jscrollpane.getViewport.getView;

        % Now turn the JIDE sorting on
%        jtable.setSortable(true);		% or: set(jtable,'Sortable','on');
%        jtable.setAutoResort(true);
%        jtable.setMultiColumnSortable(true);
%         jtable.setPreserveSelectionsAfterSorting(true);
        jtable.setAutoResizeMode(jtable.AUTO_RESIZE_SUBSEQUENT_COLUMNS);
        
%         % Modify column widths of fractional coordinate columns
%         set(UITable,'ColumnWidth',{'auto' 'auto' 'auto' 'auto' 'auto' 'auto' ...
%             'auto' 'auto' 'auto' 266 'auto' 'auto' 'auto' 'auto' 'auto' ...
%             'auto' 'auto' 266 'auto'})
    else % Save table to sheet
        try
            writetable(Data_Table,[Table_Filename '.xlsx'],'FileType','spreadsheet',...
                'WriteRowNames',true,'Sheet',Sheet_name)
            RemoveSheet123([Table_Filename '.xlsx']);
        catch
            disp(['Unable to save "' Table_Filename ...
                '.xlsx": permission denied.'])
        end
    end
    
    close(LoadbarObj);
    return
end

%% Plot options
waitbar(0.25,LoadbarObj,'Setting Plot Style...');
if strcmp(plot_style,'Scatter') || strcmp(plot_style,'Line')
    p = SetPlotColours(p,Ctype,Color_Scheme);
    title(axes1,['Comparison of Alkali Halide Salt ' energy_type ' Energy Surfaces' ],...
        'Interpreter','latex','fontsize',fs)

    set(axes1,'box','on','TickLabelInterpreter','latex');
    set(axes1,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    if X_Axis == "a"
        if X_align
            xlim(axes1,[-3 7]);
            xlabel(axes1,'Lattice Parameter $a - a_{min}$ (\AA)','fontsize',fs,'Interpreter','latex');
        else
            xlim(axes1,[2 9]);
            xlabel(axes1,'Lattice Parameter a (\AA)','fontsize',fs,'Interpreter','latex');
        end
    elseif X_Axis == "bl"
        if X_align
            xlabel(axes1,'Metal-Halide Bond Length $BL - BL_{min}$ (\AA)','fontsize',fs,'Interpreter','latex');
            xlim(axes1,[-2 4]);
        else
            xlabel(axes1,'Metal-Halide Smallest Bond Length (\AA)','fontsize',fs,'Interpreter','latex');
            xlim(axes1,[0 6]);
        end
    elseif X_Axis == "V"
        if X_align
            xlabel(axes1,'Crystal Volume $V - V_{min}$ (cm$^{3}$ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');
            xlim(axes1,[-20 40]);
        else
            xlabel(axes1,'Crystal Volume (cm$^{3}$ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');
            xlim(axes1,[0 60]);
        end
    elseif X_Axis == "P"
        if X_align
            xlabel(axes1,'Zero-T Pressure $P - P_{min}$ (GPa)','fontsize',fs,'Interpreter','latex');
        else
            xlabel(axes1,'Zero-T Pressure (GPa)','fontsize',fs,'Interpreter','latex');
        end
    elseif X_Axis == "d"
        if X_align
            xlabel(axes1,'Density $d - d_{min}$ (g cm$^{-3}$)','fontsize',fs,'Interpreter','latex');
        else
            xlabel(axes1,'Density (g cm$^{-3}$)','fontsize',fs,'Interpreter','latex');
        end
    end

    Legend_Info = Legend_Info(~cellfun('isempty',Legend_Info));
    Legend_Info = cellfun(@(x) strrep(x,'_H','/'),Legend_Info,'UniformOutput',false);
    Legend_Info = cellfun(@(x) strrep(x,[newline 'M'],newline),Legend_Info,'UniformOutput',false);
    
    legend1 = clickableLegend(axes1,Legend_Info,...
        'Interpreter','latex','X_Axis',X_Axis,'Y_Axis',Y_Axis);

    ylim(axes1,'auto');
    YL = ylim(axes1);
    
    if strcmp(Y_Axis,'L')
        if Y_align
            ylabel(axes1,[energy_title ' $E - E_{min}$ (kJ mol$^{-1}$)'],'fontsize',fs,'Interpreter','latex');
        else
            ylabel(axes1,[energy_title ' (kJ mol$^{-1}$)'],'fontsize',fs,'Interpreter','latex');
        end
        
        if YL(1) < -2000
            YL(1) = -2000;
        end
        if YL(2) > 1500
            YL(2) = 1500;
        end
        
    elseif strcmp(Y_Axis,'C')
        if Y_align
            ylabel(axes1,[energy_title ' $E - E_{min}$ (kJ mol$^{-1}$)'],'fontsize',fs,'Interpreter','latex');
        else
            ylabel(axes1,[energy_title ' (kJ mol$^{-1}$)'],'fontsize',fs,'Interpreter','latex');
        end
        
        if YL(1) < -2000
            YL(1) = -2000;
        end
        if YL(2) > 1500
            YL(2) = 1500;
        end

    elseif strcmp(Y_Axis,'P')
        if Y_align
            ylabel(axes1,'Zero-T Pressure $P - P_{min}$ (GPa)','fontsize',fs,'Interpreter','latex');
            %ylim(axes1,[-30 100]);
        else
            ylabel(axes1,'Zero-T Pressure (GPa)','fontsize',fs,'Interpreter','latex');
            %ylim(axes1,[-30 80]);
        end
        
        if YL(1) < -50
            YL(1) = -50;
        end
        if YL(2) > 150
            YL(2) = 150;
        end
    end
    ylim(axes1,sort(YL));
    
elseif strcmp(plot_style,'Bar')
    axis(axes1,'manual')
    
    title(axes1,['Comparison of Minimum ' energy_title],...
        'Interpreter','latex','fontsize',fs)
    
    set(axes1,'box','on','TickLabelInterpreter','latex');
    set(axes1,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    
    set(axes1,'XLim',[0.5 N_Col+0.5])
    set(axes1,'YLim',[barmin*0.99 barmax*1.01])
    ylabel(axes1,[energy_title ' (kJ mol$^{-1}$)'],'fontsize',fs,'Interpreter','latex');
    xticks(axes1,1:q)
    xticklabels(axes1,Bar_Label_Sorted)
end

%% Bar plot for time
waitbar(0.5,LoadbarObj,'Finishing...');
if plot_Time
    axes2 = axes(Mainfigh,'Position',[0.047395833333333,0.100726895119418,0.9515625,0.321910695742471]); %[X Y Width Height]
    hold(axes2,'on');
    
    [Time_list_Trimmed,~,TimeCData,Time_Label] = SetPlotColoursBar(Time_list,...
        Bar_Label,Bar_Info,Ctype,Color_Scheme);
    hold on
    bar(axes2,Time_list_Trimmed./(3600),'FaceColor','flat','CData',TimeCData);
    
    set(axes2,'XTick',1:1:length(Time_Label))
    set(axes2,'xticklabel',Time_Label)

    title(['Total CPU Walltime for Lithium Halide ' ...
        energy_title ' Calculations.'],...
        'Interpreter','latex','fontsize',fs)

    set(axes2,'box','on','TickLabelInterpreter','latex');
    set(axes2,'XMinorTick','off','YMinorTick','off','FontSize',fs);
    ylabel('CPU Time (hours)','fontsize',fs,'Interpreter','latex');
elseif plot_BSSE
    if strcmp(plot_style,'Bar')
        axes2 = axes(Mainfigh,'Position',[0.05625 0.0893042575285566 0.7828125 0.33]);
    else
        axes2 = axes(Mainfigh,'Position',[0.05625 0.0893042575285566 0.7828125 0.355140186915888]);
    end
    hold(axes2,'on');

    % Get list of available ions
    Ions = cell(1,2*Nsalt);
    for i=1:Nsalt
        CurSalt = Salts{i};
        [Ions{i},Ions{i+Nsalt}] = Separate_Metal_Halide(CurSalt);
    end
    IonsList = unique(Ions(~cellfun('isempty',Ions)),'stable');
    
    % Generate title text with dropdown menu
    txt2 = uicontrol(Mainfigh,'Style','text','Units','Normalized','FontSize',fs+1,...
        'Position',[0.210416666666666 0.451751817763903 0.344270833333334 0.038383177043988],...
        'String','Comparison of Estimated Basis Set Superposition for',...
        'FontName','CMU Serif Roman');

    ddm = uicontrol(Mainfigh,'Style','popupmenu','Units','Normalized','FontSize',fs,...
        'Position',[0.542 0.455 0.0604 0.038383177043988],...
        'String',IonsList,'FontName','CMU Serif Roman','Callback',@bgselection);

    
    txt3 = uicontrol(Mainfigh,'Style','text','Units','Normalized','FontSize',fs+1,...
        'Position',[0.6069791666666667 0.451751817763903 0.2020833333333334 0.038383177043988],...
        'String','Ion','FontName','CMU Serif Roman',...
        'HorizontalAlignment','left');
    
    % Modify positions for bar plot
    if strcmp(plot_style,'Bar')
        txt2.Position = [0.210416666666666 0.421637591387994 0.344270833333334 0.038383177043988];
        txt3.Position = [0.6069791666666667 0.421637591387994 0.2020833333333334 0.038383177043988];
        ddm.Position = [0.5435625,0.42696261682243,0.0604,0.038383177043988];
    end
    
    % Initialize plot
    initialize = struct('String',{IonsList},'Value',1);
    legend2 = [];
    bgselection(initialize)
end

% Initialize with all plots hidden if selected
if (strcmp(plot_style,'Scatter') || strcmp(plot_style,'Line')) && ~Auto_Plot
    for i=1:length(p)
        set(p{i},'Visible','off');
        set(findobj('Tag',p{i}.DisplayName),'Visible','off');
    end
end
axis(axes1,'manual')
waitbar(0.75,LoadbarObj,'Setting Plot Style...');

% Update crosshair button
set(CrosshairButton,'OnCallback',{@crosshairfcn,X_Axis,Y_Axis,plot_extra,fs});
CrosshairButton.Enable = 'on';

if strcmp(plot_style,'Bar')
    set(CrosshairButton,'OnCallback',{@crosshairfcnbar,X_Axis,plot_extra});
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

% Make figure visible after fully loading everything
drawnow
if strcmp(plot_style,'Scatter') || strcmp(plot_style,'Line')
    if plot_extra
        set(legend1,'Position',[0.842708333333333 0.496365524402907 0.156770833333333 0.502596053997923]);
        if plot_BSSE
            set(legend2,'Position',[0.842708333333333 0 0.156770833333333 0.496365524402908]);
        end
    else
        set(axes1,'Position',[0.0625 0.11318795430945 0.730042016806722 0.81619937694704]); %#ok<*UNRCH>
        set(legend1,'Position',[0.797395833333333 0.00103842159916926 0.202083333333334 0.99896157840083]);
    end
elseif strcmp(plot_style,'Bar')
    if plot_extra
        if plot_BSSE
            set(legend2,'Position',[0.847 0.00207684319833845 0.152 0.45]); % [ X Y Width Height]
        end
    else
        set(axes1,'Position',[0.0473958333333333 0.11 0.952083333333333 0.83392523364486]);
    end
end
waitbar(1,LoadbarObj,'Setting Plot Style.');
close(LoadbarObj);

function bgselection(SelInfo,~)
    cla(axes2) %#ok<*USENS>
    Current_ion = SelInfo.String{SelInfo.Value};
    
    if contained_in_cell(Current_ion,{'Li' 'Na' 'K' 'Rb'})
        MetalorHalide = 'Metal'; % True = metal; False = halide
    elseif contained_in_cell(Current_ion,{'F' 'Cl' 'Br' 'I'})
        MetalorHalide = 'Halide';
    end
    
    % Load pure BSSE correction plot data
    X = load([quantum_data_dir filesep Current_ion '_BSSE_Energies.mat'],[MetalorHalide '_BSSE_Data']);
    BSSE_Data = X.([MetalorHalide '_BSSE_Data']);
    clearvars X
    
    BSSE_fields = fieldnames(BSSE_Data);

    Theories_BSSE = unique(regexprep(BSSE_fields,'.+?_._',''));

    if strcmp(plotTheories,'all')
        DFTmethods_BSSE = Theories_BSSE;
    else
        DFTmethods_BSSE = intersect(Theories_BSSE,regexprep(plotTheories,'_D[2-3].*$',''));
    end
    NN = length(DFTmethods_BSSE);

    % Pre-assign colors and sub_type arrays
    sub_type_BSSE = cell(NN*NBasis,1);
    time_list_BSSE = zeros(NN*NBasis,1);
    
    qq=0;
    pp = cell(NN*NBasis,1);    
    for aidx = 1:Nsalt % Loop through salt types
        CurSalt_BSSE = Salts{aidx};
        [BMetal,BHalide] = Separate_Metal_Halide(CurSalt_BSSE);
        
        for bidx = 1:NN % Loop through levels of theory
            DFTmethod_BSSE = DFTmethods_BSSE{bidx};
            
            for cidx = 1:Nstrct % Loop through structure types
                CurStructure = Structures{1,cidx};
                CurSymbol = Structures{2,cidx};
                CurMethod = [CurSalt_BSSE '_' CurSymbol '_' DFTmethod_BSSE];
                
                % Check if this DFT method is available in data for specified
                % structure type
                if isfield(BSSE_Data,CurMethod)
                    
                    % Number of available basis sets
                    BB = size(BSSE_Data.(CurMethod),1);
                    
                    % Loop through chosen basis sets
                    for didx=1:NBasis
                        % Current chosen basis set
                        Basis_set_BSSE = Basis_Sets{didx};
                        
                        % Loop through available basis sets
                        for eidx=1:BB
                            
                            % Current Basis set
                            Basis_BSSE = BSSE_Data.(CurMethod){eidx,1};

                            % Match basis set
                            if strcmp(Basis_BSSE,Basis_set_BSSE)
                                
                                % Choose font size for basis set
                                BAL = length(Basis_BSSE);
                                if BAL <= 10
                                    basfont = '10';
                                elseif BAL <= 15
                                    basfont = '9';
                                elseif BAL <= 20
                                    basfont = '8';
                                elseif BAL <= 25
                                    basfont = '7.5';
                                elseif BAL <= 30
                                    basfont = '7';
                                else
                                    basfont = '6';
                                end
                                
                                % Matched
                                qq=qq+1;
                                time_list_BSSE(qq) = BSSE_Data.(CurMethod){eidx,3};
                                sub_type_BSSE{qq} = [CurSalt_BSSE '(' CurSymbol ') ' ...
                                    DFTmethod_BSSE newline '\fontsize{' basfont '}{0}\selectfont{' regexprep(regexprep(Basis_BSSE,'^M',''),'_H','/') '}'];

                                % Get structural info
                                a_BSSE = BSSE_Data.(CurMethod){eidx,2}{2};
                                %b_BSSE = BSSE_Data.(CurMethod){eidx,2}{3};
                                %c_BSSE = BSSE_Data.(CurMethod){eidx,2}{4};
                                BL_BSSE = BSSE_Data.(CurMethod){eidx,2}{5};
                                Vol_BSSE = BSSE_Data.(CurMethod){eidx,2}{6};
                                %FC_BSSE = BSSE_Data.(CurMethod){eidx,2}{7}; % Frac coords
                                Energy_BSSE = BSSE_Data.(CurMethod){eidx,2}{8}.*Ha_kJ_NA; % Energy
                                %DT_BSSE = BSSE_Data.(CurMethod){eidx,2}{10}; % data type
                                
                                if Interpolation_Active
                                    if X_Axis == "a" || X_Axis == "P"
                                        % Generate finer grid
                                        a_BSSE_interp = a_BSSE(1):InterpStepSize:a_BSSE(end);
                                        %Energy_BSSE_interp = interp1(a_BSSE,Energy_BSSE,a_BSSE_interp,Interpolation_Method);
                                        XB_uncor = a_BSSE_interp;
                                        %YB_uncor = Energy_BSSE_interp;
                                    elseif X_Axis == "bl"
                                        % Generate finer grid
                                        BL_BSSE_interp = BL_BSSE(1):InterpStepSize:BL_BSSE(end);
                                        %Energy_BSSE_interp = interp1(BL_BSSE,Energy_BSSE,BL_BSSE_interp,Interpolation_Method);
                                        XB_uncor = BL_BSSE_interp;
                                        %YB_uncor = Energy_BSSE_interp;
                                    elseif X_Axis == "V"
                                        % Unit conversion
                                        Vol_BSSE_cm3_mol = (Vol_BSSE.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(CurStructure)).*NA; % cm^3 per mol

                                        % Generate finer grid
                                        Vol_BSSE_cm3_mol_interp = Vol_BSSE_cm3_mol(1):InterpStepSize:Vol_BSSE_cm3_mol(end);
                                        %Energy_BSSE_interp = interp1(Vol_BSSE_cm3_mol,Energy_BSSE,Vol_BSSE_cm3_mol_interp,Interpolation_Method);
                                        XB_uncor = Vol_BSSE_cm3_mol_interp;
                                        %YB_uncor = Energy_BSSE_interp;
                                    elseif X_Axis == "d"
                                        % Mass
                                        BMassMetal = elements('Symbol',BMetal,'atomic_mass'); % g/mol
                                        BMassHalide = elements('Symbol',BHalide,'atomic_mass'); % g/mol
                                        BMolar_Mass = BMassMetal + BMassHalide; % g/mol
                                        % Volume
                                        Vol_BSSE_cm3 = (Vol_BSSE.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(CurStructure)).*NA; % cm^3 / mol

                                        Density_BSSE = BMolar_Mass./Vol_BSSE_cm3; % g / cm^3

                                        % Generate finer grid
                                        Density_BSSE_interp = Density_BSSE(1):-InterpStepSize:Density_BSSE(end);
                                        %Energy_BSSE_interp = interp1(Density_BSSE,Energy_BSSE,Density_BSSE_interp,Interpolation_Method);
                                        XB_uncor = Density_BSSE_interp;
                                        %YB_uncor = Energy_BSSE_interp;
                                    end
                                % No interpolation
                                else
                                    if X_Axis == "a" || X_Axis == "P"
                                        XB_uncor = a_BSSE;
                                        %YB_uncor = Energy_BSSE;
                                    elseif X_Axis == "bl"  
                                        XB_uncor = BL_BSSE;
                                        %YB_uncor = Energy_BSSE;
                                    elseif X_Axis == "V"
                                        % Unit conversion
                                        Vol_BSSE_cm3_mol = (Vol_BSSE.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(CurStructure)).*NA; % cm^3 per mol
                                        XB_uncor = Vol_BSSE_cm3_mol;
                                        %YB_uncor = Energy_BSSE;
                                    elseif X_Axis == "d"
                                        % Mass
                                        BMassMetal = elements('Symbol',BMetal,'atomic_mass'); % g/mol
                                        BMassHalide = elements('Symbol',BHalide,'atomic_mass'); % g/mol
                                        BMolar_Mass = BMassMetal + BMassHalide; % g/mol

                                        % Volume
                                        Vol_BSSE_cm3 = (Vol_BSSE.*cm3_per_A3./Atom_Pairs_Per_Unit_Cell(CurStructure)).*NA; % cm^3 / mol

                                        Density_BSSE = BMolar_Mass./Vol_BSSE_cm3; % g / cm^3
                                        XB_uncor = Density_BSSE;
                                        %YB_uncor = Energy_BSSE;
                                    end
                                end
                                
                                % Realign data before plotting
                                if X_align && X_Axis ~= "P"
                                    idx_type = find(strcmp(struct_data(:,2),CurStructure));
                                    idx_method = find(strcmp(struct_data(:,1),DFTmethod_BSSE));
                                    idx_match = intersect(idx_type,idx_method);
                                    if isempty(idx_match)
                                        continue
                                    end
                                    
                                    XBmin = struct_data{idx_match,3};
                                    XB = XB_uncor - XBmin;
                                else
                                    XB = XB_uncor;
                                end
                                
                                % Plot Data
                                hold on
                                if strcmp(plot_style,'Scatter') || strcmp(plot_style,'Bar')
                                    pp{qq} = scatter(axes2,XB,Energy_BSSE,40,'k',...
                                        'Marker','o','Visible','on');
                                    pp{qq}.UserData = {CurSalt_BSSE DFTmethod_BSSE CurStructure Basis_BSSE};
                                elseif strcmp(plot_style,'Line')
                                    pp{qq} = plot(axes2,XB,Energy_BSSE,'Color','k',...
                                        'LineWidth',lw,'LineStyle','-',...
                                        'Visible','on');
                                    pp{qq}.UserData = {CurSalt_BSSE DFTmethod_BSSE CurStructure Basis_BSSE};
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    pp = SetPlotColours(pp,Ctype,Color_Scheme);
    %% Plot options
    set(axes2,'box','on','TickLabelInterpreter','latex');
    set(axes2,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    if X_Axis == "a"
        if X_align 
            xlabel(axes2,'Lattice Parameter $a - a_{min}$ (\AA)','fontsize',fs,'Interpreter','latex');
            xlim(axes2,[-3 7])
        else
            xlabel(axes2,'Lattice Parameter $a$ (\AA)','fontsize',fs,'Interpreter','latex');
            xlim(axes2,[2 9]);
        end
    elseif X_Axis == "P"
            xlabel(axes2,'Lattice Parameter $a$ (\AA)','fontsize',fs,'Interpreter','latex');
            xlim(axes2,[2 9]);
    elseif X_Axis == "bl"
        if X_align
            xlabel(axes2," Bond Length $BL - BL_{min}$ (\AA)",'fontsize',fs,'Interpreter','latex');
            xlim(axes2,[-2 4]);
        else
            xlabel(axes2," Bond Length (\AA)",'fontsize',fs,'Interpreter','latex');
            xlim(axes2,[0 6]);
        end
    elseif X_Axis == "V"
        if X_align
            xlabel(axes2,'Crystal Volume $V - V_{min}$ (cm$^{3}$ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');
            xlim(axes2,[-20 40]);
        else
            xlabel(axes2,'Crystal Volume (cm$^{3}$ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');
            xlim(axes2,[0 60]);
        end
    elseif X_Axis == "d"
        if X_align
            xlabel(axes2,'Density $d - d_{min}$ (g cm$^{-3}$)','fontsize',fs,'Interpreter','latex');
            %xlim([-20 40]);
        else
            xlabel(axes2,'Density (g cm$^{-3}$)','fontsize',fs,'Interpreter','latex');
            %xlim([0 60]);
        end
    end

    ylabel('Est. BSSE Correction Energy (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');

    legend2 = clickableLegend(axes2,sub_type_BSSE(~cellfun('isempty',sub_type_BSSE)),...
        'Interpreter','latex','Parent',Mainfigh);

    set(axes2,'YLimMode','auto')

    axis(axes2,'manual');
    
    if strcmp(plot_style,'Scatter') || strcmp(plot_style,'Line')
        set(legend2,'Position',[0.842708333333333 0 0.156770833333333 0.496365524402908]);
    elseif strcmp(plot_style,'Bar')
    	set(legend2,'Position',[0.847 0.00207684319833845 0.152 0.45]); % [ X Y Width Height]
    end
    
    
    
    % Initialize with all plots hidden if selected
    if ~Auto_Plot
        for ii=1:length(pp)
            if ~isempty(pp{ii})
                set(pp{ii},'UserData', false);
                set(pp{ii},'Visible','off');
                set(findobj('Tag',pp{ii}.DisplayName),'Visible','off');
            end
        end
    end
    set(axes2,'Visible','on')
end

end