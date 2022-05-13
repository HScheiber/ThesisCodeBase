% {'PBE' 'optB88' 'optB88-vdW' 'PBE-D3' 'PBE-rVV10' 'PBE-rVV10L' 'PBEsol' 'PBEsol-rVV10' 'optPBE' 'optPBE-vdW' ...
%         'B86r' 'rev-vdW-DF2' 'PW86r' 'vdW-DF2' 'BEEF-vdW' ...
%         'SG4' 'SG4-rVV10m' 'SCAN' 'SCAN-rVV10' 'TMTPSS' 'TMTPSS-rVV10L'};
% {'optB88-vdW' 'PBE' 'PBE-D3' 'PBE-rVV10' 'PBE-rVV10L' 'PBEsol-rVV10' ...
%         'rev-vdW-DF2' 'SG4-rVV10m' 'SCAN' 'SCAN-rVV10' 'TMTPSS-rVV10L'};

function [Energy,cdir_temp] = Calculate_CP2K_Energy(Params,Salt,Structure,Theory,MPIcores,varargin)

% Parse optional inputs
p = inputParser;
p.FunctionName = 'Calculate_CP2K_Energy';
addOptional(p,'Cutoff',600,@(x)validateattributes(x,...
    {'numeric'},{'nonempty'}))
addOptional(p,'kpoints',6,@(x)validateattributes(x,...
    {'numeric'},{'nonempty'}))
addOptional(p,'Eps_default',1e-10,@(x)validateattributes(x,...
    {'numeric'},{'nonempty'}))
addOptional(p,'vdw_cutoff',600,@(x)validateattributes(x,...
    {'numeric'},{'nonempty'}))
addOptional(p,'Eps_SCF',1e-6,@(x)validateattributes(x,...
    {'numeric'},{'nonempty'}))
addOptional(p,'Eps_pgf_orb',1e-5,@(x)validateattributes(x,...
    {'numeric'},{'nonempty'}))
addOptional(p,'Metal_Basis','pob-TZVP',@(x)validateattributes(x,...
    {'char'},{'nonempty'}))
addOptional(p,'Smoothing',true,@(x)validateattributes(x,...
    {'logical'},{'nonempty'}))
addOptional(p,'Finer_Grid',false,@(x)validateattributes(x,...
    {'logical'},{'nonempty'}))
addOptional(p,'Increase_Radial_Grid',0,@(x)validateattributes(x,...
    {'numeric'},{'nonempty'}))
addOptional(p,'Broden_SCF',false,@(x)validateattributes(x,...
    {'logical'},{'nonempty'}))

addOptional(p,'hfx_eps_schwarz',1e-9,@(x)validateattributes(x,...
    {'numeric'},{'nonempty'}))
addOptional(p,'hfx_max_memory',3600,@(x)validateattributes(x,...
    {'numeric'},{'nonempty'}))
addOptional(p,'hfx_cutoff_radius',5,@(x)validateattributes(x,...
    {'numeric'},{'nonempty'}))

switch Salt
    case {'LiBr' 'LiI' 'LiAt'}
        addOptional(p,'Halide_Basis','Sapporo-DKH3-QZP',@(x)validateattributes(x,...
            {'char'},{'nonempty'}))
        addOptional(p,'Relativistic',true,@(x)validateattributes(x,...
            {'logical'},{'nonempty'}))
    otherwise
        addOptional(p,'Halide_Basis','Sapporo-QZP',@(x)validateattributes(x,...
            {'char'},{'nonempty'}))
        addOptional(p,'Relativistic',false,@(x)validateattributes(x,...
            {'logical'},{'nonempty'}))
end

addOptional(p,'Basis_Experiment','Sapporo-QZP',@(x)validateattributes(x,...
    {'char'},{'nonempty'}))
addOptional(p,'Remove_Temp',false,@(x)validateattributes(x,...
    {'logical'},{'nonempty'}))

parse(p,varargin{:});

%% Calculation parameters
Cutoff = p.Results.Cutoff;
kpoints = p.Results.kpoints; % does not apply to hybrids
Eps_default = p.Results.Eps_default;
Eps_pgf_orb = p.Results.Eps_pgf_orb;
Eps_SCF = p.Results.Eps_SCF;
vdw_cutoff = p.Results.vdw_cutoff; % Only applies when using vdW or rVV10 corrections
Metal_Basis = p.Results.Metal_Basis;
Halide_Basis = p.Results.Halide_Basis;
Basis_Experiment = p.Results.Basis_Experiment;
Relativistic = p.Results.Relativistic;
smoothing = p.Results.Smoothing; % uses SPLINE3 for XC derivative method
finer_grid = p.Results.Finer_Grid; % Use finer grid for the XC
Broden_SCF = p.Results.Broden_SCF;
Increase_Radial_Grid = p.Results.Increase_Radial_Grid; % Increase radial grid size by this amount
Remove_Temp = p.Results.Remove_Temp; % Removes temporary calculation folders

% HFX options (hybrids only)
hfx.eps_schwarz = p.Results.hfx_eps_schwarz;
hfx.max_memory = p.Results.hfx_max_memory;
hfx.cutoff_radius = p.Results.hfx_cutoff_radius;

% Relativistic options
Rel.Method = 'DKH';
Rel.DKH_Order = 3;
Rel.Potential = 'FULL';
Rel.Transformation = 'ATOM';
Rel.ZORA_Type = 'MP'; % only applies to ZORA calculations
Rel.Z_Cutoff = 1;

%% Grab input templates
if ispc
    template_loc = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\CP2K_scripts';
    ALL_Potentials_loc = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\CP2K_scripts\ALL_POTENTIALS';
    pdir = pwd;
    cp2k_exe = 'cp2k';
else
    [~,Server] = system('hostname');
    if ~isempty(regexp(Server(1:3),'se[0-9]','ONCE')) || strcmpi(Server(1:3),'log')
        template_loc = '/home/haydensc/ThesisCodeBase/CP2K_scripts';
        ALL_Potentials_loc = '/home/haydensc/ThesisCodeBase/CP2K_scripts/ALL_POTENTIALS';
        pdir = ['/home/haydensc/scratch/CP2K_Calcs/' Basis_Experiment];
        if ~exist(pdir, 'dir')
           mkdir(pdir)
        end
        cd(pdir)
        cp2k_exe = 'cp2k';
        setenv('CP2K_DATA_DIR','/home/haydensc/cp2k-master/data')
    elseif strncmpi(Server,'pat',3)
        template_loc = '/home/user/ThesisCodeBase/CP2K_scripts';
        ALL_Potentials_loc = '/home/user/ThesisCodeBase/CP2K_scripts/ALL_POTENTIALS';
        pdir = ['/home/user/project/CP2K_Calcs/' Basis_Experiment];
        if ~exist(pdir, 'dir')
           mkdir(pdir)
        end
        cd(pdir)
        cp2k_exe = 'cp2k';
        setenv('CP2K_DATA_DIR','/home/user/cp2k/data')
    else
        template_loc = '/home/scheiber/ThesisCodeBase/CP2K_scripts';
        ALL_Potentials_loc = '/home/scheiber/ThesisCodeBase/CP2K_scripts/ALL_POTENTIALS';
        pdir = ['/project/6001647/scheiber/CP2K_Calcs/' Basis_Experiment];
        if ~exist(pdir, 'dir')
           mkdir(pdir)
        end
        cd(pdir)
        cp2k_exe = ['mpiexec -np ' num2str(MPIcores) ' cp2k.popt'];
    end
end

cp2k_input_template = fileread([template_loc filesep 'input_template.inp']);

% Add in some universal parameters
cp2k_input_template = strrep(cp2k_input_template,'##CUTOFF##',num2str(Cutoff));
cp2k_input_template = strrep(cp2k_input_template,'##EPSDEFAULT##',num2str(Eps_default));
cp2k_input_template = strrep(cp2k_input_template,'##EPSPGFORB##',num2str(Eps_pgf_orb));
cp2k_input_template = strrep(cp2k_input_template,'##EPSSCF##',num2str(Eps_SCF));
cp2k_input_template = strrep(cp2k_input_template,'##ALL_POTENTIAL_FILE##',ALL_Potentials_loc);
cp2k_input_template = strrep(cp2k_input_template,'##RUNTYPE##','ENERGY');

% Relativistic parameters
if Relativistic
    Reltxt = cp2k_relativistic(Rel); %#ok<*UNRCH>
    cp2k_input_template = strrep(cp2k_input_template,'##RELATIVISTIC##',Reltxt);
    func_add = ['-' Rel.Method];
else
    cp2k_input_template = regexprep(cp2k_input_template,'\s*##RELATIVISTIC##','','once');
    func_add = '';
end

if Broden_SCF
    Broden = ['		&MIXING' newline ...
            '			ALPHA 0.5' newline ...
            '			METHOD BROYDEN_MIXING' newline ...
            '			NBUFFER 4' newline ...
            '			BROY_W0 0.0001' newline ...
            '		&END MIXING'];
    cp2k_input_template = strrep(cp2k_input_template,'##BRODEN##',Broden);
else
    cp2k_input_template = regexprep(cp2k_input_template,'\s*##BRODEN##','','once');
end

[XC_Func_input,Hybrid,~] = cp2k_XC_Func(Theory,vdw_cutoff,hfx,smoothing,finer_grid);

cp2k_input_template = strrep(cp2k_input_template,'##XCFUNC##',XC_Func_input);

% Remove MOTION section, not required for ENERGY calculations
cp2k_input_template = regexprep(cp2k_input_template,'\s*&MOTION.+&END MOTION','','once');
cp2k_input_template = regexprep(cp2k_input_template,'\s*STRESS_TENSOR ##STRESSTENSOR##','','once');

%% Deal with basis sets
[Metal,Halide] = Separate_Metal_Halide(Salt);
Basis_Set_M = cp2k_basis_set(Metal_Basis,Metal,Increase_Radial_Grid);
Basis_Set_X = cp2k_basis_set(Halide_Basis,Halide,Increase_Radial_Grid);
Basis_Set = [Basis_Set_M newline Basis_Set_X];
cp2k_input_template = strrep(cp2k_input_template,'##BASISSET##',Basis_Set);

%% Deal with crystal geometry
[Structure_Coord,Structure_Angle,ABC,Structure_Symmetry] = ...
    cp2k_Structures_fixed(Structure,Metal,Halide,Params);
[kpoints_txt,nrep] = cp2k_kpoints(kpoints,Hybrid,Structure,ABC,hfx.cutoff_radius);

cp2k_input_template = strrep(cp2k_input_template,'##ANGLES##',Structure_Angle);
cp2k_input_template = strrep(cp2k_input_template,'##COORDS##',Structure_Coord);
cp2k_input_template = strrep(cp2k_input_template,'##SYMS##',Structure_Symmetry);
cp2k_input_template = strrep(cp2k_input_template,'##ABCPAR##',ABC);
if isempty(kpoints_txt)
    cp2k_input_template = regexprep(cp2k_input_template,'\s*##KPOINTS##','','once');
else
    cp2k_input_template = strrep(cp2k_input_template,'##KPOINTS##',kpoints_txt);
end
cp2k_input_template = strrep(cp2k_input_template,'##REPA##',nrep.a);
cp2k_input_template = strrep(cp2k_input_template,'##REPB##',nrep.b);
cp2k_input_template = strrep(cp2k_input_template,'##REPC##',nrep.c);

cdir_temp = tempname([pdir filesep Theory func_add filesep Structure filesep Salt]);
cdir = [pdir filesep Theory func_add filesep Structure filesep Salt];
if ~exist(cdir_temp, 'dir')
   mkdir(cdir_temp)
end
Job_name = [Theory func_add '_' Structure '_' Salt];
Outfile = [Job_name '.out'];

% Check for wavefunction restart file
d_wf = dir([cdir filesep '*RESTART.kp']);

if ~isempty(d_wf)
    % If multiple restart files exist, choose most recent
    if size(d_wf,1) > 1
        Times = zeros(1,size(d_wf,1));
        for y = 1:size(d_wf,1)
            Times(y) = datenum(d_wf(y).date,...
                'dd-mmm-yyyy HH:MM:SS');
        end
        [~,Idx] = max(Times);
        d_wf = d_wf(Idx);
    end

    restxt = fullfile(d_wf.folder,d_wf.name);

    cp2k_input_template = regexprep(cp2k_input_template,'POTENTIAL_FILE_NAME (.+?)\n',...
        ['POTENTIAL_FILE_NAME $1' newline '    WFN_RESTART_FILE_NAME ' restxt newline],'once');
end

cp2k_input_template = strrep(cp2k_input_template,'##JOBNAME##',Job_name);

inp_text_filename = fullfile(cdir_temp,[Job_name '.inp']);
fid = fopen(inp_text_filename,'w');
fprintf(fid,'%s',cp2k_input_template);
fclose(fid);

cd(cdir_temp)
% Run the cp2k calculation
system([cp2k_exe ' ' Job_name '.inp' ' > ' Outfile]);

% Grab the energy
if isempty(Outfile)
    disp('Note: No cp2k output file found!');
    Energy = nan;
    if Remove_Temp
        rmdir(cdir_temp, 's')
    end
    return
end

jobtxt = fileread(Outfile);
comptxt = regexp(jobtxt,'SCF run converged in.+DBCSR STATISTICS','match','once');
if isempty(comptxt)
    disp('Note: Cp2k calculation did not Converge!');
end

Entxt = regexp(comptxt,'ENERGY\| Total FORCE_EVAL \( QS \) energy .a\.u\..: +(-|\.|[0-9]|E)+','tokens','once');
if isempty(Entxt)
    disp('Note: Unable to find total energy output!');
    Energy = nan;
    if Remove_Temp
        rmdir(cdir_temp, 's')
    end
    return
end

switch Structure
    case {'Rocksalt' 'CsCl' 'Sphalerite'}
        NF_per_Cell = 1;
    case {'NiAs' 'Wurtzite' 'FiveFive' 'AntiNiAs'}
        NF_per_Cell = 2;
    case 'BetaBeO'
        NF_per_Cell = 4;
end
num_cells = str2double(nrep.a)*str2double(nrep.b)*str2double(nrep.c);
Energy = str2double(Entxt{1})/(NF_per_Cell*num_cells);

%% Save the final energy to file
e_filename = fullfile(cdir,'lowest_energy.mat');

% Check for lowest-energy file
if isfile(e_filename)
    % If it exists, load its energy
    Prev = load(e_filename);
else
    % Otherwise set the energy to Infinity
    Prev.Energy = Inf;
end

% If a new lowest energy state is found, copy over all files to the outer
% calculation directory
if Prev.Energy > Energy
    if isfile(e_filename)
        delete(e_filename)
    end
    save(e_filename,'Energy')
    if Remove_Temp
        movefile([cdir_temp filesep '*'], cdir)
    else
        copyfile([cdir_temp filesep '*'], cdir)
    end
end

% Delete the temporary calculation directory when finished
if Remove_Temp
    rmdir(cdir_temp, 's')
end

end