% CP2K_deriv_free_optimizer
%% System and calculation information
Salt = '##SALT##';
Structure = '##STRUCTURE##';
Theory = '##THEORY##';
cdir = '##CDIR##';
MPIcores = ##MPICORES##;
OMPcores = ##OMPCORES##;
Parallel_search = false;
display_option = 'final';
MaxCycles = 100;
Cutoff = ##CUTOFF##;
vDW_Cutoff = ##VDWCUTOFF##;
Finer_Grid = ##FINERGRID##;
Remove_Temp = false;
kpoints = ##KPOINTS##;
Basis_Experiment = '##BASIS_EXPERIMENT##';
Eps_pgf_orb = ##EPS_PGF_ORB##;
Eps_default = ##EPS_DEFAULT##;
Eps_SCF = ##EPS_SCF##;
restart = true; % Restart from lowest energy state if available

%% Grab the initial geometry
home = find_home;
Data_Directory = [home filesep 'data'] ;
CP2K_Data_Obj = load(fullfile(Data_Directory,'CP2K_Data.mat'));

switch Salt
    case {'LiF' 'LiCl' 'NaCl'}
        Ref_Theory = 'TMTPSS_rVV10L';
        Ref_Basis_Experiment = 'Sapporo_QZP';
    case {'LiAt'}
        Ref_Theory = 'PBE_D3_DKH';
        Ref_Basis_Experiment = 'Sapporo_DZP';
    otherwise
        Ref_Theory = 'TMTPSS_rVV10L_DKH';
        Ref_Basis_Experiment = 'Sapporo_QZP';
end

if restart
    is_complete = cp2k_check_complete(cdir,true);
    if is_complete
        disp([Job_name ': Calculation completed and restart is ON - haulting calculation.'])
        return
    end
    
    d_inp = dir([cdir filesep '*.inp']);
    if isempty(d_inp)
        a = CP2K_Data_Obj.Data.(Ref_Basis_Experiment).(Ref_Theory).(Salt).(Structure).a;
        b = CP2K_Data_Obj.Data.(Ref_Basis_Experiment).(Ref_Theory).(Salt).(Structure).b;
        c = CP2K_Data_Obj.Data.(Ref_Basis_Experiment).(Ref_Theory).(Salt).(Structure).c;
        conv_to_prim = true;
    else
        inptxt = fileread(fullfile(cdir,d_inp(1).name));
        comptxt = regexp(inptxt,'ABC +(-|\.|[0-9]|E)+ +(-|\.|[0-9]|E)+ +(-|\.|[0-9]|E)+','tokens','once');
        LatPar.a = str2double(comptxt{1});
        LatPar.b = str2double(comptxt{2});
        LatPar.c = str2double(comptxt{3});
        conv_to_prim = false;
        a = LatPar.a;
    end
else
    a = CP2K_Data_Obj.Data.(Ref_Basis_Experiment).(Ref_Theory).(Salt).(Structure).a;
    b = CP2K_Data_Obj.Data.(Ref_Basis_Experiment).(Ref_Theory).(Salt).(Structure).b;
    c = CP2K_Data_Obj.Data.(Ref_Basis_Experiment).(Ref_Theory).(Salt).(Structure).c;
    conv_to_prim = true;
end

if isnan(a)
    Settings = Initialize_MD_Settings;
    Settings.Salt = Salt;
    Settings.Structure = Structure;
    Settings.Use_Conv_cell = false;
    Default_Geometry = Default_Crystal(Settings);
    LatPar.a = Default_Geometry.a;
    LatPar.b = Default_Geometry.b;
    LatPar.c = Default_Geometry.c;
    conv_to_prim = false;
end

% Convert to primitive cell
if conv_to_prim
    switch lower(Structure)
        case 'rocksalt'
            LatPar.a = a/sqrt(2);
            LatPar.b = b/sqrt(2);
            LatPar.c = c/sqrt(2);
            x0 = LatPar.a;
        case 'wurtzite'
            LatPar.a = a;
            LatPar.b = b;
            LatPar.c = c;
            x0 = [LatPar.a LatPar.c];
        case 'fivefive'
            LatPar.a = c;
            LatPar.b = (sqrt(3)/3)*b;
            LatPar.c = a;
            x0 = [LatPar.a LatPar.c];
        case 'cscl'
            LatPar.a = a;
            LatPar.b = b;
            LatPar.c = c;
            x0 = LatPar.a;
        case 'betabeo'
            LatPar.a = a;
            LatPar.b = b;
            LatPar.c = c;
            x0 = [LatPar.a LatPar.b LatPar.c];
        case 'sphalerite'
            LatPar.a = a/sqrt(2);
            LatPar.b = b/sqrt(2);
            LatPar.c = c/sqrt(2);
            x0 = LatPar.a;
        case {'nias' 'antinias'}
            LatPar.a = a;
            LatPar.b = b;
            LatPar.c = c;
            x0 = [LatPar.a LatPar.c];
        otherwise
            error(['Unknown structure: ' Structure])
    end
else
    switch lower(Structure)
        case {'rocksalt' 'cscl' 'sphalerite'}
            x0 = LatPar.a;
        case {'wurtzite' 'fivefive' 'nias' 'antinias' 'betabeo'}
            x0 = [LatPar.a LatPar.c];
        otherwise
            error(['Unknown structure: ' Structure])
    end
end

[Metal,Halide] = Separate_Metal_Halide(Salt);

switch Basis_Experiment
    case 'pob-TZVP'
        Metal_Basis = 'pob-TZVP';
        switch Halide
            case {'F' 'Cl'}
                Halide_Basis = 'pob-TZVP';
            otherwise
                Halide_Basis = 'Sapporo-DKH3-TZP';
        end
    case 'Sapporo-TZP'
        Metal_Basis = 'pob-TZVP'; % 'Sapporo-DKH3-TZP' 'Sapporo-TZP' 'pob-TZVP'
        switch Halide
            case {'F' 'Cl'}
                Halide_Basis = 'Sapporo-TZP';
            otherwise
                Halide_Basis = 'Sapporo-DKH3-TZP';
        end
    case 'Sapporo-QZP'
        Metal_Basis = 'pob-TZVP';
        switch Halide
            case {'F' 'Cl'}
                Halide_Basis = 'Sapporo-QZP';
            otherwise
                Halide_Basis = 'Sapporo-DKH3-QZP';
        end
    case 'Sapporo-DZP'
        Metal_Basis = 'pob-TZVP';
        Halide_Basis = 'Sapporo-DKH3-DZP-2012';
end

setenv('OMP_NUM_THREADS',num2str(OMPcores))

%% Define the loss function
fun = @(x)Calculate_CP2K_Energy(x,Salt,Structure,Theory,MPIcores,...
    'Cutoff',Cutoff,'vdw_cutoff',vDW_Cutoff,'Finer_Grid',Finer_Grid,...
    'kpoints',kpoints,'Remove_Temp',Remove_Temp,'Halide_Basis',Halide_Basis,...
    'Metal_Basis',Metal_Basis,'Eps_pgf_orb',Eps_pgf_orb,'Eps_default',Eps_default,...
    'Eps_SCF',Eps_SCF);

%% Define an output display function
switch Salt
    case {'LiF' 'LiCl' 'NaCl'}
        Theory_str = strrep(Theory,'-','_');
    otherwise
        Theory_str = [strrep(Theory,'-','_') '_DKH'];
end
Basis_Experiment_str = strrep(Basis_Experiment,'-','_');

try
    M_En = CP2K_Data_Obj.Data.(Basis_Experiment_str).(Theory_str).(Metal).Energy;
    X_En = CP2K_Data_Obj.Data.(Basis_Experiment_str).(Theory_str).(Halide).Energy;
catch
    M_En = nan;
    X_En = nan;
end
clear('CP2K_Data_Obj');

display_function = @(optimvalues,flag)check_LE(optimvalues,flag,Structure,M_En,X_En);

% Assign optimizer parameters
options = optimoptions(@patternsearch,'Display',display_option,'MaxIterations',MaxCycles,...
    'UseParallel',Parallel_search,'UseVectorized',false,'PlotFcn',display_function,...
    'InitialMeshSize',0.05,'StepTolerance',1e-5,'FunctionTolerance',1e-3,...
    'PollOrderAlgorithm','Success','Cache','off','UseCompletePoll',false,...
    'PollMethod','GPSPositiveBasis2N','MaxMeshSize',0.1,'MeshContractionFactor',0.1,...
    'AccelerateMesh',false,'MeshTolerance',1e-8);

disp(['Beginning ' Salt ' ' Structure ' ' Theory '/' Metal_Basis '/' Halide_Basis ' Optimization...'])
[lattice_params,E_opt_param] = patternsearch(fun,x0,[],[],[],[],[],[],[],options);


% disp(['Beginning ' Salt ' ' Structure ' ' Theory '/' Metal_Basis '/' Halide_Basis ' Optimization...'])
% optionsNM = optimset('Display',display_option,'MaxIter',MaxCycles,...
%     'TolFun',1e-6,'TolX',1e-5,'MaxFunEvals',Inf);
% 
% [lattice_params,E_opt_param] = fminsearch(fun,x0,optionsNM);

outfile = fullfile(cdir,'CALCULATION_CONVERGED.mat');
save(outfile,'lattice_params','E_opt_param');
