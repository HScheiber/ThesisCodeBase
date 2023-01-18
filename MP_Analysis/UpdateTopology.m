function UpdateTopology(Settings)

[Settings.Metal,Settings.Halide] = Separate_Metal_Halide(Settings.Salt);
Metal_Info = elements('Sym',Settings.Metal);
Halide_Info = elements('Sym',Settings.Halide);

% Load Model parameters
if isempty(Settings.Model)
    % Default setting of model parameters if a particular model is not selected
    Settings.S = Init_Scaling_Object;
    Settings.CR_Damp = Init_CRDamping_Object; % note: b = steepness of damping, r_d = position in nm. liCl = 0.21, LiI = 0.26
    Settings.C6_Damp = Init_C6Damping_Object;
    [Settings.GAdjust_MX,Settings.GAdjust_MM,Settings.GAdjust_XX] = Init_GAdjust_Object;
else
    [Settings,~] = Load_Model_Params(Settings);
end

% Check if parameters can be input without a table;
Settings.Table_Req = IsGmxTableRequired(Settings);

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

    % Define the repulsive parameters as 1.0 (already included in potentials)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETA##','1.0');
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALA##','1.0');
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALA##','1.0');
    
    % Dispersion coefficients
    [~,C6] = MakeTables(Settings,'SaveTables',false);
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETC##',num2str(C6.MM,'%.10e'));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALC##',num2str(C6.XX,'%.10e'));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALC##',num2str(C6.MX,'%.10e'));    
elseif strcmp(Settings.Theory,'LJ')
    
    % Definte the function type as 1 (LJ)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##NBFUNC##','1');

    % Define the combination rules (Lorenz-berthelot in sigma-epsilon form)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##COMBR##','2');
    
    % Add parameters to topology text
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETC##',pad(num2str(Settings.S.S.MM,'%10.8e'),10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALC##',pad(num2str(Settings.S.S.XX,'%10.8e'),10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALC##',pad(num2str(Settings.S.S.MX,'%10.8e'),10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETA##',num2str(Settings.S.E.MM,'%10.8e'));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALA##',num2str(Settings.S.E.XX,'%10.8e'));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALA##',num2str(Settings.S.E.MX,'%10.8e'));
    
elseif contains(Settings.Theory,'JC')
    
    switch Settings.Theory
        case 'JC'
            Settings.WaterModel = 'SPC/E';
        case 'JC3P'
            Settings.WaterModel = 'TIP3P';
        case 'JC4P'
            Settings.WaterModel = 'TIP4PEW';
    end

    % Definte the function type as 1 (LJ)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##NBFUNC##','1');

    % Define the combination rules (Lorenz-berthelot in sigma-epsilon form)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##COMBR##','2');

    % Get JC parameters
    [U_MX,U_MM,U_XX] = JC_Potential_Parameters(Settings);
    
    % Add parameters to topology text
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETC##',pad(num2str(U_MM.sigma,'%10.8e'),10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALC##',pad(num2str(U_XX.sigma,'%10.8e'),10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALC##',pad(num2str(U_MX.sigma,'%10.8e'),10));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETA##',num2str(U_MM.epsilon,'%10.8e'));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALA##',num2str(U_XX.epsilon,'%10.8e'));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALA##',num2str(U_MX.epsilon,'%10.8e'));
else
    rmdir(Settings.WorkDir)
    error(['Warning: Unknown theory type: "' Settings.Theory '".'])
end

Settings.Topology_Text = strrep(Settings.Topology_Text,'##N##x##N##x##N##','');
Settings.Topology_Text = strrep(Settings.Topology_Text,'##GEOM##','Liquid-Crystal Interface');

Settings.Topology_File = fullfile(Settings.WorkDir,[Settings.JobName '.top']);
Atomlist = copy_atom_order(Settings.SuperCellFile);
Settings.Topology_Text = strrep(Settings.Topology_Text,'##LATOMS##',Atomlist);
fidTOP = fopen(Settings.Topology_File,'wt');
fwrite(fidTOP,regexprep(Settings.Topology_Text,'\r',''));
fclose(fidTOP);
    
end