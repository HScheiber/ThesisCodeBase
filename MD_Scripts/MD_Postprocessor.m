% Script to run any post-processing to MD simulations
function MD_Postprocessor(Settings)

if ~isstruct(Settings)
    Settings = load(Settings,'-mat').Settings;
end

py.LiXStructureDetector.Check_Structures(Settings.WorkDir, Settings.Salt,...
            pyargs('SystemName',Settings.JobName,...
            'FileType',Settings.CoordType,...
            'SavePredictionsImage',Settings.SavePredictionsImage,...
            'ML_TimeLength',Settings.ML_TimeLength,...
            'ML_TimeStep',Settings.ML_TimeStep,...
            'TimePerFrame',Settings.TimePerFrame,...
            'SaveTrajectory',Settings.SaveTrajectory,...
            'SaveFeatures', Settings.SaveFeatures,...
            'SavePredictions',Settings.SavePredictions,...
            'Qlm_Average',Settings.Qlm_Average,...
            'Voronoi',Settings.Voronoi));

if ~isempty(Settings.RunEnergyAnalysis)
    Settings = Update_MD_Settings(Settings);
    Energy_file = fullfile(Settings.WorkDir,[Settings.JobName '.edr']);
    
    % Check energy options
    gmx_command = [Settings.wsl 'echo "0" ' Settings.pipe ...
        ' ' strrep(Settings.gmx_loc,Settings.wsl,'') Settings.g_energy ...
        ' -f ' windows2unix(Energy_file)];
    [~,outpt] = system(gmx_command);
    
    en_opts = regexp(outpt,'End your selection with an empty line or a zero.\n-+(.+?)\n\n','tokens','once');
    en_opts = en_opts{1};
    
    En_set = '';
    for idx = 1:numel(Settings.RunEnergyAnalysis)
        En_set = [En_set char(regexpi(en_opts,['([0-9]{1,2})  ' Settings.RunEnergyAnalysis{idx}],'tokens','once')) ' ']; %#ok<AGROW>
    end
    En_set = regexprep([En_set '0'],' +',' ');
    
    % Convert the energy file into readable form
    Energy_output = fullfile(Settings.WorkDir,[Settings.JobName '_energies.xvg']);
    
    eneconv_cmd = [Settings.wsl 'echo ' En_set ' ' Settings.pipe ...
        ' ' strrep(Settings.gmx_loc,Settings.wsl,'') Settings.g_energy ...
        ' -f ' windows2unix(Energy_file) ' -o ' windows2unix(Energy_output)];
    
    % Run it
    [state,eneconv_output] = system(eneconv_cmd);
    
    if state ~= 0
        disp(eneconv_cmd)
        disp(eneconv_output)
        error('Unable to collect energy.')
    end
    
    % Get total atoms filename
    grofile = fullfile(Settings.WorkDir,[Settings.JobName '.' Settings.CoordType]);
    Data = load_gro_file(grofile);

    % Open the mat file to get total number of atoms/molecules
    if Settings.Polarization
        Energy.NF = Data.N_atoms/4;
    else
        Energy.NF = Data.N_atoms/2;
    end
    
    % Import energy file as structure
    Energy_imp = import_xvg_as_struc(Energy_output);
    Energy = cell2struct([struct2cell(Energy_imp);struct2cell(Energy)],[fieldnames(Energy_imp);fieldnames(Energy)]);
    
    Output_file = fullfile(Settings.WorkDir,[Settings.JobName '_Energy.mat']);
    save(Output_file,'Energy');
end
        
fclose(fopen(fullfile(Settings.WorkDir,'POSTPROCESS_COMPLETE'), 'w'));
end