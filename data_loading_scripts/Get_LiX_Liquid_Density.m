% Returns density in units of [molecules]/[nm^3]
function density = Get_LiX_Liquid_Density(Settings)
    warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')
    home = find_home;
    datdir = fullfile(home,'data','ThermExp');
    Exp_Model = [Settings.Salt '_L_Exp.mat']; % Use this model if no matching model is available
    
    dirs = dir(datdir);
    models = {};
	for i=1:length(dirs)
        model = dirs(i).name;
        if ~strncmp(model,'.',1) && ~strcmp(model,Exp_Model) && contains(model,[Settings.Salt '_L'])
            models{end+1} = model;
        end
	end
    
    % If multiple results: Remove any that are the incorrect theory
    if length(models) > 1
        for i=length(models):-1:1
            model = models{i};
            if ~contains(model,Settings.Theory)
                models(i) = [];
            end
        end
    end
    
    % If multiple results: remove any that are not the correct model
    if length(models) > 1
        for i=length(models):-1:1
            model = models{i};
            if ~contains(model,Settings.Model)
                models(i) = [];
            end
        end
    end
    
    if isempty(models)
        model = Exp_Model;
        Scale_density = 0.95;
    else
        model = models{1};
        Scale_density = 1;
    end
    
    try
        dat = load(fullfile(datdir,model));
    catch
        dat = load(fullfile(datdir,Exp_Model));
    end
    
    Volume = dat.mdl.predict(Settings.Target_T); % Units of: [Angstrom^3]/[formula unit]
    
    % Convert to [molecules]/[nm^3]
    density = Scale_density/(Volume.*(0.001));
    warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle')
end