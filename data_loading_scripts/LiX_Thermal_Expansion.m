% Returns expansion scale
function Exp_Coeff = LiX_Thermal_Expansion(Settings)
    
    home = find_home;
    datdir = fullfile(home,'data','ThermExp');
    if contained_in_cell(Settings.Salt,{'CsCl' 'CsBr' 'CsI'})
        Exp_Model = [Settings.Salt '_C_Exp.mat']; % Use this model if no matching model is available
    else
        Exp_Model = [Settings.Salt '_R_Exp.mat']; % Use this model if no matching model is available
    end
    
    dirs = dir(datdir);
    models = {};
	for i=1:length(dirs)
        model = dirs(i).name;
        if ~strncmp(model,'.',1) && contains(model,[Settings.Salt '_' Settings.Geometry.Label])
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
    else
        model = models{1};
    end
    
    dat = load(fullfile(datdir,model));
    Exp_Coeff = dat.mdl.predict(Settings.Target_T)/dat.mdl.predict(0);
end