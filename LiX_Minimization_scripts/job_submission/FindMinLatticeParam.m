function [Settings,Updated] = FindMinLatticeParam(Settings,varargin)

% Optional input settings
p = inputParser;
p.FunctionName = 'FindMinLatticeParam';
addOptional(p,'Find_Similar_Params',true,@(x)validateattributes(x,{'logical'},{'nonempty'}))
addOptional(p,'Center_Coordinates',false,@(x)validateattributes(x,{'logical'},{'nonempty'}))
parse(p,varargin{:});

Find_Similar_Params = p.Results.Find_Similar_Params; % Gathers components of energy at the end in addition to the total energy.
Center_Coordinates = p.Results.Center_Coordinates;

D = load(fullfile(Settings.home,'data','MX_JCTF_Min_Data.mat'),'Data');
Data = D.Data;

if ~isempty(Settings.Model)
    FullModelName = [Settings.Theory '_Model_' Settings.Model];
else
    FullModelName = Settings.Theory;
end

if isfield(Data,Settings.Salt) && isfield(Data.(Settings.Salt),FullModelName) && ...
    isfield(Data.(Settings.Salt).(FullModelName),Settings.Structure) && ...
    isfield(Data.(Settings.Salt).(FullModelName).(Settings.Structure),'a')
    
    a = Data.(Settings.Salt).(FullModelName).(Settings.Structure).a;
    b = Data.(Settings.Salt).(FullModelName).(Settings.Structure).b;
    c = Data.(Settings.Salt).(FullModelName).(Settings.Structure).c;
    Updated = true;
elseif Find_Similar_Params && isfield(Data,Settings.Salt) && isfield(Data.(Settings.Salt),Settings.Theory) && ...
    isfield(Data.(Settings.Salt).(Settings.Theory),Settings.Structure) && ...
    isfield(Data.(Settings.Salt).(Settings.Theory).(Settings.Structure),'a')
    
    a = Data.(Settings.Salt).(Settings.Theory).(Settings.Structure).a;
    b = Data.(Settings.Salt).(Settings.Theory).(Settings.Structure).b;
    c = Data.(Settings.Salt).(Settings.Theory).(Settings.Structure).c;
    Updated = false;
else
    Updated = false;
    return
end

% Place into primitive cell if needed!
switch lower(Settings.Structure)
    case 'rocksalt'
        if Settings.Use_Conv_cell
            Settings.Geometry.a = a;
            Settings.Geometry.b = b;
            Settings.Geometry.c = c;
        else
            Settings.Geometry.a = a/sqrt(2);
            Settings.Geometry.b = b/sqrt(2);
            Settings.Geometry.c = c/sqrt(2);
            Settings.Geometry.FC_Metal  = [0.0 0.0 0.0];
            Settings.Geometry.FC_Halide = [1/2 1/2 1/2];
            if Center_Coordinates
                Settings.Geometry = CenterCoordinates(Settings.Geometry,Settings.Structure);
            end
            Settings.Geometry.FC = [Settings.Geometry.FC_Metal; Settings.Geometry.FC_Halide];
        end
    case 'wurtzite'
        Settings.Geometry.a = a;
        Settings.Geometry.b = b;
        Settings.Geometry.c = c;
    case 'fivefive'
        if Settings.Use_Conv_cell
            Settings.Geometry.a = a;
            Settings.Geometry.b = b;
            Settings.Geometry.c = c;
        else
            Settings.Geometry.a = (a/2)*1/(sqrt((1/4) + (1/9)*(sind(60)^2)));
            Settings.Geometry.b = (a/2)*1/(sqrt((1/4) + (1/9)*(sind(60)^2)));
            Settings.Geometry.c = 2*(a/2);
            Settings.Geometry.FC_Metal  = [1/3 2/3 1/2;...
                                2/3 1/3 0.0];
            Settings.Geometry.FC_Halide = [1/3 2/3 0.0;...
                                2/3 1/3 1/2];
            if Center_Coordinates
                Settings.Geometry = CenterCoordinates(Settings.Geometry,Settings.Structure);
            end
            Settings.Geometry.FC = [Settings.Geometry.FC_Metal; Settings.Geometry.FC_Halide];
        end
    case 'cscl'
        Settings.Geometry.a = a;
        Settings.Geometry.b = b;
        Settings.Geometry.c = c;
    case 'betabeo'
        Settings.Geometry.a = a;
        Settings.Geometry.b = b;
        Settings.Geometry.c = c;
    case 'sphalerite'
        if Settings.Use_Conv_cell
            Settings.Geometry.a = a;
            Settings.Geometry.b = b;
            Settings.Geometry.c = c;
        else
            Settings.Geometry.a = a/sqrt(2);
            Settings.Geometry.b = b/sqrt(2);
            Settings.Geometry.c = c/sqrt(2);
            Settings.Geometry.FC_Metal  = [0.0 0.0 0.0];
            Settings.Geometry.FC_Halide = [1/4 1/4 1/4];
            if Center_Coordinates
                Settings.Geometry = CenterCoordinates(Settings.Geometry,Settings.Structure);
            end
            Settings.Geometry.FC = [Settings.Geometry.FC_Metal; Settings.Geometry.FC_Halide];
        end
    case 'nias'
        Settings.Geometry.a = a;
        Settings.Geometry.b = b;
        Settings.Geometry.c = c;
    case 'antinias'
        Settings.Geometry.a = a;
        Settings.Geometry.b = b;
        Settings.Geometry.c = c;
end

end


