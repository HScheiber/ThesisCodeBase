% Bug fixing tool
Cry.BetaBeO.N = 8;
Cry.CsCl.N = 2;
Cry.FiveFive.N = 8;
Cry.NiAs.N = 4;
Cry.Rocksalt.N = 8;
Cry.Sphalerite.N = 8;
Cry.Wurtzite.N = 4;
Directory = pwd;


if ispc
    Parcores = 12;
    if isempty(gcp('nocreate'))
        parpool(Parcores);
    end
elseif isunix
    if isempty(gcp('nocreate'))
        [~,coretxt]= system('echo $SLURM_CPUS_ON_NODE');
        Parcores = str2double(coretxt);
        if isnan(Parcores)
            Parcores = 6;
        end    
        MyProfile = parcluster('local');
        MyProfile.NumWorkers = Parcores;
        MyProfile.NumThreads = Parcores;
        MyProfile.saveProfile
        poolobj = parpool(MyProfile,Parcores);
    end
end

Salts = find_directories(dir(Directory));
N = length(Salts);
% Loop through folders
for i = 1:N
    Salt = Salts{i};
    Salt_Directory = fullfile(Directory,Salt);

    % Find directories in current folder, except . and ..
    Structures = find_directories(dir(Salt_Directory));
    M = length(Structures);
    
    for j = 1:M
        Structure = Structures{j};
        Structure_Directory = fullfile(Directory,Salt,Structure);
        
        % Find directories in current folder, except . and ..
        Models = find_directories(dir(Structure_Directory));
        R = length(Models);
        
        Current_Str_Dat = Cry.(Structure).N;
        
        for k = 1:R
            Model = Models{k};
            Model_Directory = fullfile(Directory,Salt,Structure,Model);

            % Find directories in current folder, except . and ..
            LatPars = find_directories(dir(Model_Directory));
            P = length(LatPars);
            
            parfor m = 1:P
                Latpar = LatPars{m};
                Latpar_Directory = fullfile(Directory,Salt,Structure,Model,Latpar);
                
                % Find mat file
                matfile = strtrim(ls([Latpar_Directory filesep '*.mat']));
                
                X = load(matfile,'N_Cell','N_total');
                
                N_Cell = Current_Str_Dat;
                N_total = (X.N_total/8)*N_Cell;
                
                parallelsave(matfile,{'N_Cell' 'N_total'},N_Cell,N_total);
                disp([Salt ' ' Structure ' ' Model ' ' Latpar ' Fixed.'])
            end
        end
    end
end