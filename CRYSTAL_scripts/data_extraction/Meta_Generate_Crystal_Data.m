% Short script to generate, and save data
% Meta_Generate_Crystal_Data
%% Parameters
Salts = {'NaCl'}; %{'LiLi' 'FF' 'ClCl' 'BrBr' 'II'}; %{'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
Structures = {'Wurtzite'}; %{'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive' 'Pair'};
ReCalc_Salts = true; % Set to false to prevent re-scraping of salt total energies from raw data
ReCalc_RefAtom = false; % Set to false to prevent re-scraping of Lone ion/atom energies from raw data
ReCalc_BSSE = false; % Set to false to prevent re-scraping of BSSE Data from raw data
parallel = false; %#ok<*UNRCH>

R = length(Salts);

analysis_folder_PC = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\analysis';
analysis_folder_unix = '/home/scheiber/project/CRYSTAL_LiHalides';
analysis_folder_LabPC = '/home/user/project/CRYSTAL_LiHalides';
analysis_folder_Sockeye = '/home/haydensc/scratch/CRYSTAL_LiHalides';

datadir_PC = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\CRYSTAL';
datadir_unix = '/home/scheiber/ThesisCodeBase/data/CRYSTAL';
datadir_LabPC = '/home/user/ThesisCodeBase/data/CRYSTAL';
datadir_Sockeye = '/home/haydensc/ThesisCodeBase/data/CRYSTAL';

if ispc
    analysis_folder = analysis_folder_PC;
	datadir = datadir_PC;
    profile = 'local';
    addpath(analysis_folder);
elseif isunix
    [~,Servertxt] = system('hostname -s | cut -c 1-3');
    Server = strtrim(Servertxt);
    if strcmp(Server,'pat')
        analysis_folder = analysis_folder_LabPC;
        datadir = datadir_LabPC;
        profile = 'MATLAB_DCS_PROFILE';
    elseif ~isempty(regexp(Server,'se[0-9]','ONCE')) || strcmpi(Server,'log') % sockeye
        analysis_folder = analysis_folder_Sockeye;
        datadir = datadir_Sockeye;
        profile = 'MATLAB_DCS_PROFILE';
    else
        analysis_folder = analysis_folder_unix;
        datadir = datadir_unix;
        profile = 'MATLAB_DCS_PROFILE';
    end
end

cd(analysis_folder);

% Start parallel pool
if parallel
    if isempty(gcp('nocreate'))
        [~,coretxt]= system('echo $SLURM_CPUS_ON_NODE');
        Parcores = str2double(coretxt);
        if isnan(Parcores) || Parcores > R
            Parcores = R;
        end    
        MyProfile = parcluster('local');
        MyProfile.NumWorkers = Parcores;
        MyProfile.NumThreads = Parcores;
        MyProfile.saveProfile
        while true
            try
                poolobj = parpool(MyProfile,Parcores);
                break
            catch
                if Parcores < 1
                    error('Too few cores selected')
                else
                    Parcores = Parcores-1;
                end
            end
        end
    end
end

% To prevent re-doing elements
Completed_Elements = cell(1,0);
if parallel
    parfor i = 1:length(Salts)
        GD_Innerloop(Salts{i},Structures,datadir,ReCalc_Salts,ReCalc_RefAtom,ReCalc_BSSE,{});
    end
    delete(poolobj)
else
    for i = 1:length(Salts)
        Completed_Elements = GD_Innerloop(Salts{i},Structures,datadir,ReCalc_Salts,ReCalc_RefAtom,...
            ReCalc_BSSE,Completed_Elements);
    end
end

