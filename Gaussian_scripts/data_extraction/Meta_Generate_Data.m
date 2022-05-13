% Short script to generate, and save data
%% Parameters
Salt = 'LiCl';

[match,~] = regexp(Salt,{'Li' 'F' 'Cl' 'Br' 'I'},'match','tokens');
matches = match(~cellfun('isempty',match));
metal = matches{1}{1};
halide = matches{2}{1};
datadir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data';

%% Get ion energies
Data = Extract_Gaussian_Energies(metal);
filename = [metal '_Ion_Total_Energies_DFT.mat'];
save([datadir filesep filename],'Data')
clearvars Data filename

Data = Extract_Gaussian_Energies(halide);
filename = [halide '_Ion_Total_Energies_DFT.mat'];
save([datadir filesep filename],'Data')
clearvars Data filename

%% Get Rocksalt energies
cd([Salt filesep 'Rocksalt'])
Label = [Salt '_R'];
crystal_type = 'Rocksalt';
Data = Extract_Gaussian_DFT_PES(Label);
filename = [Salt '_Rocksalt_Total_Energies.mat'];
save([datadir filesep filename],'Data')
clearvars Data filename

filename = [Salt '_Rocksalt_Lattice_Energies.mat'];
Data = Generate_Lattice_energy(Salt,crystal_type,datadir);
save([datadir filesep filename],'Data')
clearvars Data filename
cd(['..' filesep '..'])
    
%% Get Wurtzite energies
cd([Salt filesep 'Wurtzite'])
Label = [Salt '_W'];
crystal_type = 'Wurtzite';
Data = Extract_Gaussian_DFT_PES(Label);
filename = [Salt '_Wurtzite_Total_Energies.mat'];
save([datadir filesep filename],'Data')
clearvars Data filename

filename = [Salt '_Wurtzite_Lattice_Energies.mat'];
Data = Generate_Lattice_energy(Salt,crystal_type,datadir);
save([datadir filesep filename],'Data')
cd(['..' filesep '..'])
clear