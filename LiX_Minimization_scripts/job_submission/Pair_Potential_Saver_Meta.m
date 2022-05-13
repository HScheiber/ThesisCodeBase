Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
DataDir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\GROMACS';
Dispersion_Scale = 0:0.05:5;
Epsilon_Scale = 1.0;

for i=1:length(Salts)
    Salt = Salts{i};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    MetMet = [Metal Metal];
    HalHal = [Halide Halide];
    
    Data = Pair_Potential_Saver(0.1,0.5,0.001,Salt,'SPC/E',Dispersion_Scale,Epsilon_Scale);
    
    SaveDirSalt = fullfile(DataDir,[Salt '_Pair_Lattice_Energies.mat']);
    SaveDirMetMet = fullfile(DataDir,[MetMet '_Pair_Lattice_Energies.mat']);
    SaveDirHalHal = fullfile(DataDir,[HalHal '_Pair_Lattice_Energies.mat']);
    
    save(SaveDirSalt,'Data')
    
    FullData = Data;
    
    % Met Met Data
    Data = rmfield(Data,Salt);
    Data = rmfield(Data,HalHal);
    save(SaveDirMetMet,'Data')
    
    % Hal Hal Data
    Data = FullData;
    Data = rmfield(Data,Salt);
    Data = rmfield(Data,MetMet);
    save(SaveDirHalHal,'Data')
    
end