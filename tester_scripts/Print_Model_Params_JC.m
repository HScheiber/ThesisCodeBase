Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
Watermodels = {'SPC/E' 'TIP3P' 'TIP4PEW' 'SD'};

kj_per_erg = 1e-10; % kJ per erg
nm_per_cm = 1e+7; % nm per cm
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
C_unit = 1e-60; % erg cm^6
D_unit = 1e-76; % erg cm^8
nm_per_Ang = 0.1; % nm per Angstrom

Settings = Initialize_MD_Settings;

for jdx = 1:length(Watermodels)
    Settings.WaterModel = Watermodels{jdx};
    for idx = 1:length(Salts)
        Settings.Salt = Salts{idx};
        
        if strcmp(Settings.WaterModel,'SD') && ~strcmpi(Settings.Salt,'NaCl')
            continue
        end
        
        [Metal,Halide] = Separate_Metal_Halide(Settings.Salt);

        [OutputMX,OutputMM,OutputXX] = JC_Potential_Parameters(Settings);

        disp(repmat('*',1,30))
        disp(['JC (' Settings.WaterModel ') Parameters for ' Settings.Salt])
        disp(repmat('*',1,30))
        disp('Sigma_MM      Epsilon_MM      Sigma_XX        Epsilon_XX')

        sigma_MM = num2str(OutputMM.sigma/nm_per_Ang,'%.3f'); % Angstrom
        sigma_XX = num2str(OutputXX.sigma/nm_per_Ang,'%.3f'); % Angstrom
        epsilon_MM = num2str(OutputMM.epsilon,'%.6f'); % Angstrom
        epsilon_XX = num2str(OutputXX.epsilon,'%.6f'); % Angstrom


        disp([sigma_MM '         ' epsilon_MM '         ' sigma_XX '           ' epsilon_XX])
    end
end