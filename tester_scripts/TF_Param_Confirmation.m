Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
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

for idx = 1:length(Salts)
    Settings.Salt = Salts{idx};

    [OutputMX,OutputMM,OutputXX] = TF_Potential_Parameters(Settings);
    
    disp(repmat('*',1,30))
    disp(['TF: Parameters for ' Settings.Salt])
    disp(repmat('*',1,30))
    disp('pij	A_MX	A_MM	A_XX    C_MX    C_MM    C_XX    D_MX    D_MM    D_XX')
    
    pij = num2str((1/OutputMM.alpha)/nm_per_Ang,'%.3f'); % Angstrom
    A_MX = num2str(OutputMX.B,'%.1f'); % kJ/mol
    A_MM = num2str(OutputMM.B,'%.1f');
    A_XX = num2str(OutputXX.B,'%.1f');
    
    C_MX = num2str((OutputMX.C)/(nm_per_Ang^6),'%.1f'); % kJ/mol
    C_MM = num2str(OutputMM.B,'%.1f');
    C_XX = num2str(OutputXX.B,'%.1f');

    disp(['Aij: ' num2str(OutputMM.B,'%.1f') ' kJ/mol'])
    disp(['pij: ' num2str((1/OutputMM.alpha)/nm_per_Ang,'%.3f')])
    disp(['Cij: ' num2str((OutputMM.C)/(nm_per_Ang^6),'%.1f')])
    disp(['Dij: ' num2str((OutputMM.D)/(nm_per_Ang^8),'%.1f')])

    disp(repmat('*',1,30))
    disp(['TF: MX Interaction for ' Settings.Salt])
    disp(repmat('*',1,30))
    disp(['Aij: ' num2str(OutputMX.B,'%.1f') ' kJ/mol'])
    disp(['pij: ' num2str((1/OutputMX.alpha)/nm_per_Ang,'%.3f')])
    disp(['Cij: ' num2str((OutputMX.C)/(nm_per_Ang^6),'%.1f')])
    disp(['Dij: ' num2str((OutputMX.D)/(nm_per_Ang^8),'%.1f')])

    disp(repmat('*',1,30))
    disp(['TF: XX Interaction for ' Settings.Salt])
    disp(repmat('*',1,30))
    disp(['Aij: ' num2str(OutputXX.B,'%.1f') ' kJ/mol'])
    disp(['pij: ' num2str((1/OutputXX.alpha)/nm_per_Ang,'%.3f')])
    disp(['Cij: ' num2str((OutputXX.C)/(nm_per_Ang^6),'%.1f')])
    disp(['Dij: ' num2str((OutputXX.D)/(nm_per_Ang^8),'%.1f')])
end