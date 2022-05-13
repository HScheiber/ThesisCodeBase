function [sigma_cross,epsilon_cross] = grab_clayff(Atom_A,Atom_B,Salt)
Settings = Initialize_MD_Settings;
Settings.Salt = Salt;
[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);
Atoms = {Atom_A Atom_B};
kJ_per_kcal = 4.184;
nm_per_Ang = 0.1;
R0_to_Sigma = (2^(-1/6)); % Conversion factor from R0 to sigma
% S.D.MM = 910;
% S.D.XX = 0.23;


switch Settings.Theory
    case 'JC'
        Settings.WaterModel = 'SPC/E';
    case 'JC3P'
        Settings.WaterModel = 'TIP3P';
    case 'JC4P'
        Settings.WaterModel = 'TIP4PEW';
    case 'JCSD'
        Settings.WaterModel = 'SD';
end

[OutputMM,OutputXX,~] = JC_Potential_Parameters(Settings);



sigma = zeros(2,1);
epsilon = zeros(2,1);
charge = zeros(2,1);
for idx = 1:2
    Atom = Atoms{idx};
    
    switch lower(Atom)
        case 'o'
            sigma(idx) = 3.5532*nm_per_Ang*R0_to_Sigma; % nm
            epsilon(idx) = 0.1554*kJ_per_kcal; % kJ/mol
            charge(idx) = -1.0500; % Charge in electron charge units
        case 'al'
            sigma(idx) = 4.7943*nm_per_Ang*R0_to_Sigma; % nm
            epsilon(idx) = 1.3298e-6*kJ_per_kcal; % kJ/mol
            charge(idx) = 1.5750; % Charge in electron charge units
        case lower(Metal)
            sigma(idx) = OutputMM.sigma;
            epsilon(idx) = OutputMM.epsilon;
            charge(idx) = 1;
        case lower(Halide)
            sigma(idx) = OutputXX.sigma;
            epsilon(idx) = OutputXX.epsilon;
            charge(idx) = -1;
            
    end
    
    if strcmp(Atom_A,Atom_B)
        sigma_cross = sigma(idx);
        epsilon_cross = epsilon(idx);
        return
    end
end

sigma_cross = sum(sigma)/2;
epsilon_cross = sqrt(prod(epsilon));

end
