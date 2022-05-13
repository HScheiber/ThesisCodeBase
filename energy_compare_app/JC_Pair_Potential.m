function Output = JC_Pair_Potential(Startpoint,Endpoint,Spacing,Salt,Watermodel,Scaling_Params,Potential_Type)

S_D = Scaling_Params(1);
S_R = Scaling_Params(2);
S_E = Scaling_Params(3);
S_S = Scaling_Params(4);
S_MMD = Scaling_Params(5);
S_XXD = Scaling_Params(6);
S_MXD = Scaling_Params(7);
S_C = Scaling_Params(8); % Scaling charges
S_D6 = Scaling_Params(9); % Scaling R6

S_MMD = S_D*S_D6*S_MMD;
S_XXD = S_D*S_D6*S_XXD;
S_MXD = S_D*S_D6*S_MXD;

S_MMR = S_R;
S_XXR = S_R;
S_MXR = S_R;

%% Conversion factors and fundamental constants
Ang_per_m = 1e+10; % Angstroms per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(Ang_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 Angstrom^-1
k_0 = S_C/(4*pi*epsilon_0); % Coulomb constant in kJ Angstrom C^-2 mol^-1
kJ_per_kcal = 4.184; % Exact conversion factor by definition

%% Generate range (r) in Angstroms
r = Startpoint:Spacing:Endpoint;

%% Split Salt Into Component Metal and Halide
[Metal,Halide] = Separate_Metal_Halide(Salt);

%% JC Ion Parameters in SPC/E water
if strcmp(Watermodel,'SPC/E')
    switch Metal
        case 'Li'
            Metal_sigma = (0.791)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.3367344)*kJ_per_kcal; % kJ/mol
        case 'Na'
            Metal_sigma = (1.212)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.3526418)*kJ_per_kcal; % kJ/mol
        case 'K'
            Metal_sigma = (1.593)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.4297054)*kJ_per_kcal; % kJ/mol
        case 'Rb'
            Metal_sigma = (1.737)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.4451036)*kJ_per_kcal; % kJ/mol
        case 'Cs'
            Metal_sigma = (2.021)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.0898565)*kJ_per_kcal; % kJ/mol
    end

    switch Halide
        case 'F'
            Halide_sigma = (2.257)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0074005)*kJ_per_kcal; % kJ/mol
        case 'Cl'
            Halide_sigma = (2.711)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0127850)*kJ_per_kcal; % kJ/mol
        case 'Br'
            Halide_sigma = (2.751)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0269586)*kJ_per_kcal; % kJ/mol
        case 'I'
            Halide_sigma = (2.919)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0427845)*kJ_per_kcal; % kJ/mol
    end

elseif strcmp(Watermodel,'TIP3P')
    %% JC Ion Parameters in TIP3P water
    
    switch Metal
        case 'Li'
            Metal_sigma = (1.025)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.0279896)*kJ_per_kcal; % kJ/mol
        case 'Na'
            Metal_sigma = (1.369)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.0874393)*kJ_per_kcal; % kJ/mol
        case 'K'
            Metal_sigma = (1.705)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.1936829)*kJ_per_kcal; % kJ/mol
        case 'Rb'
            Metal_sigma = (1.813)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.3278219)*kJ_per_kcal; % kJ/mol
        case 'Cs'
            Metal_sigma = (1.976)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.4065394)*kJ_per_kcal; % kJ/mol
    end
    
    switch Halide
        case 'F'
            Halide_sigma = (2.303)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0033640)*kJ_per_kcal; % kJ/mol
        case 'Cl'
            Halide_sigma = (2.513)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0355910)*kJ_per_kcal; % kJ/mol
        case 'Br'
            Halide_sigma = (2.608)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0586554)*kJ_per_kcal; % kJ/mol
        case 'I'
            Halide_sigma = (2.860)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0536816)*kJ_per_kcal; % kJ/mol
    end

elseif strcmp(Watermodel,'TIP4PEW')
    %% JC Ion Parameters in TIP4P water
    
    switch Metal
        case 'Li'
            Metal_sigma = (0.808)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.1039884)*kJ_per_kcal; % kJ/mol
        case 'Na'
            Metal_sigma = (1.226)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.1684375)*kJ_per_kcal; % kJ/mol
        case 'K'
            Metal_sigma = (1.590)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.2794651)*kJ_per_kcal; % kJ/mol
        case 'Rb'
            Metal_sigma = (1.709)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.4331494)*kJ_per_kcal; % kJ/mol
        case 'Cs'
            Metal_sigma = (1.888)*(2^(5/6)); % Angstrom
            Metal_epsilon = (0.3944318)*kJ_per_kcal; % kJ/mol
    end
    
    switch Halide
        case 'F'
            Halide_sigma = (2.538)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0015752)*kJ_per_kcal; % kJ/mol
        case 'Cl'
            Halide_sigma = (2.760)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0116615)*kJ_per_kcal; % kJ/mol
        case 'Br'
            Halide_sigma = (2.768)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0303773)*kJ_per_kcal; % kJ/mol
        case 'I'
            Halide_sigma = (2.952)*(2^(5/6)); % Angstrom
            Halide_epsilon = (0.0417082)*kJ_per_kcal; % kJ/mol
    end
else
    error(['Unknown Water Model: "' Watermodel ...
        '". Please choose one of "SPC/E", "TIP3P", or "TIP4PEW".'])
end

% Scale stuff
Metal_epsilon = S_E*Metal_epsilon;
Halide_epsilon = S_E*Halide_epsilon;
Halide_sigma = S_S*Halide_sigma;
Metal_sigma = S_S*Metal_sigma;

% Use combining rules
Salt_epsilon = sqrt(Metal_epsilon*Halide_epsilon); % kJ/mol
Salt_Sigma = (1/2)*(Metal_sigma + Halide_sigma);

q_M = 1*e_c; % charge of alkali metal in C
q_H = -1*e_c; % charge of halide in C

%% Build JC PES for LiF with dispersion scaling
if strcmpi(Potential_Type,'full')
    Output.MX = k_0.*q_M.*q_H./(r) + 4*Salt_epsilon.*( S_MXR.*( Salt_Sigma./r ).^12 - S_MXD*( Salt_Sigma./r ).^6 );
    Output.MM = k_0.*q_M.*q_M./(r) + 4*Metal_epsilon.*( S_MMR.*( Metal_sigma./r ).^12 - S_MMD*( Metal_sigma./r ).^6 );
    Output.XX = k_0.*q_H.*q_H./(r) + 4*Halide_epsilon.*( S_XXR.*( Halide_sigma./r ).^12 - S_XXD*( Halide_sigma./r ).^6 );
elseif strcmpi(Potential_Type,'vdw')
    Output.MX = 4*Salt_epsilon.*( S_MXR.*( Salt_Sigma./r ).^12 - S_MXD*( Salt_Sigma./r ).^6 );
    Output.MM = 4*Metal_epsilon.*( S_MMR.*( Metal_sigma./r ).^12 - S_MMD*( Metal_sigma./r ).^6 );
    Output.XX = 4*Halide_epsilon.*( S_XXR.*( Halide_sigma./r ).^12 - S_XXD*( Halide_sigma./r ).^6 );
elseif strcmpi(Potential_Type,'dispersion')
    Output.MX = 4*Salt_epsilon.*( - S_MXD*( Salt_Sigma./r ).^6 );
    Output.MM = 4*Metal_epsilon.*( - S_MMD*( Metal_sigma./r ).^6 );
    Output.XX = 4*Halide_epsilon.*( - S_XXD*( Halide_sigma./r ).^6 );
elseif strcmpi(Potential_Type,'repulsion')
    Output.MX = S_MXR*4*Salt_epsilon.*( ( Salt_Sigma./r ).^12 );
    Output.MM = S_MMR*4*Metal_epsilon.*( ( Metal_sigma./r ).^12 );
    Output.XX = S_XXR*4*Halide_epsilon.*( ( Halide_sigma./r ).^12 );
end
Output.r = r;
end