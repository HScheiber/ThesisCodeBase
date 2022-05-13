% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETERS AND kJ/mol

% C6Damping:
% 0 = no (default) damping
% 1 = BJ/rational damping (same as in D3(BJ), damps to a constant. Fairly
% weak damping)
% 2 = Tang Damping (Mid strength damping, damps to zero)
% 3 = MMDRE Damping function (very weak damping)
% 4 = PAMoC Damping function (weak damping)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)

% C6Damping adds close-range damping

% Parameter sets for C6/C8 coefficients
% 0 = default TF Parameters
% 1 = D3 values
% 2 = Best literature values available
% 3 = D4 with C6 and C8 generated on the fly

% GAdjust are N x 3 arrays of gaussian parameters
% (i , 1) is the Gaussian height of the ith adjustment (may be negative or
% positive)
% (i , 2) is the center point of the ith Gaussian (should be positive)
% (i , 3) is the standard deviation or width (negative and positive values
% are the same)

function [U_MX_out, U_MM_out, U_XX_out,d4fail] = TF_Potential_Generator(Settings,varargin)

% Optional inputs
p = inputParser;
p.FunctionName = 'TF_Potential_Generator';
addOptional(p,'PlotType','full')
addOptional(p,'ReturnAsStructure',false);
addOptional(p,'Startpoint',0);
addOptional(p,'Plotswitch',false);
addOptional(p,'MDP_Minimize',false);
addOptional(p,'FactorOutParams',false);
parse(p,varargin{:});
PlotType = p.Results.PlotType;
ReturnAsStructure = p.Results.ReturnAsStructure;
Startpoint = p.Results.Startpoint;
Plotswitch = p.Results.Plotswitch;
FactorOutParams = p.Results.FactorOutParams;

% Allowed plot types: 'full', 'lj', 'full-derivative', 'lj-derivative',
% 'dispersion', 'dispersion-derivative', 'repulsive',
% 'repulsive-derivative'

if p.Results.MDP_Minimize
    MDP = 'MinMDP';
else
    MDP = 'MDP';
end

d4fail = false;

%% Gaussian adjustments
G_a.MM = Settings.GAdjust_MM(:,1);
G_b.MM = Settings.GAdjust_MM(:,2);
G_c.MM = Settings.GAdjust_MM(:,3);

G_a.XX = Settings.GAdjust_XX(:,1);
G_b.XX = Settings.GAdjust_XX(:,2);
G_c.XX = Settings.GAdjust_XX(:,3);

G_a.MX = Settings.GAdjust_MX(:,1);
G_b.MX = Settings.GAdjust_MX(:,2);
G_c.MX = Settings.GAdjust_MX(:,3);


%% Conversion factors and fundamental constants
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

%% Split Salt Into Component Metal and Halide
[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);

if Settings.TF_Paramset == 1 % D3 dispersion Parameters
    %% C6 Coefficients
    C.LiF.MM  = 4918.9042107878*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiF.MX  = 990.0857326400*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiF.XX  = 411.2995536808*(nm_per_Ang^6); % kJ/mol nm^6
   
    C.LiCl.MM = 4918.9042107878*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiCl.MX = 4450.9148362278*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiCl.XX = 5211.7103353493*(nm_per_Ang^6); % kJ/mol nm^6
    
    C.LiBr.MM = 4918.9042107878*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiBr.MX = 6238.2901762131*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiBr.XX = 9756.9852141205*(nm_per_Ang^6); % kJ/mol nm^6
    
    C.LiI.MM = 4918.9042107878*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiI.MX = 9306.8796821688*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiI.XX = 20668.4353099615*(nm_per_Ang^6); % kJ/mol nm^6

    C.NaCl.MM  = 10729.4523062025*(nm_per_Ang^6); % kJ/mol nm^6
    C.NaCl.MX  = 6456.1075384258*(nm_per_Ang^6); % kJ/mol nm^6
    C.NaCl.XX  = 5211.710335349*(nm_per_Ang^6); % kJ/mol nm^6
    
    %% C8 coefficients
    D.LiF.MM  = 104130.2216951388*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiF.MX  = 9971.6966373607*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiF.XX  = 1970.7989468780*(nm_per_Ang^8); % kJ/mol nm^8

    D.LiCl.MM = 104130.2216951388*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiCl.MX = 69999.5686578376*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiCl.XX = 60892.5274692634*(nm_per_Ang^8); % kJ/mol nm^8

    D.LiBr.MM = 104130.2216951388*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiBr.MX = 120775.5671985948*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiBr.XX = 172756.4400846956*(nm_per_Ang^8); % kJ/mol nm^8

    D.LiI.MM = 104130.2216951388*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiI.MX = 217169.0376048321*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiI.XX = 531602.2210887696*(nm_per_Ang^8); % kJ/mol nm^8

    D.NaCl.MM  = 390953.8836355778*(nm_per_Ang^8); % kJ/mol nm^8
    D.NaCl.MX  = 133209.9331112395*(nm_per_Ang^8); % kJ/mol nm^8
    D.NaCl.XX  = 60892.5274692634*(nm_per_Ang^8); % kJ/mol nm^8

elseif Settings.TF_Paramset == 2 % Best Literature Parameters from doi:10.1080/00268978600102091

    C.LiF.MM  = 4.44*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiF.MX  = 61.46*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiF.XX  = 1103.47*(nm_per_Ang^6); % kJ/mol nm^6

    C.LiCl.MM = 4.44*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiCl.MX = 146.67; % kJ/mol nm^6
    C.LiCl.XX = 8446.11; % kJ/mol nm^6

    C.LiBr.MM = 4.44*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiBr.MX = (2.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiBr.XX = (185*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.LiI.MM = 4.44*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiI.MX = (3.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiI.XX = (378*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.NaCl.MM  = 87.29*(nm_per_Ang^6); % kJ/mol nm^6
    C.NaCl.MX  = 716.62*(nm_per_Ang^6); % kJ/mol nm^6
    C.NaCl.XX  = 9195.59*(nm_per_Ang^6); % kJ/mol nm^6

    %% C8 Coefficients
    D.LiF.MM  = 4.05*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiF.MX  = 184.21*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiF.XX  = 5695.74*(nm_per_Ang^8); % kJ/mol nm^8

    D.LiCl.MM = 4.05*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiCl.MX = 760.56*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiCl.XX = 75265.10*(nm_per_Ang^8); % kJ/mol nm^8

    D.LiBr.MM = 4.05*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiBr.MX = (3.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiBr.XX = (423*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.LiI.MM = 4.05*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiI.MX = (5.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiI.XX = (1060*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.NaCl.MM  = 198.90*(nm_per_Ang^8); % kJ/mol nm^8
    D.NaCl.MX  = 4388.04*(nm_per_Ang^8); % kJ/mol nm^8
    D.NaCl.XX  = 87244.23*(nm_per_Ang^8); % kJ/mol nm^8
    
elseif Settings.TF_Paramset == 3 % D4 on-the-fly parameters
    % Conversion factors for D4
    Bohr_nm = 0.05291772108; % a_0 - > nm
    c6conv = 1e-3/2625.4999/((0.052917726)^6); % J/mol nm^6 - > au (from D3 sourcecode)
    J_kJ = 1e-3; % J - > kJ
    Ha_kJmol = 2625.4999; % Ha - > kJ/mol
    c6units = (1/c6conv)*J_kJ; % au - > kJ/mol nm^6
    c8units = (Ha_kJmol)*(Bohr_nm^8); % au - > kJ/mol nm^8

    % Create vasp format geometry file in given Directory
    N = Settings.Geometry.N;   
    geom_txt = ['New Structure' newline '1.0' newline];
    TM = Settings.Geometry.Transform*[Settings.Geometry.a 0 0; 0 Settings.Geometry.b 0; 0 0 Settings.Geometry.c];
    for i = 1:3
        for j = 1:3
            geom_txt = [geom_txt pad(num2str(TM(i,j),'%10.10f'),20,'left')];
        end
        geom_txt = [geom_txt newline];
    end
    Nd2 = num2str(N/2);

    geom_txt = [geom_txt pad(Metal,5,'left') pad(Halide,5,'left') newline];
    geom_txt = [geom_txt pad(Nd2,5,'left') pad(Nd2,5,'left') newline 'Direct' newline];
    
    for i = 1:(N/2)
        for j = 1:3
            geom_txt = [geom_txt pad(num2str(Settings.Geometry.FC_Metal(i,j),'%10.10f'),16,'left')];
        end
        geom_txt = [geom_txt newline];
    end
    for i = 1:(N/2)
        for j = 1:3
            geom_txt = [geom_txt pad(num2str(Settings.Geometry.FC_Halide(i,j),'%10.10f'),16,'left')]; %#ok<*AGROW>
        end
        geom_txt = [geom_txt newline];
    end
    
    % save to file
    filename = [Settings.WorkDir filesep 'tempstruct.vasp'];
    fidPM = fopen(filename,'wt');
    fwrite(fidPM,regexprep(geom_txt,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidPM);
    
    % run dftd4 on the structure
    if ispc % for testing
        dftd4 = 'wsl source ~/.bashrc; dftd4 ';
        fnunix = windows2unix(filename);
    elseif isunix
        dftd4 = 'dftd4 ';
        fnunix = filename;
    end
    [ercode,dftd4_out] = system([dftd4 fnunix]);
    
    if ercode ~= 0 %dftd4 failed
        error(['DFTD4 module failed from input: '  dftd4 fnunix])
    end
    
    % Grab C6's from output (in AU)
    C6s = regexp(dftd4_out,'# +Z +covCN +q +C6AA.+\n\n\n','match','ONCE');
    MM_C6 = regexp(C6s,[Metal ' +[-.0-9]+ +[-.0-9]+ +([-.0-9]+) +'],'tokens','once');
    if isempty(MM_C6) % check for D4 fail
        d4fail = true;
        U_MX_out = [];
        U_MM_out = [];
        U_XX_out = [];
        return
    else
        C6_MM = str2double(MM_C6{1});
    end
    XX_C6 = regexp(C6s,[Halide ' +[-.0-9]+ +[-.0-9]+ +([-.0-9]+) +'],'tokens','once');
    if isempty(XX_C6) % check for D4 fail
        d4fail = true;
        U_MX_out = [];
        U_MM_out = [];
        U_XX_out = [];
        return
    else
        C6_XX = str2double(XX_C6{1});
    end
    
    Mol_C6 = regexp(C6s,'Mol\. C6AA.+? + : +([-.0-9]+)','tokens','once');
    C6_Mol = str2double(Mol_C6{1});
    
    % Calculate cross term C6
    C6_MX = (C6_Mol - ((N/2)^2)*C6_MM - ((N/2)^2)*C6_XX)/(2*(N/2)^2);
    
    % Load r2r4 from disc (for calculation of C8)
    loadr2r4 = load('r2r4.mat','r2r4');
    r2r4 = loadr2r4.r2r4;

    % Atomic numbers of atoms in salt
    Z_M = elements('Symbol',Metal,'atomic_number');
    Z_X = elements('Symbol',Halide,'atomic_number');

    %% Calculate C8 for each possible pair of atoms in unit cell, as well as R0AB
    sqrt_Q_M = r2r4(Z_M); % Factor used to calculate C8 for Metal
    sqrt_Q_X = r2r4(Z_X); % Factor used to calculate C8 for Halide

    % C6 coefficients
    C.(Settings.Salt).MX = C6_MX*c6units; % in kJ/mol nm^6
    C.(Settings.Salt).MM = C6_MM*c6units; % in kJ/mol nm^6
    C.(Settings.Salt).XX = C6_XX*c6units; % in kJ/mol nm^6
    
    % C8 coefficients
    D.(Settings.Salt).MX = 3.0*(C6_MX)*sqrt_Q_M*sqrt_Q_X*c8units; % in kJ/mol nm^8
    D.(Settings.Salt).MM = 3.0*(C6_MM)*sqrt_Q_M*sqrt_Q_M*c8units; % in kJ/mol nm^8
    D.(Settings.Salt).XX = 3.0*(C6_XX)*sqrt_Q_X*sqrt_Q_X*c8units; % in kJ/mol nm^8
    
else % Default parameters
    %% Huggins-Mayer Dipole-Dipole Dispersion Parameter C: MX = +-   MM = ++     XX = --
    C.LiF.MM  = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiF.MX  = (0.8*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiF.XX  = (14.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.LiCl.MM = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiCl.MX = (2.0*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiCl.XX = (111*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.LiBr.MM = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiBr.MX = (2.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiBr.XX = (185*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.LiI.MM = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiI.MX = (3.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiI.XX = (378*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.NaF.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaF.MX  = (4.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaF.XX  = (16.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.NaCl.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaCl.MX  = (11.2*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaCl.XX  = (116*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.NaBr.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaBr.MX  = (14.0*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaBr.XX  = (196*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.NaI.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaI.MX  = (19.1*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaI.XX  = (392*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.KF.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KF.MX  = (19.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KF.XX  = (18.6*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.KCl.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KCl.MX  = (48*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KCl.XX  = (124.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.KBr.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KBr.MX  = (60*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KBr.XX  = (206*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.KI.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KI.MX  = (82*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KI.XX  = (403*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.RbF.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbF.MX  = (31*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbF.XX  = (18.9*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.RbCl.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbCl.MX  = (79*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbCl.XX  = (130*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.RbBr.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbBr.MX  = (99*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbBr.XX  = (215*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.RbI.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbI.MX  = (135*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbI.XX  = (428*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.CsF.MM  = (152*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.CsF.MX  = (52*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.CsF.XX  = (19.1*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    
    C.CsCl.MM  = (152*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.CsCl.MX  = (129*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.CsCl.XX  = (129*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    
    C.CsBr.MM  = (152*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.CsBr.MX  = (163*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.CsBr.XX  = (214*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    
    C.CsI.MM  = (152*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.CsI.MX  = (224*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.CsI.XX  = (424*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    %% Huggins-Mayer Dipole-Quadrupole Dispersion Parameter D: MX = +-   MM = ++     XX = --
    D.LiF.MM  = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiF.MX   = (0.6*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiF.XX    = (17*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.LiCl.MM = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiCl.MX = (2.4*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiCl.XX   = (223*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.LiBr.MM = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiBr.MX = (3.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiBr.XX = (423*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.LiI.MM = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiI.MX = (5.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiI.XX = (1060*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.NaF.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaF.MX  = (3.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaF.XX  = (20*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.NaCl.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaCl.MX  = (13.9*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaCl.XX  = (233*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.NaBr.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaBr.MX  = (19*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaBr.XX  = (450*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.NaI.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaI.MX  = (31*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaI.XX  = (1100*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.KF.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KF.MX  = (21*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KF.XX  = (22*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.KCl.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KCl.MX  = (73*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KCl.XX  = (250*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.KBr.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KBr.MX  = (99*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KBr.XX  = (470*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.KI.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KI.MX  = (156*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KI.XX  = (1130*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.RbF.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbF.MX  = (40*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbF.XX  = (23*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.RbCl.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbCl.MX  = (134*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbCl.XX  = (260*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.RbBr.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbBr.MX  = (180*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbBr.XX  = (490*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.RbI.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbI.MX  = (280*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbI.XX  = (1200*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.CsF.MM  = (278*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.CsF.MX  = (78*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.CsF.XX  = (23*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    
    D.CsCl.MM  = (278*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.CsCl.MX  = (250*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.CsCl.XX  = (260*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    
    D.CsBr.MM  = (278*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.CsBr.MX  = (340*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.CsBr.XX  = (490*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    
    D.CsI.MM  = (278*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.CsI.MX  = (520*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.CsI.XX  = (1190*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
end

%% TF Repulsive Size Parameter sigma (AKA r+/-): P = +   M = -
% Metals
sigma.Li = 0.816*nm_per_Ang; % nm
sigma.Na = 1.170*nm_per_Ang; % nm
sigma.K  = 1.463*nm_per_Ang; % nm
sigma.Rb = 1.587*nm_per_Ang; % nm
sigma.Cs = 1.720*nm_per_Ang; % nm

% Halides
sigma.F  = 1.179*nm_per_Ang; % nm
sigma.Cl = 1.585*nm_per_Ang; % nm
sigma.Br = 1.716*nm_per_Ang; % nm
sigma.I  = 1.907*nm_per_Ang; % nm

%% TF Parameter: Number of Valence electrons (for Pauling Coefficient Calculation)
% Metals
valence.Li = 2;
valence.Na = 8;
valence.K = 8;
valence.Rb = 8;
valence.Cs = 8;

% Halides
valence.F = 8;
valence.Cl = 8;
valence.Br = 8;
valence.I = 8;

%% TF Hardness Parameter Rho
rho.LiF = 0.299*nm_per_Ang; % nm
rho.LiCl = 0.342*nm_per_Ang; % nm
rho.LiBr = 0.353*nm_per_Ang; % nm
rho.LiI = 0.430*nm_per_Ang; % nm

rho.NaF = 0.330*nm_per_Ang; % nm
rho.NaCl = 0.317*nm_per_Ang; % nm
rho.NaBr = 0.340*nm_per_Ang; % nm
rho.NaI = 0.386*nm_per_Ang; % nm

rho.KF = 0.338*nm_per_Ang; % nm
rho.KCl = 0.337*nm_per_Ang; % nm
rho.KBr = 0.335*nm_per_Ang; % nm
rho.KI = 0.355*nm_per_Ang; % nm

rho.RbF = 0.328*nm_per_Ang; % nm
rho.RbCl = 0.318*nm_per_Ang; % nm
rho.RbBr = 0.335*nm_per_Ang; % nm
rho.RbI = 0.337*nm_per_Ang; % nm

rho.CsF = 0.282*nm_per_Ang; % nm
rho.CsCl = 0.272*nm_per_Ang; % nm THIS IS ESTIMATED
rho.CsBr = 0.289*nm_per_Ang; % nm THIS IS ESTIMATED
rho.CsI = 0.291*nm_per_Ang; % nm THIS IS ESTIMATED

%% TF Parameter: q (charge)
q.Li =  Settings.S.Q; % atomic
q.Na =  Settings.S.Q; % atomic
q.K  =  Settings.S.Q; % atomic
q.Rb =  Settings.S.Q; % atomic
q.Cs =  Settings.S.Q; % atomic

q.F  = -Settings.S.Q; % atomic
q.Cl = -Settings.S.Q; % atomic
q.Br = -Settings.S.Q; % atomic
q.I  = -Settings.S.Q; % atomic

%% Huggins-Mayer potential parameter b (same for all salts)
b = (0.338e-12)*kj_per_erg*NA; % kJ/mol

%% Generate range (r) in nm
r = Startpoint:Settings.Table_StepSize:Settings.Table_Length;

%% Calculate Pauling Coefficients beta: MX = +-   MM = ++     XX = --
beta.MM = 1 + 2*q.(Metal)/valence.(Metal); % Unitless
beta.MX = 1 + q.(Metal)/valence.(Metal) + q.(Halide)/valence.(Halide); % Unitless
beta.XX = 1 + 2*q.(Halide)/valence.(Halide); % Unitless

%% Calculate TF Repulsive Exponential Parameter alpha: MX = +-   MM = ++     XX = --
alpha.MM = Settings.S.A.All*Settings.S.A.MM/rho.(Settings.Salt); % nm^-1
alpha.MX = Settings.S.A.All*Settings.S.A.MX/rho.(Settings.Salt); % nm^-1
alpha.XX = Settings.S.A.All*Settings.S.A.XX/rho.(Settings.Salt); % nm^-1

%% Calculate TF Repulsive Scaling Parameter B: MX = +-   MM = ++     XX = -- (Including scaling)
B.MM = Settings.S.R.All*Settings.S.R.MM*beta.MM*b*exp(2*sigma.(Metal)/rho.(Settings.Salt));
B.MX = Settings.S.R.All*Settings.S.R.MX*beta.MX*b*exp((sigma.(Metal) + sigma.(Halide))/rho.(Settings.Salt));
B.XX = Settings.S.R.All*Settings.S.R.XX*beta.XX*b*exp(2*sigma.(Halide)/rho.(Settings.Salt));

%% Scale Dispersion
C.(Settings.Salt).MM = Settings.S.D6D.All*Settings.S.D6D.MM*Settings.S.D.All*Settings.S.D.MM.*C.(Settings.Salt).MM;
D.(Settings.Salt).MM = Settings.S.D8D.All*Settings.S.D8D.MM*Settings.S.D.All*Settings.S.D.MM.*D.(Settings.Salt).MM;

C.(Settings.Salt).XX = Settings.S.D6D.All*Settings.S.D6D.XX*Settings.S.D.All*Settings.S.D.XX.*C.(Settings.Salt).XX;
D.(Settings.Salt).XX = Settings.S.D8D.All*Settings.S.D8D.XX*Settings.S.D.All*Settings.S.D.XX.*D.(Settings.Salt).XX;

C.(Settings.Salt).MX = Settings.S.D6D.All*Settings.S.D6D.MX*Settings.S.D.All*Settings.S.D.MX.*C.(Settings.Salt).MX;
D.(Settings.Salt).MX = Settings.S.D8D.All*Settings.S.D8D.MX*Settings.S.D.All*Settings.S.D.MX.*D.(Settings.Salt).MX;

%% If Damping at close range, affects all attractive interactions
for interaction = {'MX' 'XX' 'MM'}
    int = interaction{1};
    if Settings.CR_Damp.(int).r_d >= 0 && Settings.CR_Damp.(int).b >= 0
        r_d = Settings.CR_Damp.(int).r_d;
        sb  = Settings.CR_Damp.(int).b;
        
        f_r.(int) = 1./(1 + exp(-sb.*(r - r_d))); % sigmoid damping function
        df_r.(int) = (sb.*exp(-sb.*(r - r_d)))./((1 + exp(-sb.*(r - r_d))).^2); % sigmoid damping function derivative
        f_cutoff.(int) = 1/(1 + exp(-sb*(Settings.(MDP).RVDW_Cutoff - r_d))); % damping function value at vdw cutoff
    else
        f_r.(int) = 1; % No damping
        df_r.(int) = 0; % Damping derivative is zero
        f_cutoff.(int) = 1; % damping function value at vdw cutoff
    end

    switch int
        case 'MX'
            Y1 = Metal;
            Y2 = Halide;
        case 'MM'
            Y1 = Metal;
            Y2 = Metal;
        case 'XX'
            Y1 = Halide;
            Y2 = Halide;
    end
    
    %% No Damping
    if Settings.C6_Damp.(int) == 0
        % No-damping function
        f6.(int) = f_r.(int);
        f8.(int) = f_r.(int);

        % Derivative parts of no-damping function   
        df6.(int) = df_r.(int);
        df8.(int) = df_r.(int);

        % Values at vdw cutoff
        f6_cutoff.(int) = f_cutoff.(int);
        f8_cutoff.(int) = f_cutoff.(int);

    %% BJ-Type Rational Damping
    elseif Settings.C6_Damp.(int) == 1

        % Damping distances
        if Settings.C6_Damp.input_rvdw
            R0.(int) = Settings.C6_Damp.rvdw.(Y1) + Settings.C6_Damp.rvdw.(Y2);
        else
            R0.(int) = sqrt(D.(Settings.Salt).(int)/C.(Settings.Salt).(int)); % in nm
        end

        % Damping functions (unitless)
        f6.(int) = f_r.(int)./( 1 + ( R0.(int) ./ r ).^6 );

        f8.(int) = f_r.(int)./( 1 + ( R0.(int) ./ r ).^8 );

        % Values of damping function at vdw cutoff
        f6_cutoff.(int) = f_cutoff.(int)/( 1 + ( R0.(int) / Settings.(MDP).RVDW_Cutoff )^6 );
        f8_cutoff.(int) = f_cutoff.(int)/( 1 + ( R0.(int) / Settings.(MDP).RVDW_Cutoff )^8 );

        % Derivative of damping functions
        df6.(int) = f_r.(int).*6.*(R0.(int).^6).*(r.^5)./(((r.^6) + (R0.(int).^6)).^2) + df_r.(int)./( 1 + ( R0.(int) ./ r ).^6 );
        df8.(int) = f_r.(int).*8.*(R0.(int).^8).*(r.^7)./(((r.^8) + (R0.(int).^8)).^2) + df_r.(int)./( 1 + ( R0.(int) ./ r ).^8 );

    %% Tang and Toennies Damping function. Cite:
    % "An improved simple model for the van der Waals potential based on universal damping functions for the dispersion coefficients."
    % K. T. Tang, J. P. Toennies
    % J. Chem. Phys. 1984, 80, 3726-3741.
    elseif Settings.C6_Damp.(int) == 2

        % C6 damping functions
        f6sum = 0;
        for k = 0:6
            f6sum = f6sum + ((alpha.(int).*r).^k)./factorial(k); 
        end
        f6.(int) = f_r.(int).*(1 - f6sum.*exp(-alpha.(int).*r));

        % C8 damping functions
        f8sum = 0;
        for k = 0:8
            f8sum = f8sum + ((alpha.(int).*r).^k)./factorial(k);
        end
        f8.(int) = f_r.(int).*(1 - f8sum.*exp(-alpha.(int).*r));

        % Values of damping function at vdw cutoff
        f6_cutoff.(int) = f_cutoff.(int)*(1 - sum(((alpha.(int).*Settings.(MDP).RVDW_Cutoff).^(0:6))./factorial(0:6)).*exp(-alpha.(int).*Settings.(MDP).RVDW_Cutoff));
        f8_cutoff.(int) = f_cutoff.(int)*(1 - sum(((alpha.(int).*Settings.(MDP).RVDW_Cutoff).^(0:8))./factorial(0:8)).*exp(-alpha.(int).*Settings.(MDP).RVDW_Cutoff));

        % Calculate C6 damping derivatives
        df6sum = 0;
        for k = 1:6
            df6sum = df6sum + k.*(alpha.(int).^k).*(r.^(k-1))./factorial(k);
        end
        df6.(int) = f_r.(int).*((alpha.(int).*exp(-alpha.(int).*r)).*f6sum ...
            - (exp(-alpha.(int).*r)).*df6sum) + df_r.(int).*(1 - f6sum.*exp(-alpha.(int).*r));

        %% Calculate C8 dispersion derivative with damping
        df8sum = 0;
        for k = 1:8
            df8sum = df8sum + k.*(alpha.(int).^k).*(r.^(k-1))./factorial(k);
        end
        df8.(int) = f_r.(int).*((alpha.(int).*exp(-alpha.(int).*r)).*f8sum ...
            - (exp(-alpha.(int).*r)).*df8sum) + df_r.(int).*(1 - f8sum.*exp(-alpha.(int).*r));

    %% Family of other damping functions related by a single formula, requiring van der waals radii
    elseif Settings.C6_Damp.(int) >= 3 && Settings.C6_Damp.(int) <= 6

        % Define crystal radii source:
        % https://en.wikipedia.org/wiki/Ionic_radius
        R0.Li = 0.09; % nm
        R0.Na = 0.116; % nm
        R0.K  = 0.152; % nm
        R0.Rb = 0.166; % nm
        R0.Cs = 0.181; % nm
        R0.F  = 0.119; % nm
        R0.Cl = 0.167; % nm
        R0.Br = 0.182; % nm
        R0.I  = 0.206; % nm

        %% MMDRE Damping function. Citation:
        % "Transferable ab initio intermolecular potentials. 1. Derivation from methanol dimer and trimer calculations"
        % W.T.M. Mooij, F.B. van Duijneveldt, J.G.C.M. van Duijneveldt-van de Rijdt, B.P. van Eijck
        % J. Phys. Chem. A 1999, 103, 9872-9882.
        if Settings.C6_Damp.(int) == 3
            a = -1;
            b = 7.19;
            m = 3;
            n = 2;

        %% PAMoC Damping function. Citation:
        % "Empirical correction to density functional theory for van der Waals interactions"
        % Q. Wu, W. Yang
        % J. Chem. Phys. 2002, 116, 515-524.
        elseif Settings.C6_Damp.(int) == 4
            a = -1;
            b = 3.54;
            m = 3;
            n = 2;

        %% EHFSK Damping function. Citation:
        % "Hydrogen bonding and stacking interactions of nucleic acid base pairs: A density-functional-theory based treatment"
        % M. Elstner, P. Hobza, T. Frauenheim, S. Suhai, E. Kaxiras
        % J. Chem. Phys. 2001, 114, 5149-5155.
        elseif Settings.C6_Damp.(int) == 5
            a = -1;
            b = 3;
            m = 7;
            n = 4;

        %% WY damping function. Citation:
        % "Empirical correction to density functional theory for van der Waals interactions"
        % Q. Wu, W. Yang
        % J. Chem. Phys. 2002, 116, 515-524.
        elseif Settings.C6_Damp.(int) == 6
            a = exp(23);
            b = 23;
            m = 1;
            n = -1;
        end

        R0.(int) = R0.(Y1) + R0.(Y2);

        % Damping functions
        f6.(int) = f_r.(int).*((1 + a.*exp(-b.*(r./R0.(int)).^m) ).^n);
        f8.(int) = f6.(int); % These are not separately defined

        % Derivative parts of damping functions
        df6.(int) = f_r.(int).*(-(a.*b.*m.*n.*(r./R0.(int)).^(m - 1).*exp(-b.*(r./R0.(int)).^m).*(a.*exp(-b.*(r/R0.(int)).^m) + 1).^(n - 1))./R0.(int)) ...
            + df_r.(int).*((1 + a.*exp(-b.*(r./R0.(int)).^m) ).^n);
        df8.(int) = df6.(int);

        % Values at vdw cutoff
        f6_cutoff.(int) = f_cutoff.(int).*((1 + a*exp(-b*(Settings.(MDP).RVDW_Cutoff/R0.(int))^m) )^n);
        f8_cutoff.(int) = f6_cutoff.(int);
    end
    
    % Apply Damping for extra functions
    %% No Damping
    if Settings.C6_Damp.N.(int) == 0
        % No-damping function
        fN.(int) = f_r.(int);

        % Derivative parts of no-damping function   
        dfN.(int) = df_r.(int);
        
        % Values at vdw cutoff
        fN_cutoff.(int) = f_cutoff.(int);

    %% BJ-Type Rational Damping
    elseif Settings.C6_Damp.N.(int) == 1

        % Generate conversion factors
        Bohr_nm = 0.0529177; % a_0 - > Angstrom
        c6conv = 1e-3/2625.4999/((0.052917726)^6); % J/mol nm^6 - > au (from sourcecode)
        J_kJ = 1e-3; % J - > kJ
        Ha_kJmol = 2625.4999; % Ha - > kJ/mol
        c6units = (1/c6conv)*J_kJ; % au - > kJ/mol nm^6
        c8units = (Ha_kJmol)*(Bohr_nm^8); % au - > kJ/mol nm^8

        % Factor used to calculate C8
        sqrt_Q.Li = 5.019869340000000;
        sqrt_Q.Na = 6.585855360000000;
        sqrt_Q.K  = 7.977627530000000;
        sqrt_Q.Rb = 9.554616980000000;
        sqrt_Q.Cs = 11.02204549000000;
        sqrt_Q.F  = 2.388252500000000;
        sqrt_Q.Cl = 3.729323560000000;
        sqrt_Q.Br = 4.590896470000000;
        sqrt_Q.I  = 5.533218150000000;

        % Calculate C8 (needed for cutoff radius)
        C8.(int) = 3.0*(C.(int)/c6units)*sqrt_Q.(Y1)*sqrt_Q.(Y2)*c8units; % in kJ/mol nm^8

        % Damping distances (no C8 term so define wrt crystal radii)
        if Settings.C6_Damp.input_rvdw
            R0.(int) = Settings.C6_Damp.rvdw.(Y1) + Settings.C6_Damp.rvdw.(Y2);
        else
            R0.(int) = sqrt(C8.(int)/C.(int)); % in nm
        end

        % Damping functions (unitless)
        NV = Settings.S.N.(int).Value;
        fN.(int) = f_r.(int)./( 1 + ( R0.(int) ./ r ).^NV );

        % Values of damping function at vdw cutoff
        fN_cutoff.(int) = f_cutoff.(int)/( 1 + ( R0.(int) / Settings.(MDP).RVDW_Cutoff )^NV );

        % Derivative of damping functions
        dfN.(int) = f_r.(int).*(NV.*(R0.(int).^NV).*(r.^(NV-1))./(((r.^NV) + (R0.(int).^NV)).^2)) + df_r.(int)./( 1 + ( R0.(int) ./ r ).^NV );

    %% Tang and Toennies Damping function. Cite:
    % "An improved simple model for the van der Waals potential based on universal damping functions for the dispersion coefficients."
    % K. T. Tang, J. P. Toennies
    % J. Chem. Phys. 1984, 80, 3726-3741.
    elseif Settings.C6_Damp.N.(int) == 2
        
        % CN damping functions
        NV = Settings.S.N.(int).Value;
        
        fNsum = 0;
        for k = 0:NV
            fNsum = fNsum + ((alpha.(int).*r).^k)./factorial(k); 
        end
        fN.(int) = f_r.(int).*(1 - fNsum.*exp(-alpha.(int).*r));

        % Calculate C6 damping derivatives
        dfNsum = 0;
        for k = 1:NV
            dfNsum = dfNsum + k.*(alpha.(int).^k).*(r.^(k-1))./factorial(k);
        end
        dfN.(int) = f_r.(int).*((alpha.(int).*exp(-alpha.(int).*r)).*fNsum ...
            - (exp(-alpha.(int).*r)).*dfNsum) + df_r.(int).*(1 - fNsum.*exp(-alpha.(int).*r));

        % Values of damping function at vdw cutoff
        fN_cutoff.(int) = f_cutoff.(int)*(1 - sum(((alpha.(int).*Settings.(MDP).RVDW_Cutoff).^(0:NV))./factorial(0:NV)).*exp(-alpha.(int).*Settings.(MDP).RVDW_Cutoff));

    %% Family of other damping functions related by a single formula, requiring van der waals radii
    elseif Settings.C6_Damp.N.(int) >= 3 && Settings.C6_Damp.N.(int) <= 6

        % Define crystal radii source:
        % https://en.wikipedia.org/wiki/Ionic_radius
        R0.Li = 0.09; % nm
        R0.Na = 0.116; % nm
        R0.K  = 0.152; % nm
        R0.Rb = 0.166; % nm
        R0.Cs = 0.181; % nm
        R0.F  = 0.119; % nm
        R0.Cl = 0.167; % nm
        R0.Br = 0.182; % nm
        R0.I  = 0.206; % nm

        %% MMDRE Damping function. Citation:
        % "Transferable ab initio intermolecular potentials. 1. Derivation from methanol dimer and trimer calculations"
        % W.T.M. Mooij, F.B. van Duijneveldt, J.G.C.M. van Duijneveldt-van de Rijdt, B.P. van Eijck
        % J. Phys. Chem. A 1999, 103, 9872-9882.
        if Settings.C6_Damp.N.(int) == 3
            a = -1;
            b = 7.19;
            m = 3;
            n = 2;

        %% PAMoC Damping function. Citation:
        % "Empirical correction to density functional theory for van der Waals interactions"
        % Q. Wu, W. Yang
        % J. Chem. Phys. 2002, 116, 515-524.
        elseif Settings.C6_Damp.N.(int) == 4
            a = -1;
            b = 3.54;
            m = 3;
            n = 2;

        %% EHFSK Damping function. Citation:
        % "Hydrogen bonding and stacking interactions of nucleic acid base pairs: A density-functional-theory based treatment"
        % M. Elstner, P. Hobza, T. Frauenheim, S. Suhai, E. Kaxiras
        % J. Chem. Phys. 2001, 114, 5149-5155.
        elseif Settings.C6_Damp.N.(int) == 5
            a = -1;
            b = 3;
            m = 7;
            n = 4;

        %% WY damping function. Citation:
        % "Empirical correction to density functional theory for van der Waals interactions"
        % Q. Wu, W. Yang
        % J. Chem. Phys. 2002, 116, 515-524.
        elseif Settings.C6_Damp.N.(int) == 6
            a = exp(23);
            b = 23;
            m = 1;
            n = -1;
        end

        R0.(int) = R0.(Y1) + R0.(Y2);

        % Damping functions
        fN.(int) = f_r.(int).*((1 + a.*exp(-b.*(r./R0.(int)).^m) ).^n);

        % Derivative parts of damping functions
        dfN.(int) = f_r.(int).*(-(a.*b.*m.*n.*(r./R0.(int)).^(m - 1).*exp(-b.*(r./R0.(int)).^m).*(a.*exp(-b.*(r/R0.(int)).^m) + 1).^(n - 1))./R0.(int)) ...
            + df_r.(int).*((1 + a.*exp(-b.*(r./R0.(int)).^m) ).^n);

        % Values at vdw cutoff
        fN_cutoff.(int) = f_cutoff.(int)*((1 + a*exp(-b*(Settings.(MDP).RVDW_Cutoff/R0.(int))^m) )^n);
    end
    
    %% Modify potential with Gaussian Adjustments
    G_r.(int) = zeros(1,length(r));
    dG_r.(int) = zeros(1,length(r));
    G_r_cutoff.(int) = 0;
    for i = 1:length(G_a.(int))
        G_r.(int) = G_r.(int) + G_a.(int)(i).*exp((-(r - G_b.(int)(i)).^2)./(2.*(G_c.(int)(i).^2)));
        G_r_cutoff.(int) = G_r_cutoff.(int) + G_a.(int)(i)*exp((-(Settings.(MDP).RVDW_Cutoff - G_b.(int)(i))^2)/(2*(G_c.(int)(i)^2)));
        dG_r.(int) = dG_r.(int) - (G_a.(int)(i).*(r - G_b.(int)(i))).*(exp((-(r - G_b.(int)(i)).^2)./(2.*(G_c.(int)(i).^2))))/(G_c.(int)(i).^2);
    end
    
    %% Build PES
    if strcmp(int,'MX')
        % Place a CR damping function for the attractive coulomb potential
        U_Qdamp = - k_0*(e_c^2)*q.(Y1)*q.(Y2)./(r) ...
                  + f_r.(int).*k_0*(e_c^2)*q.(Y1)*q.(Y2)./(r);

        % Negative Derivative of CR damping function
        dU_Qdamp = - k_0*(e_c^2)*q.(Y1)*q.(Y2)./(r.^2) ...
                   + f_r.(int).*k_0*(e_c^2)*q.(Y1)*q.(Y2)./(r.^2) ...
                   - df_r.(int).*k_0*(e_c^2)*q.(Y1)*q.(Y2)./(r);
       
        U_Qdamp_cutoff = - k_0*(e_c^2)*q.(Y1)*q.(Y2)/(Settings.(MDP).RVDW_Cutoff) ...
                         + f_cutoff.(int)*k_0*(e_c^2)*q.(Y1)*q.(Y2)/(Settings.(MDP).RVDW_Cutoff);
    else
        % No close-range coulomb damping
        U_Qdamp = 0;
        dU_Qdamp = 0;
        U_Qdamp_cutoff = 0;
    end
    
    % Additional functions
    if Settings.S.N.(int).Value > 0
        N = Settings.S.N.(int).Value; % Exponent
        C = Settings.S.N.(int).Scale; % Scaling
        U_add  = - fN.(int).*C./(r.^N);
        dU_add = + fN.(int).*C.*N./(r.^(N+1)) ...
                 - dfN.(int).*C./(r.^N);
        U_add_cutoff = - fN_cutoff.(int)*C/(Settings.(MDP).RVDW_Cutoff^N);
    else
        U_add = 0;
        dU_add = 0;
        U_add_cutoff = 0;
    end
    
    % Components of potential
    U.(int).f = 1./r; % Electrostatics function f(r)
    U.(int).g = - f6.(int).*C.(Settings.Salt).(int)./(r.^6) ...
                - f8.(int).*D.(Settings.Salt).(int)./(r.^8) ...
                + G_r.(int) ...
                + U_add; % Dispersion g(r)
    U.(int).h = B.(int)*exp(-alpha.(int).*r) ...
                + U_Qdamp; % Short range repulsion h(r) (with possible close-range coulomb damping)

    % Negative components of derivative
    U.(int).df = 1./(r.^2); % Electrostatics function (not including Coulomb constant or charges)
    U.(int).dg = - f6.(int).*6.*C.(Settings.Salt).(int)./(r.^7) ...
                 + df6.(int).*C.(Settings.Salt).(int)./(r.^6) ...
                 - f8.(int).*8.*D.(Settings.Salt).(int)./(r.^9) ...
                 + df8.(int).*D.(Settings.Salt).(int)./(r.^8) ...
                 - dG_r.(int) ...
                 + dU_add; % Dispersion -dg(r)/dr
    U.(int).dh = alpha.(int)*B.(int)*exp(-alpha.(int).*r) ...
                 + dU_Qdamp;% Short range repulsion -dh(r)/dr

    if contains(Settings.(MDP).vdw_modifier,'potential-shift','IgnoreCase',true)
        EVDW_Cutoff = B.(int)*exp(-alpha.(int)*Settings.(MDP).RVDW_Cutoff) ...
                    + U_Qdamp_cutoff ...
                    - f6_cutoff.(int)*C.(Settings.Salt).(int)/(Settings.(MDP).RVDW_Cutoff^6) ...
                    - f8_cutoff.(int)*D.(Settings.Salt).(int)/(Settings.(MDP).RVDW_Cutoff^8) ...
                    + G_r_cutoff.(int) ...
                    + U_add_cutoff;

        % Shift by the dispersion energy at vdw cutoff radius. only affects one
        % energy component, not derivatives (i.e. forces)
        U.(int).g = U.(int).g - EVDW_Cutoff;
    end

    % remove infinities
    U.(int) = Remove_Infinities(U.(int));

    if FactorOutParams
        U.(int).g = U.(int).g./C.(Settings.Salt).(int);
        U.(int).dg = U.(int).dg./C.(Settings.Salt).(int);
        U.(int).h = U.(int).h./B.(int);
        U.(int).dh = U.(int).dh./B.(int);
    end
    
    % Print
    U_out = [r ; U.(int).f ; U.(int).df ; U.(int).g ; U.(int).dg ; U.(int).h ; U.(int).dh];
    U.(int).out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],U_out(:)) );
end

if ReturnAsStructure
	U.MX.f0 = U.MX.f;
    U.MM.f0 = U.MM.f;
    U.XX.f0 = U.XX.f;
    
	U.MX.df0 = U.MX.df;
    U.MM.df0 = U.MM.df;
    U.XX.df0 = U.XX.df;
    
	U.MX.f = k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.f;
    U.MM.f = k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.f;
    U.XX.f = k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.f;
    
	U.MX.df = k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.df;
    U.MM.df = k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.df;
    U.XX.df = k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.df;
    
    U.MX.Total = U.MX.f + U.MX.g + U.MX.h;
    U.MM.Total = U.MM.f + U.MM.g + U.MM.h;
    U.XX.Total = U.XX.f + U.XX.g + U.XX.h;
    
    U.MX.dTotal = -(U.MX.df + U.MX.dg + U.MX.dh);
    U.MM.dTotal = -(U.MM.df + U.MM.dg + U.MM.dh);
    U.XX.dTotal = -(U.XX.df + U.XX.dg + U.XX.dh);
    
    U_MX_out = U.MX;
    U_MM_out = U.MM;
    U_XX_out = U.XX;
    
    U_MX_out.r = r;
    U_MM_out.r = r;
    U_XX_out.r = r;
else
    U_MX_out = U.MX.out;
    U_MM_out = U.MM.out;
    U_XX_out = U.XX.out;
end

%% PLOT if plotswitch chosen
if Plotswitch
    figure;
    
    % Options
    lw=2;
    fs=25;

    h = cell(1,8);
    hold on
    switch lower(PlotType)
        case 'full'
            h{1} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.f + U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.f + U.MM.g + U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.f + U.XX.g + U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Full Potential';
        case 'full-derivative'
            h{1} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.df + U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.df + U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.df + U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Derivative of Full Potential';
        case 'lj'
            h{1} = plot(r.*10,U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.g + U.MM.h,'Color','g','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.g + U.XX.h,'Color','b','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Lennard-Jones Potential';
        case 'lj-derivative'
            h{1} = plot(r.*10,U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Lennard-Jones Potential';
        case 'dispersion'
            h{1} = plot(r.*10,U.MX.g,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.g,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.g,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Dispersion Potential';
        case 'dispersion-derivative'
            h{1} = plot(r.*10,U.MX.dg,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dg,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dg,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Dispersion Potential';
        case 'repulsive'
            h{1} = plot(r.*10,U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Repulsive Potential';
        case 'repulsive-derivative'
            h{1} = plot(r.*10,U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Repulsive Potential';
    end
    
    title(['Plot of ' ttxt ' for ' Settings.Salt ' TF Model'],...
       'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('r (Angstrom)','fontsize',fs,'Interpreter','latex');
    ylabel('Pair Potential (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');

    %ylim([-1000 500]);
    ylim(yl);
    xlim([Startpoint Settings.Table_Length]);

    % Blank line
    hline = refline([0 0]);
    hline.Color = 'k';
    hline.LineWidth = lw-1;
    hline.LineStyle = '--';
    leg1 = legend([h{:}],{[Metal '$^{+}$' ' - ' Halide '$^{-}$'] [Metal '$^{+}$' ' - ' Metal '$^{+}$'] [Halide '$^{-}$' ' - ' Halide '$^{-}$']});
    leg1.Interpreter = 'latex';
end

end