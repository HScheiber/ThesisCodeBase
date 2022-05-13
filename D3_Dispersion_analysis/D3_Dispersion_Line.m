%% Inputs
Salt = 'LiI';
XC_Func = 'PW1PW'; %  only for D3 Dispersion
Theories = {'JC'};% 'JC' 'TF'
%% Scaling parameters for TF and JC models
S_D = 1;
S_R = 1;
S_E = 1;
S_S = 1;
S_MMD = 1;
S_XXD = 1;
S_MXD = 1;
S_C = 1; % Scaling charges
S_D6 = 1; % Scaling R6
S_D8 = 1; % Scaling R8
%% Define your own C6 Params, Damping radii, and Scale parameters
Input_C6 = false;
C6_MM = 87.29; % in kJ/mol Ang^6
C6_MX = 716.62; % in kJ/mol Ang^6
C6_XX = 9195.59; % in kJ/mol Ang^6

C8_MM = 198.9; % in kJ/mol Ang^6
C8_MX = 4388.04; % in kJ/mol Ang^6
C8_XX = 87244.23; % in kJ/mol Ang^6

Input_Damp_Radii = false;
R0_MM = 1.16; % XC dependent function (Angstroms);
R0_XX = 1.67; % XC dependent function (Angstroms);
R0_MX = (R0_MM+R0_XX)/2;
Input_Scale_Params = false; % Modify scaling/damping parameters
a1=1;
a2=0;
s6=1;
s8=1;

%%
Structures = {'NiAs' 'Rocksalt'}; % 'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'
Plotswitch = true;
if ispc
    datadir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\';
else
    datadir = '/home/user/ThesisCodeBase/data/';
end
Data_Type = 1;
Basis_Set = 'pob-TZVP';
Stepsize = 0.01;
Rij = 0.01:Stepsize:10.1; %D3 default is 50.1  Bond Lengths in angstroms
fs = 34; % Font size
lw = 1.5; % Line width
Y_limits = [-5 0.5];
X_limits = [0 8.2];
Subtitle_X = 0.90;
Subtitle_Y = 0.1;
Include_C6 = true; % include C6 dispersion term
Include_C8 = true; % include C8 dispersion term (where applicable)
Include_R = false; % Include repulsive wall in JC and TF models
Include_Coulomb = false; % Include simple coulomb term
Text_Overlay = false; % Include a text overlay at each circle
Transparency = 0.8;
Plot3D = true;
Merge_Tol = 0.001; % merge coordination numbers closer than this
Star_Area = 'CN'; % Either CN or E or EP. Scale the area of the circle corresponding to either coordination number or energy of the shell
if strcmp(Star_Area,'CN')
    AC = 50; % Area of plotted circle when CN = 1. This is scaled by CN
elseif strcmp(Star_Area,'EP')
    AC = 10;
elseif strcmp(Star_Area,'E')
    AC = 50;
end
% Set coordination numbers
CN_M = 100;
CN_X = 100;

Panel_Position{1} = [0.10375 0.69 0.7509 0.29]; %[X Y Width Height]
Panel_Position{2} = [0.10375 0.40 0.7509 0.29];
Panel_Position{3} = [0.10375 0.11 0.7509 0.29];
Legend_Position = [0.8650 0.6600 0.1224 0.3053];
N_Struc = length(Structures);
N_Theor = length(Theories);
Colors = cbrewer('qual','Set1',max(N_Struc,3),'PCHIP');
if N_Theor == 1
    Colors_Theories = [0 0 0];
else
    Colors_Theories = cbrewer('qual','Dark2',max(N_Theor,3),'PCHIP');
end
Panel_Title = cell(3,1);
LineStyles = {'-' '-.' '--' '--' '-' '--' '-'}; %'-' (default) | '--' | ':' | '-.' | 'none')
[Metal,Halide] = Separate_Metal_Halide(Salt);
%% Global parameters of D3(BJ) and XC specific parameters
Version = 4; % D3(BJ)
TZ = false;
if ~Input_Scale_Params
    [s6,a1,s8,a2,~] = setfuncpar(XC_Func,Version,TZ);
end
k1 = 16;
k2 = 4/3;
k3 = -4;

%% Conversion factors and Cutoff radii
Bohr_Ang = 0.52917726; % a_0 - > Angstrom
CN_Cutoff = sqrt(1600)*Bohr_Ang; % Cutoff for CN calculation in Angstroms
Disp_Cutoff = sqrt(9000.0)*Bohr_Ang; % Cutoff for dispersion calculation in Angstroms
c6conv = 1e-3/2625.4999/((0.052917726)^6); % J/mol nm^6 - > au (from sourcecode)
nm_Ang = 10; % nm - > Angstrom
J_kJ = 1e-3; % J - > kJ
Ha_kJmol = 2625.4999; % Ha - > kJ/mol
c6units = (1/c6conv)*J_kJ*(nm_Ang^6); % au - > kJ/mol Ang^6
c8units = (Ha_kJmol)*(Bohr_Ang^8); % au - > kJ/mol Ang^8

%% Load figure
if Plotswitch
    figh = figure('WindowState','maximized','NumberTitle','off',...
        'Name','','Visible','On');
    axh = cell(3,1);
    set(figh,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
end
NR = length(Rij);

Scaling_Params(1) = S_D;
Scaling_Params(2) = S_R;
Scaling_Params(3) = S_E;
Scaling_Params(4) = S_S;
Scaling_Params(5) = S_MMD;
Scaling_Params(6) = S_XXD;
Scaling_Params(7) = S_MXD;
Scaling_Params(8) = S_C;
Scaling_Params(9) = S_D6;
Scaling_Params(10) = S_D8;
if ~Include_C6
    Scaling_Params(9) = 0;
end
if ~Include_C8
    Scaling_Params(10) = 0;
end
if ~Include_R
    Scaling_Params(2) = 0;
end
if ~Include_Coulomb
    Scaling_Params(8) = 0;
end

%% Loop through theories
n=0; % counter for legend
pobj = gobjects(N_Struc,1); % legend objects for circles
linobj = gobjects(N_Theor,1); % legend objects for lines
for tidx = 1:length(Theories)
    Theory = Theories{tidx};
    LineCol = Colors_Theories(tidx,:);

    if strcmp(Theory,'D3(BJ)')
        %% Load C6 function,max coordination number function (mxc) vs atomic number, and r2r4 from disc
        [C6_Func,mxc] = copyc6;
        loadr2r4 = load('r2r4.mat','r2r4');
        r2r4 = loadr2r4.r2r4;

        %% Atomic numbers of atoms in salt
        Z_M = elements('Symbol',Metal,'atomic_number');
        Z_X = elements('Symbol',Halide,'atomic_number');

        %% Calculate C6 and C8 for each possible pair of atoms in unit cell, as well as R0AB
        sqrt_Q_M = r2r4(Z_M); % Factor used to calculate C8 for Metal
        sqrt_Q_X = r2r4(Z_X); % Factor used to calculate C8 for Halide

        % C6 coefficients
        if ~Input_C6
            C6_MX = getc6(C6_Func,mxc,Z_M,Z_X,CN_M,CN_X)*c6units; % in kJ/mol Ang^6
            C6_MM = getc6(C6_Func,mxc,Z_M,Z_M,CN_M,CN_M)*c6units; % in kJ/mol Ang^6
            C6_XX = getc6(C6_Func,mxc,Z_X,Z_X,CN_X,CN_X)*c6units; % in kJ/mol Ang^6
        end
        
        disp(['C6 Coefficients when ' Metal ' FCN = ' num2str(CN_M,'%2.2f') ' and ' Halide ' FCN = ' num2str(CN_X,'%2.2f')])
        disp([Metal '-' Halide ' C6 = ' num2str(C6_MX,'%2.10f') ' kj/mol-Ang^6']);
        disp([Metal '-' Metal ' C6 = ' num2str(C6_MM,'%2.10f') ' kj/mol-Ang^6']);
        disp([Halide '-' Halide ' C6 = ' num2str(C6_XX,'%2.10f') ' kj/mol-Ang^6']);
        
        % C8 coefficients
        if ~Input_C6
            C8_MX = 3.0*(C6_MX/c6units)*sqrt_Q_M*sqrt_Q_X*c8units; % in kJ/mol Ang^8
            C8_MM = 3.0*(C6_MM/c6units)*sqrt_Q_M*sqrt_Q_M*c8units; % in kJ/mol Ang^8
            C8_XX = 3.0*(C6_XX/c6units)*sqrt_Q_X*sqrt_Q_X*c8units; % in kJ/mol Ang^8
        end
        
        disp(['C8 Coefficients when ' Metal ' FCN = ' num2str(CN_M,'%2.2f') ' and ' Halide ' FCN = ' num2str(CN_X,'%2.2f')])
        disp([Metal '-' Halide ' C8 = ' num2str(C8_MX,'%2.10f') ' kj/mol-Ang^8']);
        disp([Metal '-' Metal ' C8 = ' num2str(C8_MM,'%2.10f') ' kj/mol-Ang^8']);
        disp([Halide '-' Halide ' C8 = ' num2str(C8_XX,'%2.10f') ' kj/mol-Ang^8']);
        
        
        % vdW radius for damping
        if ~Input_Damp_Radii
            R0_MX = sqrt(C8_MX/C6_MX); % in Angstroms
            R0_MM = sqrt(C8_MM/C6_MM); % in Angstroms  
            R0_XX = sqrt(C8_XX/C6_XX); % in Angstroms
        end
        
        disp('Default D3(BJ) vdW Damping Radii:')
        disp([Metal '-' Halide ' R0 = ' num2str(R0_MX,'%2.10f') ' Angstrom']);
        disp([Metal '-' Metal ' R0 = ' num2str(R0_MM,'%2.10f') ' Angstrom']);
        disp([Halide '-' Halide ' R0 = ' num2str(R0_XX,'%2.10f') ' Angstrom']);

        % Modify damping radius based on theory
        F_R0_MX = a1*R0_MX + a2*Bohr_Ang; % XC dependent function (Angstroms);
        F_R0_MM = a1*R0_MM + a2*Bohr_Ang; % XC dependent function (Angstroms);
        F_R0_XX = a1*R0_XX + a2*Bohr_Ang; % XC dependent function (Angstroms);

        disp([XC_Func '- modified vdW Damping Radii:'])
        disp([Metal '-' Halide ' R0 = ' num2str(F_R0_MX,'%2.10f') ' Angstrom']);
        disp([Metal '-' Metal ' R0 = ' num2str(F_R0_MM,'%2.10f') ' Angstrom']);
        disp([Halide '-' Halide ' R0 = ' num2str(F_R0_XX,'%2.10f') ' Angstrom']);
        
        %% Calculate D3(BJ) pair dispersion energy vs BL   
        if Include_C6
            E6_MX = -s6.*C6_MX./((Rij.^6) + (F_R0_MX.^6));
            E6_MM = -s6.*C6_MM./((Rij.^6) + (F_R0_MM.^6));
            E6_XX = -s6.*C6_XX./((Rij.^6) + (F_R0_XX.^6));
        else
            E6_MX = zeros(1,NR);
            E6_MM = zeros(1,NR);
            E6_XX = zeros(1,NR);
        end
        if Include_C8
            E8_MX = -s8.*C8_MX./((Rij.^8) + (F_R0_MX.^8));
            E8_MM = -s8.*C8_MM./((Rij.^8) + (F_R0_MM.^8));
            E8_XX = -s8.*C8_XX./((Rij.^8) + (F_R0_XX.^8));
        else
            E8_MX = zeros(1,NR);
            E8_MM = zeros(1,NR);
            E8_XX = zeros(1,NR);
        end

        E_Tot_MX = E6_MX + E8_MX; % in kJ/mol
        E_Tot_MM = E6_MM + E8_MM; % in kJ/mol
        E_Tot_XX = E6_XX + E8_XX; % in kJ/mol
    elseif strcmp(Theory,'D4(BJ)')
        Structure = Structures{1};
        
        % Load data from file
        load(fullfile(datadir,'CRYSTAL',[Salt '_' Structure '_Lattice_Energies.mat']),'Data');

        Label = [Salt '_' Structure(1) '_' XC_Func '_D4'];
        SubData = Data.(Label);

        DOI = [];
        for i = 1:size(SubData,1)
            if strcmp(SubData{i,1},Basis_Set)
                DOI = SubData{i,2};
                break
            end
        end
        if ~isempty(DOI)
            DT = DOI{10};

            DatInd = ismember(DT,Data_Type);

            if sum(DatInd) == 0
                continue
            end

            CryStruc.a = DOI{2}(DatInd);
            CryStruc.b = DOI{3}(DatInd);
            CryStruc.c = DOI{4}(DatInd);
            CryStruc.alpha = DOI{11}(DatInd);
            CryStruc.beta = DOI{12}(DatInd);
            CryStruc.gamma = DOI{13}(DatInd);
            CryStruc.FC = DOI{7}(DatInd);        
            if sum(DatInd) > 1
                En = DOI{8}(DatInd);
                [~,midx] = min(En);
                CryStruc.a = CryStruc.a(midx);
                CryStruc.b = CryStruc.b(midx);
                CryStruc.c = CryStruc.c(midx);
                CryStruc.alpha = CryStruc.alpha(midx);
                CryStruc.beta = CryStruc.beta(midx);
                CryStruc.gamma = CryStruc.gamma(midx);
                CryStruc.FC = CryStruc.FC{midx};
            else
                CryStruc.FC = CryStruc.FC{1};
            end
        else
            continue
        end
        
        switch lower(Structure)
            case 'rocksalt'
                CryStruc.N = 8;
            case 'wurtzite'
                CryStruc.N = 4;
            case 'cscl'
                CryStruc.N = 2;
            case 'fivefive'
                CryStruc.N = 8;
            case 'nias'
                CryStruc.N = 4;
            case 'sphalerite'
                CryStruc.N = 8;
            case 'betabeo'
                CryStruc.N = 8;
        end
        
        % Get FC of entire unit cell
        [CryStruc.FC_Metal,CryStruc.FC_Halide] = UnitCell_FractionalCoords([CryStruc.FC{1,2:end}],...
        	[CryStruc.FC{2,2:end}],Structure);
        CryStruc.Transform = GenTransformMatrix(CryStruc);

        % Create vasp format geometry file in given Directory
        N = CryStruc.N;   
        geom_txt = ['New Structure' newline '1.0' newline];
        TM = CryStruc.Transform*[CryStruc.a 0 0; 0 CryStruc.b 0; 0 0 CryStruc.c];
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
                geom_txt = [geom_txt pad(num2str(CryStruc.FC_Metal(i,j),'%10.10f'),16,'left')];
            end
            geom_txt = [geom_txt newline];
        end
        for i = 1:(N/2)
            for j = 1:3
                geom_txt = [geom_txt pad(num2str(CryStruc.FC_Halide(i,j),'%10.10f'),16,'left')]; %#ok<*AGROW>
            end
            geom_txt = [geom_txt newline];
        end

        % save to file
        tmpdir = tempname;
        mkdir(tmpdir);
        filename = [tmpdir filesep 'tempstruct.vasp'];
        fidPM = fopen(filename,'wt');
        fwrite(fidPM,regexprep(geom_txt,{'\r', '\n\n+'}',{'', '\n'}));
        fclose(fidPM);

        % run dftd4 on the structure
        if ispc % for testing
            dftd4 = 'wsl source ~/.bashrc; dftd4 ';
            dftd4_func = ['wsl source ~/.bashrc; dftd4 -f ' XC_Func ' '];
            fnunix = windows2unix(filename);
        elseif isunix
            dftd4 = 'dftd4_24 ';
            dftd4_func = ['dftd4_24 -f ' XC_Func ' '];
            fnunix = filename;
        end
        [ercode,dftd4_out] = system([dftd4 fnunix]);

        if ercode ~= 0 %dftd4 failed
            error(['DFTD4 module failed from input: '  dftd4 fnunix])
        end
        
        [ercode,dftd4_func_out] = system([dftd4_func fnunix]);

        if ercode ~= 0 %dftd4 failed
            error(['DFTD4 module failed from input: '  dftd4_func fnunix])
        end

        % Grab C6's from output (in AU)
        C6s = regexp(dftd4_out,'# +Z +covCN +q +C6AA.+\n\n\n','match','ONCE');
        MM_C6 = regexp(C6s,[Metal ' +([-.0-9]+) +[-.0-9]+ +([-.0-9]+) +'],'tokens','once');
        CN_M = str2double(MM_C6{1});
        C6_MM = str2double(MM_C6{2})*c6units;
        XX_C6 = regexp(C6s,[Halide ' +([-.0-9]+) +[-.0-9]+ +([-.0-9]+) +'],'tokens','once');
        CN_X = str2double(XX_C6{1});
        C6_XX = str2double(XX_C6{2})*c6units;

        Mol_C6 = regexp(C6s,'Mol\. C6AA.+? + : +([-.0-9]+)','tokens','once');
        C6_Mol = str2double(Mol_C6{1});

        % Calculate cross term C6
        C6_MX = (C6_Mol*c6units - ((N/2)^2)*C6_MM - ((N/2)^2)*C6_XX)/(2*(N/2)^2);

        % Grab damping parameters from output
        Dampings = regexp(dftd4_func_out,'Damping.+?Results','match','ONCE');
        s6 = regexp(Dampings,'s6 +: +([-.0-9]+)','tokens','once');
        s6 = str2double(s6{1});
        s8 = regexp(Dampings,'s8 +: +([-.0-9]+)','tokens','once');
        s8 = str2double(s8{1});
        a1 = regexp(Dampings,'a1 +: +([-.0-9]+)','tokens','once');
        a1 = str2double(a1{1});
        a2 = regexp(Dampings,'a2 +: +([-.0-9]+)','tokens','once');
        a2 = str2double(a2{1});
        
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
        C.(Salt).MX = C6_MX*c6units; % in kJ/mol Ang^6
        C.(Salt).MM = C6_MM*c6units; % in kJ/mol Ang^6
        C.(Salt).XX = C6_XX*c6units; % in kJ/mol Ang^6

        % C8 coefficients
        D.(Salt).MX = 3.0*(C6_MX)*sqrt_Q_M*sqrt_Q_X*c8units; % in kJ/mol Ang^8
        D.(Salt).MM = 3.0*(C6_MM)*sqrt_Q_M*sqrt_Q_M*c8units; % in kJ/mol Ang^8
        D.(Salt).XX = 3.0*(C6_XX)*sqrt_Q_X*sqrt_Q_X*c8units; % in kJ/mol Ang^8 
        
        % Output
        disp(['C6 Coefficients when ' Metal ' FCN = ' num2str(CN_M,'%2.2f') ' and ' Halide ' FCN = ' num2str(CN_X,'%2.2f')])
        disp([Metal '-' Halide ' C6 = ' num2str(C6_MX,'%2.10f') ' kj/mol-Ang^6']);
        disp([Metal '-' Metal ' C6 = ' num2str(C6_MM,'%2.10f') ' kj/mol-Ang^6']);
        disp([Halide '-' Halide ' C6 = ' num2str(C6_XX,'%2.10f') ' kj/mol-Ang^6']);
        
        % C8 coefficients
        if ~Input_C6
            C8_MX = 3.0*(C6_MX/c6units)*sqrt_Q_M*sqrt_Q_X*c8units; % in kJ/mol Ang^8
            C8_MM = 3.0*(C6_MM/c6units)*sqrt_Q_M*sqrt_Q_M*c8units; % in kJ/mol Ang^8
            C8_XX = 3.0*(C6_XX/c6units)*sqrt_Q_X*sqrt_Q_X*c8units; % in kJ/mol Ang^8
        end
        
        disp(['C8 Coefficients when ' Metal ' FCN = ' num2str(CN_M,'%2.2f') ' and ' Halide ' FCN = ' num2str(CN_X,'%2.2f')])
        disp([Metal '-' Halide ' C8 = ' num2str(C8_MX,'%2.10f') ' kj/mol-Ang^8']);
        disp([Metal '-' Metal ' C8 = ' num2str(C8_MM,'%2.10f') ' kj/mol-Ang^8']);
        disp([Halide '-' Halide ' C8 = ' num2str(C8_XX,'%2.10f') ' kj/mol-Ang^8']);
        
        
        % vdW radius for damping
        if ~Input_Damp_Radii
            R0_MX = sqrt(C8_MX/C6_MX); % in Angstroms
            R0_MM = sqrt(C8_MM/C6_MM); % in Angstroms  
            R0_XX = sqrt(C8_XX/C6_XX); % in Angstroms
        end
        
        disp('Default D4(BJ) vdW Damping Radii:')
        disp([Metal '-' Halide ' R0 = ' num2str(R0_MX,'%2.10f') ' Angstrom']);
        disp([Metal '-' Metal ' R0 = ' num2str(R0_MM,'%2.10f') ' Angstrom']);
        disp([Halide '-' Halide ' R0 = ' num2str(R0_XX,'%2.10f') ' Angstrom']);

        % Modify damping radius based on theory
        F_R0_MX = a1*R0_MX + a2*Bohr_Ang; % XC dependent function (Angstroms);
        F_R0_MM = a1*R0_MM + a2*Bohr_Ang; % XC dependent function (Angstroms);
        F_R0_XX = a1*R0_XX + a2*Bohr_Ang; % XC dependent function (Angstroms);

        disp([XC_Func '- modified vdW Damping Radii:'])
        disp([Metal '-' Halide ' R0 = ' num2str(F_R0_MX,'%2.10f') ' Angstrom']);
        disp([Metal '-' Metal ' R0 = ' num2str(F_R0_MM,'%2.10f') ' Angstrom']);
        disp([Halide '-' Halide ' R0 = ' num2str(F_R0_XX,'%2.10f') ' Angstrom']);
        
        %% Calculate D4(BJ) pair dispersion energy vs BL   
        if Include_C6
            E6_MX = -s6.*C6_MX./((Rij.^6) + (F_R0_MX.^6));
            E6_MM = -s6.*C6_MM./((Rij.^6) + (F_R0_MM.^6));
            E6_XX = -s6.*C6_XX./((Rij.^6) + (F_R0_XX.^6));
        else
            E6_MX = zeros(1,NR);
            E6_MM = zeros(1,NR);
            E6_XX = zeros(1,NR);
        end
        if Include_C8
            E8_MX = -s8.*C8_MX./((Rij.^8) + (F_R0_MX.^8));
            E8_MM = -s8.*C8_MM./((Rij.^8) + (F_R0_MM.^8));
            E8_XX = -s8.*C8_XX./((Rij.^8) + (F_R0_XX.^8));
        else
            E8_MX = zeros(1,NR);
            E8_MM = zeros(1,NR);
            E8_XX = zeros(1,NR);
        end

        E_Tot_MX = E6_MX + E8_MX; % in kJ/mol
        E_Tot_MM = E6_MM + E8_MM; % in kJ/mol
        E_Tot_XX = E6_XX + E8_XX; % in kJ/mol
        
        
    elseif strcmp(Theory,'TF')
        TF_Model = TF_Pair_Potential(Rij(1),Rij(end),Stepsize,Salt,...
            Scaling_Params,'full'); %#ok<*UNRCH>
        E_Tot_MX = TF_Model.MX; % in kJ/mol
        E_Tot_MM = TF_Model.MM; % in kJ/mol
        E_Tot_XX = TF_Model.XX; % in kJ/mol
        
    elseif contains(Theory,'JC')
        switch Theory
            case 'JC'
                Watermodel = 'SPC/E';
            case 'JC3P'
                Watermodel = 'TIP3P';
            case 'JC4P'
                Watermodel = 'TIP4PEW';
            otherwise
                error(['Unknown JC model: ' Theory]);
        end
        
        JC_Model = JC_Pair_Potential(Rij(1),Rij(end),Stepsize,Salt,...
            Watermodel,Scaling_Params,'full');

        
        E_Tot_MX = JC_Model.MX; % in kJ/mol
        E_Tot_MM = JC_Model.MM; % in kJ/mol
        E_Tot_XX = JC_Model.XX; % in kJ/mol
    end

    % Remove infinity
    E_Tot_MX(isinf(E_Tot_MX)) = -1e10;
    E_Tot_MM(isinf(E_Tot_MM)) = -1e10;
    E_Tot_XX(isinf(E_Tot_XX)) = -1e10;
    
    %% Plot the dispersion lines
    if Plotswitch
        for ii = 1:3 % Loop over the three plots
            if ii == 1 % M-X
                E = E_Tot_MX;
                Panel_Title = [Metal '$^{+}$' Halide '$^{-}$'];
            elseif ii == 2 % M-M
                E = E_Tot_MM;
                Panel_Title = [Metal '$^{+}$' Metal '$^{+}$'];
            elseif ii == 3
                E = E_Tot_XX; % X-X
                Panel_Title = [Halide '$^{-}$' Halide '$^{-}$'];
            end
            % Generate axis
            if tidx == 1
                axh{ii} = axes('Parent',figh,'Position',Panel_Position{ii},'FontSize',fs,...
                    'box','on','TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on');
                hold(axh{ii},'on')

                % Plot zero line
                yline(axh{ii},0,'Color','k','LineStyle','--','LineWidth',lw);
            end

            % plot energy curve
            if ii == 1
                linobj(tidx) = plot(axh{ii},Rij,E,'Color',LineCol,'LineStyle','-','LineWidth',lw+2,'Marker','none');
            else
                plot(axh{ii},Rij,E,'Color',LineCol,'LineStyle','-','LineWidth',lw+2,'Marker','none');
            end

            if tidx == 1
                % Put in axis sub-title and set plot limits
                text(Subtitle_X,Subtitle_Y,Panel_Title,'fontsize',fs,'Interpreter','latex',...
                    'Parent',axh{ii},'Units','Normalized','VerticalAlignment', 'bottom',...
                    'HorizontalAlignment', 'center');

                xlim(axh{ii},X_limits);
                ylim(axh{ii},Y_limits);
                xticks(axh{ii},linspace(X_limits(1),X_limits(end),6))
                yticks(axh{ii},linspace(Y_limits(1),0,6))
                if ii ~= 3
                    xticklabels(axh{ii},'')
                end
            end
        end
    end


    %% Loop through structures to get equilibrium bond lengths
    for ii = 1:N_Struc

        Structure = Structures{ii};
        Color = Colors(ii,:); % For Circles

        % Find equilibrium lattice parameters of given Salt/Structure/Theory
        if contains(Theory,{'TF' 'JC'}) % empirical models
            load(fullfile(datadir,'GROMACS',[Salt '_' Structure '_Lattice_Energies.mat']),'Data');
            
            ModTag = Theory;
            if S_D ~= 1
                ModTag = [ModTag '_D' regexprep(num2str(S_D,'%8.5f'),{'\.' '-'},{'P' 'N'})];
            end
            if S_R ~= 1
            	ModTag = [ModTag '_R' regexprep(num2str(S_R,'%8.5f'),{'\.' '-'},{'P' 'N'})];
            end
            if S_E ~= 1
                ModTag = [ModTag '_E' regexprep(num2str(S_E,'%8.5f'),{'\.' '-'},{'P' 'N'})];
            end
            if S_S ~= 1
                ModTag = [ModTag '_S' regexprep(num2str(S_S,'%8.5f'),{'\.' '-'},{'P' 'N'})];
            end
            if S_MMD ~= 1
                ModTag = [ModTag '_MMD' regexprep(num2str(S_MMD,'%8.5f'),{'\.' '-'},{'P' 'N'})];
            end
            if S_XXD ~= 1
                ModTag = [ModTag '_XXD' regexprep(num2str(S_XXD,'%8.5f'),{'\.' '-'},{'P' 'N'})];
            end
            if S_MXD ~= 1
                ModTag = [ModTag '_MXD' regexprep(num2str(S_MXD,'%8.5f'),{'\.' '-'},{'P' 'N'})];
            end

            try
                SubData = Data.(Salt).(Structure).(ModTag);
            catch
                warning(['Missing Data for: ' Salt ' ' Structure ' ' ModTag]);
                continue
            end
            DT = [SubData{:,9}];
            DatInd = ismember(DT,Data_Type);

            if sum(DatInd) == 0
                continue
            end
            
            CryStruc.a = [SubData{DatInd,1}];
            CryStruc.b = [SubData{DatInd,2}];
            CryStruc.c = [SubData{DatInd,3}];
            CellAngles = SubData(DatInd,8);
            CryStruc.FC = SubData(DatInd,6);   
            if sum(DatInd) > 1
                En = [SubData{DatInd,7}];
                [~,midx] = min(En);
                CryStruc.a = CryStruc.a(midx);
                CryStruc.b = CryStruc.b(midx);
                CryStruc.c = CryStruc.c(midx);
                CryStruc.alpha = CellAngles{midx}(1);
                CryStruc.beta = CellAngles{midx}(2);
                CryStruc.gamma = CellAngles{midx}(3);
                CryStruc.FC = CryStruc.FC{midx};
            else
                CryStruc.alpha = CellAngles{1}(1);
                CryStruc.beta = CellAngles{1}(2);
                CryStruc.gamma = CellAngles{1}(3);
                CryStruc.FC = CryStruc.FC{1};
            end
        else % DFT data
            load(fullfile(datadir,'CRYSTAL',[Salt '_' Structure '_Lattice_Energies.mat']),'Data');
            
            Label = [Salt '_' Structure(1) '_' XC_Func '_D3'];
            SubData = Data.(Label);

            DOI = [];
            for i = 1:size(SubData,1)
                if strcmp(SubData{i,1},Basis_Set)
                    DOI = SubData{i,2};
                    break
                end
            end
            if ~isempty(DOI)
                DT = DOI{10};

                DatInd = ismember(DT,Data_Type);

                if sum(DatInd) == 0
                    continue
                end

                CryStruc.a = DOI{2}(DatInd);
                CryStruc.b = DOI{3}(DatInd);
                CryStruc.c = DOI{4}(DatInd);
                CryStruc.alpha = DOI{11}(DatInd);
                CryStruc.beta = DOI{12}(DatInd);
                CryStruc.gamma = DOI{13}(DatInd);
                CryStruc.FC = DOI{7}(DatInd);        
                if sum(DatInd) > 1
                    En = DOI{8}(DatInd);
                    [~,midx] = min(En);
                    CryStruc.a = CryStruc.a(midx);
                    CryStruc.b = CryStruc.b(midx);
                    CryStruc.c = CryStruc.c(midx);
                    CryStruc.alpha = CryStruc.alpha(midx);
                    CryStruc.beta = CryStruc.beta(midx);
                    CryStruc.gamma = CryStruc.gamma(midx);
                    CryStruc.FC = CryStruc.FC{midx};
                else
                    CryStruc.FC = CryStruc.FC{1};
                end
            else
                continue
            end
        end
        
        % Generate coordination shells using structural parameters
        [MX_BL,MM_BL,XX_BL,XM_BL] = Gen_Stars(CryStruc,Structure,Metal,Halide,Rij(end));
        
        % Collapse coordination shells if less than the specified merge tolerance
        for ix = size(MX_BL,1):-1:2
            if abs(MX_BL(ix,1) - MX_BL(ix-1,1)) < Merge_Tol
                MX_BL(ix-1,1) = (MX_BL(ix,1)*MX_BL(ix,2) + MX_BL(ix-1,1)*MX_BL(ix-1,2))/...
                    (MX_BL(ix,2)+MX_BL(ix-1,2)); % new distance = weighted mean
                MX_BL(ix-1,2) = MX_BL(ix-1,2) + MX_BL(ix,2); % New CN is sum of old
                MX_BL(ix,:) = [];
            end
        end
        for ix = size(MM_BL,1):-1:2
            if abs(MM_BL(ix,1) - MM_BL(ix-1,1)) < Merge_Tol
                MM_BL(ix-1,1) = (MM_BL(ix,1)*MM_BL(ix,2) + MM_BL(ix-1,1)*MM_BL(ix-1,2))/...
                    (MM_BL(ix,2)+MM_BL(ix-1,2)); % new distance = weighted mean
                MM_BL(ix-1,2) = MM_BL(ix-1,2) + MM_BL(ix,2); % New CN is sum of old
                MM_BL(ix,:) = [];
            end
        end
        for ix = size(XX_BL,1):-1:2
            if abs(XX_BL(ix,1) - XX_BL(ix-1,1)) < Merge_Tol
                XX_BL(ix-1,1) = (XX_BL(ix,1)*XX_BL(ix,2) + XX_BL(ix-1,1)*XX_BL(ix-1,2))/...
                    (XX_BL(ix,2)+XX_BL(ix-1,2)); % new distance = weighted mean
                XX_BL(ix-1,2) = XX_BL(ix-1,2) + XX_BL(ix,2); % New CN is sum of old
                XX_BL(ix,:) = [];
            end
        end

        %% Plot Circles   
        for jj = 1:3
            if jj == 1 % M-X
                E = E_Tot_MX;
                X = MX_BL(:,1);
                N_X = MX_BL(:,2);
            elseif jj == 2 % M-M
                E = E_Tot_MM;
                X = MM_BL(:,1);
                N_X = MM_BL(:,2);
            elseif jj == 3
                E = E_Tot_XX; % X-X
                X = XX_BL(:,1); % X-X
                N_X = XX_BL(:,2);
            end

            % Interpolate values
            eqE = interp1(Rij,E,X,'spline');

            if jj == 1
                eqE_MX = interp1(Rij,E_Tot_MX,MX_BL(:,1),'spline');
                eqE_MM = interp1(Rij,E_Tot_MM,MM_BL(:,1),'spline');
                eqE_XX = interp1(Rij,E_Tot_XX,XX_BL(:,1),'spline');
                Total_DE = sum(MX_BL(:,2).*eqE_MX)+sum(MM_BL(:,2).*eqE_MM)+sum(XX_BL(:,2).*eqE_XX);
                MX_DE = sum(MX_BL(:,2).*eqE_MX);
                MM_DE = sum(MM_BL(:,2).*eqE_MM);
                XX_DE = sum(XX_BL(:,2).*eqE_XX);
                disp([Structure ' Total ' Theory ' Dispersion Energy = ' num2str(Total_DE,'%3.2f') ' kJ/mol']);
                disp([Structure ' ' Metal '-' Halide ' ' Theory ...
                    ' Dispersion Energy = ' num2str(MX_DE,'%3.2f') ' kj/mol (' num2str(100*MX_DE/Total_DE,'%3.0f') '%)'])
                disp([Structure ' ' Metal '-' Metal ' ' Theory ...
                    ' Dispersion Energy = ' num2str(MM_DE,'%3.2f') ' kj/mol (' num2str(100*MM_DE/Total_DE,'%3.0f') '%)'])
                disp([Structure ' ' Halide '-' Halide ' ' Theory ...
                    ' Dispersion Energy = ' num2str(XX_DE,'%3.2f') ' kj/mol (' num2str(100*XX_DE/Total_DE,'%3.0f') '%)'])
            end
            
            % Plot circle at equilibrium bond lengths
            if Plotswitch
                if strcmp(Star_Area,'CN')
                    ScA = N_X;
                    prec = '%u';
                    units = '';
                elseif strcmp(Star_Area,'EP')                
                    ScA = 100.*abs(N_X.*eqE)./abs(Total_DE);
                    prec = '%2.1f';
                    units = '\%';
                elseif strcmp(Star_Area,'E')
                    ScA = abs(N_X.*eqE);
                    prec = '%2.1f';
                    units = ' kJ mol$^{-1}$';
                end
                for kk = 1:length(eqE)
                    if kk == 1 && jj == 1 && tidx == 1
                        n = n+1;
                        if Plot3D
                            pobj(n) = scatter3(axh{jj},X(kk),eqE(kk),ScA(kk),AC*ScA(kk),Color,...
                                'filled','MarkerEdgeColor','k','LineWidth',lw,'MarkerFaceAlpha',Transparency);
                        else
                            pobj(n) = scatter(axh{jj},X(kk),eqE(kk),AC*ScA(kk),Color,...
                                'filled','MarkerEdgeColor','k','LineWidth',lw,'MarkerFaceAlpha',Transparency);
                        end
                    else
                        if Plot3D
                            scatter3(axh{jj},X(kk),eqE(kk),ScA(kk),AC*ScA(kk),Color,'filled',...
                                'MarkerEdgeColor','k','LineWidth',lw,'MarkerFaceAlpha',Transparency);
                        else
                            scatter(axh{jj},X(kk),eqE(kk),AC*ScA(kk),Color,'filled',...
                                'MarkerEdgeColor','k','LineWidth',lw,'MarkerFaceAlpha',Transparency);
                        end
                    end
                    if (X(kk) < X_limits(2)) && Text_Overlay
                        text(X(kk),eqE(kk)-1.0,[num2str(ScA(kk),prec) units],'fontsize',fs-10,'Interpreter','latex',...
                        'Parent',axh{jj},'Units','data','VerticalAlignment', 'bottom',...
                        'HorizontalAlignment', 'center')
                    end
                end
            end
        end   
    end
end
if Plotswitch
    % Place legends
    if N_Theor == 1
        [~,legh] = legend(axh{1},pobj,Structures,'FontSize',fs-5,...
            'Interpreter','latex','Box','off','Orientation','vertical',...
            'Position',Legend_Position);
    else
        [~,legh] = legend(axh{1},[pobj; linobj],[Structures Theories],'FontSize',fs-5,...
                'Interpreter','latex','Box','off','Orientation','vertical',...
                'Position',Legend_Position);
    end
    patchobj = findobj(legh,'type','patch'); % Find objects of type 'patch'
    set(patchobj,'MarkerSize', sqrt(AC*6),'Linewidth',lw,'FaceAlpha',Transparency,...
        'BackFaceLighting','lit');



    % Define invisible axis
    AxesH = axes('Units', 'normalized', 'Position', [0,0,1,1], 'visible', 'off', ...
        'YLimMode', 'manual', 'YLim',  [0, 1],'XLim', [-0.357978723404255 1.357978723404255], ...
        'XTick',[],'YTick',[], ...
        'NextPlot','add', ...
        'HitTest', 'off');


    % Y label
    ylabeltext = '$E_{ij}$[D3(BJ)] (kJ mol$^{-1}$)';
    text(0.004,0.5360,ylabeltext,'fontsize',fs,'Interpreter','latex',...
        'rotation',90,'Parent',AxesH,'Units','Normalized','VerticalAlignment','top',...
        'HorizontalAlignment', 'center');

    % X label
    text(0.48,-0.0041,'Bond Length (\AA)','fontsize',fs,'Interpreter','latex',...
        'Parent',AxesH,'Units','Normalized','VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center');

    if strcmp(Star_Area,'CN')
    %     text(0.93,0.62,'Coord. N.','fontsize',fs-5,'Interpreter','latex',...
    %         'Parent',AxesH,'Units','Normalized','VerticalAlignment', 'bottom',...
    %         'HorizontalAlignment', 'center');

        Const = AC*24/(0.025^2); % empirically derived constant
        Ni = [1 2 4 6 8 12 24];
        Dy = 0.052;
        for idx = 1:length(Ni)
            r = sqrt(AC*Ni(idx)/Const); % circle radius

            viscircles(AxesH,[1.16 0.65-Dy^0.9*idx],r,'Color','k','LineWidth',lw);
            text(0.95,0.615-Dy^0.9*idx,num2str(Ni(idx)),'fontsize',fs,'Interpreter','latex',...
            'Parent',AxesH,'Units','Normalized','VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'center')
        end
    elseif strcmp(Star_Area,'EP')
    %     text(0.93,0.59,'E Contr.','fontsize',fs-5,'Interpreter','latex',...
    %         'Parent',AxesH,'Units','Normalized','VerticalAlignment', 'bottom',...
    %         'HorizontalAlignment', 'center');

        Const = AC*43.4067/(0.015^2); % empirically derived constant
        Ni = [1 2 5 10 25 50 100];
        Dy = 0.052;
        for idx = 1:length(Ni)
            r = sqrt(AC*Ni(idx)/Const); % circle radius

            viscircles(AxesH,[1.16 0.62-Dy^0.9*idx],r,'Color','k','LineWidth',lw);
            text(0.94,0.585-Dy^0.9*idx,[num2str(Ni(idx)) '\%'],'fontsize',fs,'Interpreter','latex',...
            'Parent',AxesH,'Units','Normalized','VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'center')
        end

    elseif strcmp(Star_Area,'E')
    %     text(0.93,0.59,'E Contr.','fontsize',fs-5,'Interpreter','latex',...
    %         'Parent',AxesH,'Units','Normalized','VerticalAlignment', 'bottom',...
    %         'HorizontalAlignment', 'center');

        Const = AC*43.4067/(0.032^2); % empirically derived constant
        Ni = [1 2 4 8 10 15 25];
        Dy = 0.052;
        for idx = 1:length(Ni)
            r = sqrt(AC*Ni(idx)/Const); % circle radius

            viscircles(AxesH,[1.14 0.63-Dy^0.9*idx],r,'Color','k','LineWidth',lw);
            text(0.945,0.608-Dy^0.9*idx,[num2str(Ni(idx)) ' kJ mol$^{-1}$'],'fontsize',fs-10,'Interpreter','latex',...
            'Parent',AxesH,'Units','Normalized','VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'center')
        end
    end


    %viscircles(AxesH,[0.285 0.854],0.015,'Color','y','LineWidth',lw,'Linestyle',':')
    %ScA
end
