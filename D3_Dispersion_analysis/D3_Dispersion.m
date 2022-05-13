%% Inputs
Salt = 'LiF';
Theory = 'PW1PW';
Structures = {'Rocksalt'}; %{'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'};
if ispc
    datadir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\CRYSTAL';
else
    datadir = '/home/user/ThesisCodeBase/data/CRYSTAL';
end
Data_Type = 1;
Basis_Set = 'pob-TZVP';
Scale_Factor = 0.3:0.05:3;
Interp_Step = 0.001;
Show_CN = true; % true to plot coordination number over
Fix_CN = true; % When true, forces the fractional coordination number to stay >= 1
X_Axis = 'SF'; % a, SF, or BL
fs = 35; % Font size
lw = 1.5; % Line width
Y_limits = [-75 25];
BL_Limit = [1 6];
Subtitle_X = 0.10;
Subtitle_Y = 0.80;

Panel_Position{1} = [0.10375 0.69 0.7109 0.29]; %[X Y Width Height]
Panel_Position{2} = [0.10375 0.40 0.7109 0.29];
Panel_Position{3} = [0.10375 0.11 0.7109 0.29];
Legend_Position = [0.8148 0.6390 0.1777 0.3078];
N_Struc = length(Structures);
Colors = cbrewer('qual','Set1',N_Struc,'PCHIP');
Panel_Title = cell(3,1);
LineStyles = {'-' '-.' '--' '--' '-' '--' '-'}; %'-' (default) | '--' | ':' | '-.' | 'none')
[Metal,Halide] = Separate_Metal_Halide(Salt);
%% Global parameters of D3(BJ) and XC specific parameters
Version = 4; % D3(BJ)
TZ = false;
[s6,a1,s8,a2,~] = setfuncpar(Theory,Version,TZ);
k1 = 16;
k2 = 4/3;
k3 = -4;

%% Parallel stuff
PrefCores = feature('numcores');

if ~isempty(gcp('nocreate'))
    Cur_Pool = gcp;
    Cur_Workers = Cur_Pool.NumWorkers;
    if Cur_Workers < PrefCores
        delete(Cur_Pool);
        ppool = parpool(PrefCores);
    else
        ppool = Cur_Pool;
    end
else
    ppool = parpool(PrefCores); 
end

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

%% Load C6 function,max coordination number function (mxc) vs atomic number, and r2r4 from disc
[C6_Func,mxc] = copyc6;
loadr2r4 = load('r2r4.mat','r2r4');
r2r4 = loadr2r4.r2r4;

%% Load figure
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
axh = cell(3,1);
MA = length(Scale_Factor);
set(figh,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);

%% Preallocate arrays
E6_MX = nan(MA,N_Struc);
E8_MX = nan(MA,N_Struc);
E_Tot_MX = nan(MA,N_Struc);

E6_MM = nan(MA,N_Struc);
E8_MM = nan(MA,N_Struc);
E_Tot_MM = nan(MA,N_Struc);
M_CN = nan(MA,N_Struc);

E6_XX = nan(MA,N_Struc);
E8_XX = nan(MA,N_Struc);
E_Tot_XX = nan(MA,N_Struc);
X_CN = nan(MA,N_Struc);

Plot_X = nan(MA,N_Struc);
Plot_X2 = nan(MA,N_Struc); % for 2nd neighbour bond lengths if applicable
pobj = gobjects(N_Struc,1);
Ref_X = nan(N_Struc,1);
Ref_X2 = nan(N_Struc,1); % for 2nd neighbour bond lengths

%% Loop through structures
for jj = 1:N_Struc
       
    % Find equilibrium lattice parameters of given Salt/Structure/XC functional
    Structure = Structures{jj};
    load(fullfile(datadir,[Salt '_' Structure '_Lattice_Energies.mat']),'Data');
    
    Label = [Salt '_' Structure(1) '_' Theory '_D3'];
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
        
        CellIn.a = DOI{2}(DatInd);
        CellIn.b = DOI{3}(DatInd);
        CellIn.c = DOI{4}(DatInd);
        CellIn.BL = DOI{5}(DatInd);
        CellIn.BL2 = DOI{9}(DatInd);
        CellIn.alpha = DOI{11}(DatInd);
        CellIn.beta = DOI{12}(DatInd);
        CellIn.gamma = DOI{13}(DatInd);
        eq_FC = DOI{7}(DatInd);
        
        if length(eq_FC) > 1
            En = DOI{8}(DatInd);
            [~,mind] = min(En);
            CellIn.a = CellIn.a(mind);
            CellIn.b = CellIn.b(mind);
            CellIn.c = CellIn.c(mind);
            CellIn.BL = CellIn.BL(mind);
            CellIn.BL2 = CellIn.BL2(mind);
            CellIn.alpha = CellIn.alpha(mind);
            CellIn.beta = CellIn.beta(mind);
            CellIn.gamma = CellIn.gamma(mind);
            eq_FC = eq_FC{mind};
        else
            eq_FC = eq_FC{1};
        end
    else
        continue
    end
    
    Met_Ind = strcmpi(eq_FC(:,1),Metal);
    Hal_Ind = strcmpi(eq_FC(:,1),Halide);
    CellIn.FC_Metal = mod([eq_FC{Met_Ind,2:4}],1);
    CellIn.FC_Halide = mod([eq_FC{Hal_Ind,2:4}],1);
    
    % number of atoms in unit cell
    switch lower(Structure)
        case 'rocksalt'
            CellIn.N = 8; 
        case 'wurtzite'
            CellIn.N = 4;
        case 'betabeo'
            CellIn.N = 8;
        case 'cscl'
            CellIn.N = 2;
        case 'fivefive'
            CellIn.N = 8;
        case 'nias'
            CellIn.N = 4;
        case 'sphalerite'
            CellIn.N = 8;
        otherwise
            error(['Unknown structure: ' Structure])
    end
    
    [CellIn.FC_Metal,CellIn.FC_Halide] = ...
        UnitCell_FractionalCoords(CellIn.FC_Metal,...
        CellIn.FC_Halide,Structure);

    Frac_Coordinates = [CellIn.FC_Metal; CellIn.FC_Halide];

    Identities_UC = [repmat(string(Metal),CellIn.N/2,1);
        repmat(string(Halide),CellIn.N/2,1)];
    N = CellIn.N; % Atoms per unit cell

    switch X_Axis
        case 'a'
            Ref_X(jj) = CellIn.a;
            Ref_X2(jj) = CellIn.a;
        case 'BL'
            Ref_X(jj) = CellIn.BL;
            Ref_X2(jj) = CellIn.BL2;
        case 'SF'
            Ref_X(jj) = 1.0;
            Ref_X2(jj) = 1.0;
        otherwise
            error(['Unknown X_Axis Type: ' X_Axis '.']);
    end
    
    parfor ii = 1:length(Scale_Factor) % parfor
        %% Scale cell
        CellParams = CellIn;
        CellParams = Scale_Structure(CellParams,Scale_Factor(ii));
        
        switch X_Axis
            case 'a'
                Plot_X(ii,jj) = CellParams.a;
                Plot_X2(ii,jj) = CellParams.a;
            case 'BL'
                Plot_X(ii,jj) = CellParams.BL;
                Plot_X2(ii,jj) = CellParams.BL2;
            case 'SF'
                Plot_X(ii,jj) = Scale_Factor(ii);
                Plot_X2(ii,jj) = Scale_Factor(ii);
            otherwise
                error(['Unknown X_Axis Type: ' X_Axis '.']);
        end
        
        %% Atomic numbers of atoms in unit cell
        AN_UC = nan(N,1);
        for i = 1:N
            AN_UC(i) = elements('Symbol',Identities_UC(i),'atomic_number'); %#ok<PFBNS>
        end

        %% Generate transformation matrix
        abc = [CellParams.a 0 0; 0 CellParams.b 0; 0 0 CellParams.c];
        Transform_Matrix = abc*GenTransformMatrix(CellParams);
        a_vec = Transform_Matrix(1,:);
        b_vec = Transform_Matrix(2,:);
        c_vec = Transform_Matrix(3,:);
        
        %% Calculate size of search space in each direction based on dispersion cutoff and shortest distance between planes
        Ma = ceil(max(Disp_Cutoff/(CellParams.a*cosd(CellParams.gamma-90)),...
            Disp_Cutoff/(CellParams.a*cosd(CellParams.beta-90))));
        Mb = ceil(max(Disp_Cutoff/(CellParams.b*cosd(CellParams.gamma-90)),...
            Disp_Cutoff/(CellParams.b*cosd(CellParams.alpha-90))));
        Mc = ceil(max(Disp_Cutoff/(CellParams.c*cosd(CellParams.alpha-90)),...
            Disp_Cutoff/(CellParams.c*cosd(CellParams.beta-90))));

        %% Expand into 2Ma+1 x 2Mb+1 x 2Mc+1 Supercell of cartesian coordinates       
        Search_N = Frac_Coordinates*Transform_Matrix;
        Search_a = [0 -Ma:-1 1:Ma]' * a_vec;
        Search_b = [0 -Mb:-1 1:Mb]' * b_vec;
        Search_c = [0 -Mc:-1 1:Mc]' * c_vec;
        
        Search_Combs = combvec(1:N,1:(2*Ma+1),1:(2*Mb+1),1:(2*Mc+1));
        Identities = Search_Combs(1,:);
        SCa = Search_Combs(2,:);
        SCb = Search_Combs(3,:);
        SCc = Search_Combs(4,:);
        Cartesian_Coordinates = Search_N(Identities,:) + Search_a(SCa,:) + Search_b(SCb,:) + Search_c(SCc,:);
        
        %% Calculating fractional coordination number and C6 coefficients
        Frac_CN_Cell = nan(N,1); % Fractional coordination number for each atom in unit cell
        Disp_Search_Identities = cell(N,1); % Identities of all atoms within dispersion search radius for each atom in unit cell
        Disp_Search_rAB = cell(N,1);
        for idx = 1:N % For each atom in unit cell
            %% Initialize reference
            % Coordinates of reference atom
            Reference_Atom = Cartesian_Coordinates(idx,:);
            % Atomic species of reference atom
            Reference_ID = Identities_UC(idx);
            % Covalent radius of reference atom
            RAcov = Covalent_Radii_D3(Reference_ID);

            %% Pick out atoms closer than Dispersion Cutoff from reference unit atom
            rAB = vecnorm(Reference_Atom - Cartesian_Coordinates,2,2);
            
            % Which are closer than cutoff?
            disp_idx = (rAB <= Disp_Cutoff);
            % Skip reference atom
            disp_idx(idx) = false;
            
            Disp_Search_Coordinates = Cartesian_Coordinates(disp_idx,:);
            Disp_Identities = Identities(disp_idx);
            Disp_rAB = rAB(disp_idx);

            %% Save into cell
            Disp_Search_Identities{idx} = Disp_Identities;
            Disp_Search_rAB{idx} = Disp_rAB;

            %% Pick out atoms closer than CN-Cutoff from reference atom
            CN_idx = (Disp_rAB <= CN_Cutoff);
            CN_Search_Coordinates = Disp_Search_Coordinates(CN_idx,:);
            CN_Identities = Disp_Identities(CN_idx);
            CN_rAB = Disp_rAB(CN_idx);

            %% Calculate fractional Coordination Number of reference atom
            PN = length(CN_rAB);
            CN_Part = nan(PN,1);
            
            RBcov_M = Covalent_Radii_D3(Metal);
            RBcov_X = Covalent_Radii_D3(Halide);
            M_ind = Identities_UC(CN_Identities) == Metal;
            X_ind = ~M_ind;
            
            % Split fractional coordination computation into 2 parts for
            % vectorization speed up
            CN_Part_M = 1 ./ ( 1 + exp( -k1.*(k2*(RAcov + RBcov_M)./CN_rAB(M_ind) - 1) ) );
            CN_Part_X = 1 ./ ( 1 + exp( -k1.*(k2*(RAcov + RBcov_X)./CN_rAB(X_ind) - 1) ) );
            % Fractional coordination numbers of reference atom
            Frac_CN_Cell(idx) = sum(CN_Part_M)+sum(CN_Part_X);
        end
        if Fix_CN
            Frac_CN_Cell = max(Frac_CN_Cell,1);
        end
        
        M_CN(ii,jj) = Frac_CN_Cell(1);
        X_CN(ii,jj) = Frac_CN_Cell(end);


        %% Calculate C6 and C8 for each possible pair of atoms in unit cell, as well as R0AB
        C6AB = nan(N,N);
        C8AB = nan(N,N);
        R0AB = nan(N,N);
        for i = 1:N
            Z_A = AN_UC(i); % Atomic number of atom A
            CN_A = Frac_CN_Cell(i); % Fractional coordination number of atom A
            sqrt_Q_A = r2r4(Z_A); % Factor used to calculate C8
            for j = 1:N
                Z_B = AN_UC(j); % Atomic number of atom B
                CN_B = Frac_CN_Cell(j); % Fractional coordination number of atom B
                sqrt_Q_B = r2r4(Z_B); % Factor used to calculate C8

                % Calculate C6 and C6 matrix
                C6AB(i,j) = getc6(C6_Func,mxc,Z_A,Z_B,CN_A,CN_B)*c6units; % in kJ/mol Ang^6
                C8AB(i,j) = 3.0*getc6(C6_Func,mxc,Z_A,Z_B,CN_A,CN_B)*sqrt_Q_A*sqrt_Q_B*c8units; % in kJ/mol Ang^8
                R0AB(i,j) = sqrt(C8AB(i,j)/C6AB(i,j)); % in Angstroms               
            end
        end

        %% Calculate total dispersion energy of unit cell
        E6_mx = nan(N,1);
        E8_mx = nan(N,1);
        E6_mm = nan(N,1);
        E8_mm = nan(N,1);
        E6_xx = nan(N,1);
        E8_xx = nan(N,1);
        mx = 1;
        mm = 1;
        xx = 1;
        Z_M = AN_UC(1);
        Z_X = AN_UC(N);
        for idx = 1:N

            Z_A = AN_UC(idx);
            QM = length(Disp_Search_rAB{idx});
            E6_Part = zeros(QM,3); % column 1 == MX, column 2 == MM, column 3 == XX
            E8_Part = zeros(QM,3);
            for i = 1:QM
                Z_B = AN_UC(Disp_Search_Identities{idx}(i));

                C6 = C6AB(idx,Disp_Search_Identities{idx}(i));
                C8 = C8AB(idx,Disp_Search_Identities{idx}(i));
                R_AB = Disp_Search_rAB{idx}(i); % Distance between atoms (Angstrom)
                R0 = R0AB(idx,Disp_Search_Identities{idx}(i)); % vDW radius in angstroms
                F_R0 = a1*R0 + a2*Bohr_Ang; % XC dependent function (Angstroms);

                if (Z_A ~= Z_B) % MX
                    E6_Part(i,1) = -s6*C6/((R_AB^6) + (F_R0^6));
                    E8_Part(i,1) = -s8*C8/((R_AB^8) + (F_R0^8));
                elseif Z_A == Z_M % MM
                    E6_Part(i,2) = -s6*C6/((R_AB^6) + (F_R0^6));
                    E8_Part(i,2) = -s8*C8/((R_AB^8) + (F_R0^8));
                elseif Z_A == Z_X %% XX
                    E6_Part(i,3) = -s6*C6/((R_AB^6) + (F_R0^6));
                    E8_Part(i,3) = -s8*C8/((R_AB^8) + (F_R0^8));
                end
            end

            E6_mx(idx) = sum(E6_Part(:,1)); % C6 energy
            E8_mx(idx) = sum(E8_Part(:,1)); % C8 energy
            E6_mm(idx) = sum(E6_Part(:,2)); % C6 energy
            E8_mm(idx) = sum(E8_Part(:,2)); % C8 energy
            E6_xx(idx) = sum(E6_Part(:,3)); % C6 energy
            E8_xx(idx) = sum(E8_Part(:,3)); % C8 energy

        end

        E6_MX(ii,jj) = (sum(E6_mx)/2)/(N/2);
        E8_MX(ii,jj) = (sum(E8_mx)/2)/(N/2);
        E_Tot_MX(ii,jj) = E6_MX(ii,jj) + E8_MX(ii,jj); % in kJ/mol of ion pairs

        E6_MM(ii,jj) = (sum(E6_mm)/2)/(N/2);
        E8_MM(ii,jj) = (sum(E8_mm)/2)/(N/2);
        E_Tot_MM(ii,jj) = E6_MM(ii,jj) + E8_MM(ii,jj); % in kJ/mol of ion pairs

        E6_XX(ii,jj) = (sum(E6_xx)/2)/(N/2);
        E8_XX(ii,jj) = (sum(E8_xx)/2)/(N/2);
        E_Tot_XX(ii,jj) = E6_XX(ii,jj) + E8_XX(ii,jj); % in kJ/mol of ion pairs
    end

    %% Plot stuff
    Color = Colors(jj,:);
    for j = 1:3 % Loop over the three plots
        if j == 1 % M-X
            X = Plot_X(:,jj);
            RX = Ref_X(jj);
            idx1 = 1;
            idx2 = N;
            E = E_Tot_MX(:,jj);
        elseif j == 2 % M-M
            X = Plot_X2(:,jj);
            RX = Ref_X2(jj);
            idx1 = 1;
            idx2 = 1;
            E = E_Tot_MM(:,jj);
        elseif j == 3
            X = Plot_X2(:,jj);
            RX = Ref_X2(jj);
            idx1 = N;
            idx2 = N;
            E = E_Tot_XX(:,jj);
        end

        if jj == 1
            % Generate axis
            axh{j} = axes('Parent',figh,'Position',Panel_Position{j},'FontSize',fs,...
                'box','on','TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on');
            hold(axh{j},'on')
            
            % Plot zero line
            yline(axh{j},0,'Color','k','LineStyle','--','LineWidth',lw);
        end

        if (j == 2 || j == 3) && Show_CN
            yyaxis(axh{j},'left')
        end
        
        X_v = X(1):Interp_Step:X(end);
        E_v = interp1(X,E',X_v,'spline');
        % plot energy curve
        p = plot(axh{j},X_v,E_v,'Color',Color,'LineStyle',LineStyles{jj},'LineWidth',lw+2,'Marker','none');
        
        % Plot circle at equilibrium bond length
        eqE = interp1(X_v,E_v,RX,'spline');
        scatter(axh{j},RX,eqE,50,Color,'filled','MarkerEdgeColor','k','LineWidth',lw);
        
        if j == 1
            pobj(jj) = p;
        end

        if j == 1 % M-X
            if jj == 1
                Panel_Title{j} = [Metal '$^{+}$' Halide '$^{-}$'];
            end
        elseif j == 2 % M-M
            if Show_CN
                yyaxis(axh{j},'right')
                M_CN_v = interp1(X,M_CN(:,jj),X_v,'makima');
                plot(axh{j},X_v,M_CN_v,'Color',Color,'LineStyle',':','LineWidth',lw,'Marker','none');
            end
            if jj == 1
                if Show_CN
                    ylim(axh{j},[0 2]);
                    yticks(axh{j},0:0.5:2.0);
                    yticklabels(axh{j},{'0' '' '1' '' ''})
                    % Plot one line
                    yline(axh{j},1,'Color','k','LineStyle',':','LineWidth',lw);
                end
                Panel_Title{j} = [Metal '$^{+}$' Metal '$^{+}$'];
            end
            if Show_CN
                yyaxis(axh{j},'left')
            end
        elseif j == 3
            if Show_CN
                yyaxis(axh{j},'right')
                X_CN_v = interp1(X,X_CN(:,jj),X_v,'makima');
                plot(axh{j},X_v,X_CN_v,'Color',Color,'LineStyle',':','LineWidth',lw,'Marker','none');
            end
            if jj == 1
                if Show_CN
                    ylim(axh{j},[0 2]);
                    yticks(axh{j},0:0.5:2.0);
                    yticklabels(axh{j},{'0' '' '1' '' ''})
                    % Plot one line
                    yline(axh{j},1,'Color','k','LineStyle',':','LineWidth',lw);
                end
                Panel_Title{j} = [Halide '$^{-}$' Halide '$^{-}$'];
            end
            if Show_CN
                yyaxis(axh{j},'left')
            end
        end
    end
end

X_limits = [floor(min([Plot_X Plot_X2],[],'all')) ceil(max([Plot_X Plot_X2],[],'all'))]; % Limits of x axis

switch X_Axis
    case 'a'
        xlabeltext = 'a (\AA)';
    case 'BL'
        xlabeltext = 'Bond Length (\AA)';
        X_limits = BL_Limit;
    case 'SF'
        xlabeltext = 'Scale Factor';
    otherwise
        error(['Unknown X_Axis Type: ' X_Axis '.']);
end


for j = 1:3
    % Put in axis sub-title and set plot limits
    text(Subtitle_X,Subtitle_Y,Panel_Title{j},'fontsize',fs,'Interpreter','latex',...
        'Parent',axh{j},'Units','Normalized','VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center');

    xlim(axh{j},X_limits);
    ylim(axh{j},Y_limits);
    xticks(axh{j},linspace(X_limits(1),X_limits(end),6))
    if j ~= 3
        xticklabels(axh{j},'')
    end
end

LegText = strrep(Structures,'Sphalerite','\begin{tabular}{p{1.02cm} p{0.1cm} p{1.2cm}} Sphal. & 4 & 12\end{tabular}');
LegText = strrep(LegText,'FiveFive','\begin{tabular}{p{1.02cm}p{0.1cm} p{1.2cm}} 5-5 & 5 & 6/2\end{tabular}');
LegText = strrep(LegText,'BetaBeO','\begin{tabular}{p{1.02cm} p{0.1cm} p{1.2cm}} $\beta$-BeO & 4 & 1/2/8\end{tabular}');
LegText = strrep(LegText,'Rocksalt','\begin{tabular}{p{1.02cm} p{0.1cm} p{1.2cm}} Rock. & 6 & 12\end{tabular}');
LegText = strrep(LegText,'NiAs','\begin{tabular}{p{1.02cm} p{0.1cm} p{1.2cm}} NiAs & 6 & (2)/6/(6)\end{tabular}');
LegText = strrep(LegText,'Wurtzite','\begin{tabular}{p{1.02cm} p{0.1cm} p{1.2cm}} Wurtz. & 4 & 12\end{tabular}');
LegText = strrep(LegText,'CsCl','\begin{tabular}{p{1.02cm} p{0.1cm} p{1.2cm}} CsCl & 8 & 6\end{tabular}');

% Place legend
legend(axh{1},pobj,LegText,'FontSize',fs-10,...
    'Interpreter','latex','Box','off','Orientation','vertical',...
    'Position',Legend_Position);


% Define invisible axis
AxesH = axes('Units', 'normalized', 'Position', [0,0,1,1], 'visible', 'off', ...
    'YLimMode', 'manual', 'YLim',  [0, 1], ...
    'XTick',[],'YTick',[], ...
    'NextPlot','add', ...
    'HitTest', 'off');

% Y label
ylabeltext = '$E$[D3(BJ)] (kJ mol$^{-1}$)';
text(0.004,0.5360,ylabeltext,'fontsize',fs,'Interpreter','latex',...
    'rotation',90,'Parent',AxesH,'Units','Normalized','VerticalAlignment','top',...
    'HorizontalAlignment', 'center');

% Y label 2
if Show_CN
    ylabeltext = [Metal '$^{+}$ FCN'];
    text(0.8528,0.5360,ylabeltext,'fontsize',fs,'Interpreter','latex',...
        'rotation',90,'Parent',AxesH,'Units','Normalized','VerticalAlignment','top',...
        'HorizontalAlignment', 'center');

    % Y label 3
    ylabeltext = [Halide '$^{-}$ FCN'];
    text(0.8528,0.2551,ylabeltext,'fontsize',fs,'Interpreter','latex',...
        'rotation',90,'Parent',AxesH,'Units','Normalized','VerticalAlignment','top',...
        'HorizontalAlignment', 'center');
end

% X label
text(0.48,-0.0041,xlabeltext,'fontsize',fs,'Interpreter','latex',...
    'Parent',AxesH,'Units','Normalized','VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center');

% Coordination number text
coordination = '\begin{tabular}{p{1.02cm} p{0.1cm} p{0.3cm}} \hfill CN: & 1$^{\textrm{st}}$ & 2$^{\textrm{nd}}$\\ \hline \end{tabular}';
text(0.9176,0.9455,coordination,'fontsize',fs-10,'Interpreter','latex',...
    'Parent',AxesH,'Units','Normalized','VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center');


