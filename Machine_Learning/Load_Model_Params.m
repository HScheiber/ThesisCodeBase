function [Settings,ModelFound] = Load_Model_Params(Settings)

ModelFound = true;
if strcmp(Settings.Theory,'JC') && strcmp(Settings.Salt,'LiI') && strcmp(Settings.Model,'BF1b')
    damp_MM = true;
    Settings.Model = 'BF1';
else
    damp_MM = false;
end

[~,~,~,~,Model_Data_Loc] = find_home;

Model_Filename = fullfile(Model_Data_Loc,[Settings.Salt '_' Settings.Theory '_Model_' Settings.Model '_fullopt.mat']);
BO_Filename = fullfile(Model_Data_Loc,[Settings.Salt '_' Settings.Theory '_Model_' Settings.Model '_bayesopt.mat']);

Settings.S = Init_Scaling_Object;
Settings.C6_Damp = Init_C6Damping_Object;
Settings.CR_Damp = Init_CRDamping_Object;
[Settings.GAdjust_MX,Settings.GAdjust_MM,Settings.GAdjust_XX] = Init_GAdjust_Object;
    
% Special case
if isempty(Settings.Model)
    return
elseif strcmp(Settings.Theory,'JC') && strcmp(Settings.Salt,'LiI') && strcmp(Settings.Model,'A')
    disp(['Model found. Loading ' Settings.Theory ' Model ' Settings.Model])
    Settings.CR_Damp.MM.r_d = 0.26; % liCl = 0.21, LiI = 0.26
    Settings.CR_Damp.MM.b  = 75;
    Settings.S.D.MM = 910;
    Settings.S.D.XX = 0.23;
    return
end

if isfile(Model_Filename) && isfile(BO_Filename)
    pdat = load(Model_Filename);
    
    warning('off','MATLAB:class:InvalidSuperClass')
    warning('off','MATLAB:load:classNotFound')
    warning('off','MATLAB:load:classError')
    BO = load(BO_Filename);
    warning('on','MATLAB:class:InvalidSuperClass')
    warning('on','MATLAB:load:classNotFound')
    warning('on','MATLAB:load:classError')
    
    s = functions(BO.results.ObjectiveFcn);
    
    if isfield(s.workspace{1},'Salt')
        BO_WS = s.workspace{1};
    else
        BO_WS = s.workspace{1}.Model;
    end
    
    if isfield(BO_WS,'Fix_Charge')
        Fix_Charge = BO_WS.Fix_Charge;
    else
        Fix_Charge = true;
    end
    if isfield(BO_WS,'Additivity')
        Additivity = BO_WS.Additivity;
    else
        Additivity = true;
    end
    if isfield(BO_WS,'Additional_MM_Disp')
        Additional_MM_Disp = BO_WS.Additional_MM_Disp;
    elseif isfield(BO_WS,'Additional_MM_Force')
        Additional_MM_Disp = BO_WS.Additional_MM_Force;
    else
        Additional_MM_Disp = false;
    end
    if isfield(BO_WS,'Additional_GAdjust')
        Additional_GAdjust = BO_WS.Additional_GAdjust;
    else
        Additional_GAdjust = {};
    end
    if isfield(BO_WS,'SigmaEpsilon')
        SigmaEpsilon = BO_WS.SigmaEpsilon;
    else
        SigmaEpsilon = false;
    end
    if isfield(BO_WS,'C6Damp')
        Settings.C6_Damp = BO_WS.C6Damp;
        
        def_C6Damp = Init_C6Damping_Object;
        if ~isfield(Settings.C6_Damp,'input_rvdw')
            Settings.C6_Damp.input_rvdw = def_C6Damp.input_rvdw;
        end
        if ~isfield(Settings.C6_Damp,'rvdw')
            Settings.C6_Damp.rvdw = def_C6Damp.rvdw;
        end
        if ~isfield(Settings.C6_Damp,'N')
            Settings.C6_Damp.N = def_C6Damp.N;
        end
    elseif isfield(BO_WS,'C6_Damp')
        Settings.C6_Damp = BO_WS.C6_Damp;
        
        def_C6Damp = Init_C6Damping_Object;
        if ~isfield(Settings.C6_Damp,'input_rvdw')
            Settings.C6_Damp.input_rvdw = def_C6Damp.input_rvdw;
        end
        if ~isfield(Settings.C6_Damp,'rvdw')
            Settings.C6_Damp.rvdw = def_C6Damp.rvdw;
        end
        if ~isfield(Settings.C6_Damp,'N')
            Settings.C6_Damp.N = def_C6Damp.N;
        end
    else
        Settings.C6_Damp = Init_C6Damping_Object;
    end
    
    
    if isfield(BO_WS,'CRDamp')
        Settings.CR_Damp = BO_WS.CRDamp;
    elseif isfield(BO_WS,'CR_Damp')
        Settings.CR_Damp = BO_WS.CR_Damp;
    else
        % Used as the default in older models
        Settings.CR_Damp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Settings.CR_Damp.MX.b = 100; % sigmoid "steepness" for damping
        Settings.CR_Damp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Settings.CR_Damp.MM.b  = 75; % 75
        Settings.CR_Damp.XX.r_d = 0.20; 
        Settings.CR_Damp.XX.b  = 100;
    end
    

    if isfield(BO_WS,'Additional_Function')
        Additional_Function = BO_WS.Additional_Function;
    else
        Additional_Function = Init_Additional_Function;
    end
    
    if isfield(BO_WS,'Fix_C8')
        Fix_C8 = BO_WS.Fix_C8;
    else
        Fix_C8 = false;
    end
    
    if isfield(BO_WS,'Fix_Alpha')
        Fix_Alpha = BO_WS.Fix_Alpha;
    else
        Fix_Alpha = false;
    end
    if isfield(BO_WS,'Q_value')
        Q_value = BO_WS.Q_value;
    else
        Q_value = 1;
    end
else
    warning('Model not found.')
    ModelFound = false;
    return
end

if damp_MM
    Settings.C6_Damp.MM = true;
end

if isfield(pdat,'full_opt_point')
    Param = pdat.full_opt_point;
else
    warning('Model not found.')
    ModelFound = false;
    return
end

if strcmp(Settings.Theory,'TF')

    % Loose form of exp-C6-C8 model
    if SigmaEpsilon
        % Default model parameters: all length-scale units are nm,
        % energy scale units are kJ/mol
        DefMDSettings = Initialize_MD_Settings;
        DefMDSettings.Salt = Settings.Salt;
        [defTFMX,defTFMM,defTFXX] = TF_Potential_Parameters(DefMDSettings);

        % Input parameters
        r0_MM = Param(1); % nm
        r0_XX = Param(2); % nm

        if Additivity
            r0_MX = (r0_MM + r0_XX)/2; % nm

            epsilon_MM = Param(3); % kJ/mol
            epsilon_XX = Param(4); % kJ/mol
            epsilon_MX = sqrt(epsilon_MM*epsilon_XX); % kJ/mol

            gamma_MX = Param(5); % Unitless
            gamma_MM = gamma_MX; % Unitless
            gamma_XX = gamma_MX; % Unitless
            pidx = 5;
        else
            r0_MX = Param(3); % nm

            epsilon_MM = Param(4); % kJ/mol
            epsilon_XX = Param(5); % kJ/mol
            epsilon_MX = Param(6); % kJ/mol

            gamma_MM = Param(7); % Unitless
            gamma_XX = Param(8); % Unitless
            gamma_MX = Param(9); % Unitless
            pidx = 9;
        end
        
        if gamma_MX > 48/7
            epsilon_MX = -epsilon_MX;
        end
        if gamma_MM > 48/7
            epsilon_MM = -epsilon_MM;
        end
        if gamma_XX > 48/7
            epsilon_XX = -epsilon_XX;
        end

        % Convert to Condensed form
        alpha_MM = gamma_MM/r0_MM;
        alpha_XX = gamma_XX/r0_XX;
        alpha_MX = gamma_MX/r0_MX;

        B_MM = 48*epsilon_MM*exp(gamma_MM)/(48 - 7*gamma_MM);
        B_XX = 48*epsilon_XX*exp(gamma_XX)/(48 - 7*gamma_XX);
        B_MX = 48*epsilon_MX*exp(gamma_MX)/(48 - 7*gamma_MX);

        C_MM = 4*epsilon_MM*gamma_MM*(r0_MM^6)/(48 - 7*gamma_MM);
        C_XX = 4*epsilon_XX*gamma_XX*(r0_XX^6)/(48 - 7*gamma_XX);
        C_MX = 4*epsilon_MX*gamma_MX*(r0_MX^6)/(48 - 7*gamma_MX);

        D_MM = 3*epsilon_MM*gamma_MM*(r0_MM^8)/(48 - 7*gamma_MM);
        D_XX = 3*epsilon_XX*gamma_XX*(r0_XX^8)/(48 - 7*gamma_XX);
        D_MX = 3*epsilon_MX*gamma_MX*(r0_MX^8)/(48 - 7*gamma_MX);

        % Convert to scaling w.r.t. TF (Pars output)
        % Dispersion scale
        Settings.S.D6D.MM = C_MM/defTFMM.C;
        Settings.S.D6D.XX = C_XX/defTFXX.C;
        Settings.S.D6D.MX = C_MX/defTFMX.C;
        Settings.S.D8D.MM = D_MM/defTFMM.D;
        Settings.S.D8D.XX = D_XX/defTFXX.D;
        Settings.S.D8D.MX = D_MX/defTFMX.D;

        % Repulsive wall exponent scale
        Settings.S.A.MM = alpha_MM/defTFMM.alpha;
        Settings.S.A.XX = alpha_XX/defTFXX.alpha;
        Settings.S.A.MX = alpha_MX/defTFMX.alpha;

        % Repulsive wall prefactor scale
        Settings.S.R.MM = B_MM/defTFMM.B;
        Settings.S.R.XX = B_XX/defTFXX.B;
        Settings.S.R.MX = B_MX/defTFMX.B;

        % Scaling Coulombic Charge
        if Fix_Charge
            Settings.S.Q = Q_value;
        else
            Settings.S.Q = Param(pidx+1);
        end

    % Tight form of exp-C6-C8 model
    else
        % 1/R6 Dispersion (TF only)
        Settings.S.D6D.MM = Param(1);
        Settings.S.D6D.XX = Param(2);
        Settings.S.D6D.MX = Param(3);

        if Fix_C8

            % Calculate value of C8 using recursive relations
            [Metal,Halide] = Separate_Metal_Halide(Settings.Salt);

            % Default TF params: C in units of kJ/mol nm^6, D in units of kJ/mol nm^8
            DefMDSettings = Initialize_MD_Settings;
            DefMDSettings.Salt = Settings.Salt;
            [TF_MX,TF_MM,TF_XX] = TF_Potential_Parameters(DefMDSettings);
            
            % Conversion factors
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

            % Calculate Scaled C8 using recursion relation from D3 paper
            C8.MM = 3.0*(Settings.S.D6D.MM*TF_MM.C/c6units)*sqrt_Q.(Metal)*sqrt_Q.(Metal)*c8units; % in kJ/mol nm^8
            C8.XX = 3.0*(Settings.S.D6D.XX*TF_XX.C/c6units)*sqrt_Q.(Halide)*sqrt_Q.(Halide)*c8units; % in kJ/mol nm^8
            C8.MX = 3.0*(Settings.S.D6D.MX*TF_MX.C/c6units)*sqrt_Q.(Metal)*sqrt_Q.(Halide)*c8units; % in kJ/mol nm^8

            % Update the 1/R8 scaling
            Settings.S.D8D.MM = C8.MM/TF_MM.D;
            Settings.S.D8D.XX = C8.XX/TF_XX.D;
            Settings.S.D8D.MX = C8.MX/TF_MX.D;
            pidx = 3;
        else
            % 1/R8 Dispersion (TF only)
            Settings.S.D8D.MM = Param(4);
            Settings.S.D8D.XX = Param(5);
            Settings.S.D8D.MX = Param(6);
            pidx = 6;
        end
        
        % Alpha (TF exponential steepness repulsive parameter)
        if Fix_Alpha
            Settings.S.A.MM = 1;
            Settings.S.A.XX = 1;
            Settings.S.A.MX = 1;
        else
            Settings.S.A.MM = Param(pidx+1);
            Settings.S.A.XX = Param(pidx+2);
            Settings.S.A.MX = Param(pidx+3);
            pidx = pidx+3;
        end

        % Repulsive wall prefactor
        Settings.S.R.MM = Param(pidx+1);
        Settings.S.R.XX = Param(pidx+2);
        Settings.S.R.MX = Param(pidx+3);
        pidx = pidx+3;

        % Scaling Coulombic Charge
        if Fix_Charge
            Settings.S.Q = Q_value;
        else
            Settings.S.Q = Param(pidx+1);
            pidx = pidx+1;
        end
    end
    
elseif strcmp(Settings.Theory,'BH') % Buckingham model

    % Loose form of exp-C6 model
    if SigmaEpsilon

        % Default model parameters: all length-scale units are nm, energy scale units are kJ/mol
        DefMDSettings = Initialize_MD_Settings;
        DefMDSettings.Salt = Settings.Salt;
        [defBHMX,defBHMM,defBHXX] = BH_Potential_Parameters(DefMDSettings);

        % Input parameters
        r0_MM = Param(1); % nm
        r0_XX = Param(2); % nm

        if Additivity
            r0_MX = (r0_MM + r0_XX)/2; % nm

            epsilon_MM = Param(3); % kJ/mol
            epsilon_XX = Param(4); % kJ/mol
            epsilon_MX = sqrt(epsilon_MM*epsilon_XX); % kJ/mol

            gamma_MX = Param(5); % Unitless
            gamma_MM = gamma_MX; % Unitless
            gamma_XX = gamma_MX; % Unitless
            pidx = 5;
        else
            r0_MX = Param(3); % nm

            epsilon_MM = Param(4); % kJ/mol
            epsilon_XX = Param(5); % kJ/mol
            epsilon_MX = Param(6); % kJ/mol

            gamma_MM = Param(7); % Unitless
            gamma_XX = Param(8); % Unitless
            gamma_MX = Param(9); % Unitless
            pidx = 9;
        end
        
        if gamma_MX < 6
            epsilon_MX = -epsilon_MX;
        end
        if gamma_MM < 6
            epsilon_MM = -epsilon_MM;
        end
        if gamma_XX < 6
            epsilon_XX = -epsilon_XX;
        end

        % Convert to Condensed form
        alpha_MM = gamma_MM/r0_MM;
        alpha_XX = gamma_XX/r0_XX;
        alpha_MX = gamma_MX/r0_MX;

        B_MM = 6*epsilon_MM*exp(gamma_MM)/(gamma_MM - 6);
        B_XX = 6*epsilon_XX*exp(gamma_XX)/(gamma_XX - 6);
        B_MX = 6*epsilon_MX*exp(gamma_MX)/(gamma_MX - 6);

        C_MM = epsilon_MM*gamma_MM*(r0_MM^6)/(gamma_MM - 6);
        C_XX = epsilon_XX*gamma_XX*(r0_XX^6)/(gamma_XX - 6);
        C_MX = epsilon_MX*gamma_MX*(r0_MX^6)/(gamma_MX - 6);

        % Convert to scaling w.r.t. default BH
        Settings.S.D.MM = C_MM/defBHMM.C;
        Settings.S.D.XX = C_XX/defBHXX.C;
        Settings.S.D.MX = C_MX/defBHMX.C;

        Settings.S.A.MM = alpha_MM/defBHMM.alpha;
        Settings.S.A.XX = alpha_XX/defBHXX.alpha;
        Settings.S.A.MX = alpha_MX/defBHMX.alpha;

        Settings.S.R.MM = B_MM/defBHMM.B;
        Settings.S.R.XX = B_XX/defBHXX.B;
        Settings.S.R.MX = B_MX/defBHMX.B;

        % Scaling Coulombic Charge
        if Fix_Charge
            Settings.S.Q = Q_value;
        else
            Settings.S.Q = Param(pidx+1);
        end

    % Tight form of exp-C6 model
    else

        if Additivity
            Settings.S.D.MM = Param(1);
            Settings.S.D.XX = Param(2);
            Settings.S.D.MX = sqrt(Settings.S.D.MM*Settings.S.D.XX);

            % Repulsion
            Settings.S.R.MM = Param(3);
            Settings.S.R.XX = Param(4);
            Settings.S.R.MX = sqrt(Settings.S.R.MM*Settings.S.R.XX);
            pidx = 4;
        else
            Settings.S.D.MM = Param(1);
            Settings.S.D.XX = Param(2);
            Settings.S.D.MX = Param(3);

            % Repulsion
            Settings.S.R.MM = Param(4);
            Settings.S.R.XX = Param(5);
            Settings.S.R.MX = Param(6);
            pidx = 6;
        end

        % Alpha (TF exponential steepness repulsive parameter)
        if ~Fix_Alpha
            Settings.S.A.MM = Param(pidx+1);
            Settings.S.A.XX = Param(pidx+2);
            Settings.S.A.MX = Param(pidx+3);
            pidx = pidx+3;
        end

        if ~Fix_Charge
            Settings.S.Q = Param(pidx+1);
            pidx = pidx+1;

            if Additional_MM_Disp
                Settings.S.D.MM = Settings.S.D.MM + Param(pidx+1);
            end
        else
            Settings.S.Q = Q_value;
            if Additional_MM_Disp
                Settings.S.D.MM = Settings.S.D.MM + Param(pidx+1);
            end
        end
    end
else % JC models
    % sigma/epsilon form (cast in terms of sigma/epsilon scaling internally)
    if SigmaEpsilon
        
        % Default JC params
        DefMDSettings = Initialize_MD_Settings;
        DefMDSettings.WaterModel = 'SPC/E';
        DefMDSettings.Salt = Settings.Salt;
        [MXParams,MMParams,XXParams] = JC_Potential_Parameters(DefMDSettings);
        def_S_MX = MXParams.sigma;
        def_E_MX = MXParams.epsilon;

        % D <-> sigma
        % R <-> epsilon
        if Additivity
            % Sigma scaling
            Settings.S.S.MM = Param(1)/MMParams.sigma;
            Settings.S.S.XX = Param(2)/XXParams.sigma;

            % Epsilon scaling
            Settings.S.E.MM = Param(3)/MMParams.epsilon;
            Settings.S.E.XX = Param(4)/XXParams.epsilon;

            Sigma_MX = (Param(1) + Param(2))/2;
            Epsilon_MX = sqrt(Param(3)*Param(4));

            Settings.S.S.MX = Sigma_MX/def_S_MX;
            Settings.S.E.MX = Epsilon_MX/def_E_MX;

            % Scaling Coulombic Charge
            if Fix_Charge
                Settings.S.Q = Q_value;

                if Additional_MM_Disp
                    Full_MM_Epsilon = Param(3) + Param(5);
                    Settings.S.E.MM = Full_MM_Epsilon/MMParams.epsilon;
                    pidx = 5;
                else
                    pidx = 4;
                end
            else
                Settings.S.Q = Param(5);

                if Additional_MM_Disp
                    Full_MM_Epsilon = Param(3) + Param(6);
                    Settings.S.E.MM = Full_MM_Epsilon/MMParams.epsilon;
                    pidx = 6;
                else
                    pidx = 5;
                end
            end
        else
            % Sigma scaling
            Settings.S.S.MM = Param(1)/MMParams.sigma;
            Settings.S.S.XX = Param(2)/XXParams.sigma;
            Settings.S.S.MX = Param(3)/def_S_MX;

            % Epsilon scaling
            Settings.S.E.MM = Param(4)/MMParams.epsilon;
            Settings.S.E.XX = Param(5)/XXParams.epsilon;
            Settings.S.E.MX = Param(6)/def_E_MX;

            % Scaling Coulombic Charge
            if Fix_Charge
                Settings.S.Q = Q_value;
                pidx = 6;
            else
                Settings.S.Q = Param(7);
                pidx = 7;
            end
        end
    % Scaled dispersion/repulsion form
    else
        if Additivity
            % Dispersion
            Settings.S.D.MM = Param(1);
            Settings.S.D.XX = Param(2);

            % Repulsion
            Settings.S.R.MM = Param(3);
            Settings.S.R.XX = Param(4);
            
            DefMDSettings = Initialize_MD_Settings;
            DefMDSettings.WaterModel = 'SPC/E';
            DefMDSettings.Salt = Settings.Salt;
            [MXParams,MMParams,XXParams] = JC_Potential_Parameters(DefMDSettings);
            
            % Unscaled
            MX_Epsilon = MXParams.epsilon;
            MX_Sigma   = MXParams.sigma;

            MX_R = 4*MX_Epsilon*MX_Sigma^12;
            MX_D = 4*MX_Epsilon*MX_Sigma^6;

            % Scaled
            MM_Epsilon = MMParams.epsilon*(Settings.S.D.MM^2)*(1/Settings.S.R.MM);
            MM_Sigma = MMParams.sigma*(1/(Settings.S.D.MM^(1/6)))*(Settings.S.R.MM^(1/6));

            XX_Epsilon = XXParams.epsilon*(Settings.S.D.XX^2)*(1/Settings.S.R.XX);
            XX_Sigma = XXParams.sigma*(1/(Settings.S.D.XX^(1/6)))*(Settings.S.R.XX^(1/6));

            MX_Epsilon = sqrt(MM_Epsilon*XX_Epsilon);
            MX_Sigma   = (MM_Sigma + XX_Sigma)/2;

            MX_R_scaled = 4*MX_Epsilon*MX_Sigma^12;
            MX_D_scaled = 4*MX_Epsilon*MX_Sigma^6;

            Settings.S.D.MX = MX_D_scaled/MX_D;
            Settings.S.R.MX = MX_R_scaled/MX_R;

            % Scaling Coulombic Charge
            if Fix_Charge
                Settings.S.Q = Q_value;

                if Additional_MM_Disp
                    Settings.S.D.MM = Settings.S.D.MM + Param(5);
                    pidx = 5;
                else
                    pidx = 4;
                end
            else
                Settings.S.Q = Param(5);

                if Additional_MM_Disp
                    Settings.S.D.MM = Settings.S.D.MM + Param(6);
                    pidx = 6;
                else
                    pidx = 5;
                end
            end
        else
            % Dispersion
            Settings.S.D.MM = Param(1);
            Settings.S.D.XX = Param(2);
            Settings.S.D.MX = Param(3);

            % Repulsion
            Settings.S.R.MM = Param(4);
            Settings.S.R.XX = Param(5);
            Settings.S.R.MX = Param(6);

            % Scaling Coulombic Charge
            if Fix_Charge
                Settings.S.Q = Q_value;
                pidx = 6;
            else
                Settings.S.Q = Param(7);
                pidx = 7;
            end
        end
    end
end

% Perturb the potential with Gaussians
if isempty(Additional_GAdjust)
    Settings.GAdjust_MX = [0 0 1];  %[0 0 1];
    Settings.GAdjust_MM = [0 0 1];
    Settings.GAdjust_XX = [0 0 1];
else
    Settings.GAdjust_MX = [0 0 1];
    Settings.GAdjust_MM = [0 0 1];
    Settings.GAdjust_XX = [0 0 1];
    mx = 1;
    mm = 1;
    xx = 1;
    for idx = 1:length(Additional_GAdjust)
        int = [Additional_GAdjust{idx} '_' num2str(idx)];
        switch int
            case ['MM' '_' num2str(idx)]
                Settings.GAdjust_MM(mm,1) = Param(pidx + 1);
                Settings.GAdjust_MM(mm,2) = Param(pidx + 2);
                Settings.GAdjust_MM(mm,3) = Param(pidx + 3);
                mm = mm+1;
            case ['XX' '_' num2str(idx)]
                Settings.GAdjust_XX(xx,1) = Param(pidx + 1);
                Settings.GAdjust_XX(xx,2) = Param(pidx + 2);
                Settings.GAdjust_XX(xx,3) = Param(pidx + 3);
                xx = xx+1;
            case ['MX' '_' num2str(idx)]
                Settings.GAdjust_MX(mx,1) = Param(pidx + 1);
                Settings.GAdjust_MX(mx,2) = Param(pidx + 2);
                Settings.GAdjust_MX(mx,3) = Param(pidx + 3);
                mx = mx+1;
        end
        pidx = pidx + 3;
    end
end

% Additional functions
ints = {'MM' 'XX' 'MX'};
for idx = 1:length(ints)
    int = ints{idx};
    
    if Additional_Function.(int).N >= 0
        Settings.S.N.(int).Value = Additional_Function.(int).N;
        Settings.S.N.(int).Scale = Param(pidx + 1);
        pidx = pidx + 1;
    end
end

end