function Settings = Potential_Scaling(Settings,Param)

% Defaults
if ~isfield(Settings,'Q_value')
    Settings.Q_value = 1;
end
Settings.S = Init_Scaling_Object;

% Potential Scaling
if istable(Param) || isstruct(Param)
    if strcmp(Settings.Theory,'TF')
        
        % Loose form of exp-C6-C8 model
        if Settings.SigmaEpsilon
                                                    
            % Default model parameters: all length-scale units are nm,
            % energy scale units are kJ/mol
            
            [defTFMX,defTFMM,defTFXX] = TF_Potential_Parameters(Settings.Salt,Init_Scaling_Object);
            
            % Input parameters
            r0_MM = Param.r0_MM; % nm
            r0_XX = Param.r0_XX; % nm
            
            epsilon_MM = Param.epsilon_MM; % kJ/mol
            epsilon_XX = Param.epsilon_XX; % kJ/mol
            
            gamma_MX = Param.gamma_MX; % Unitless
            
            if Settings.Additivity
                r0_MX = (r0_MM + r0_XX)/2; % nm
                epsilon_MX = sqrt(epsilon_MM*epsilon_XX); % kJ/mol
                gamma_MM = gamma_MX; % Unitless
                gamma_XX = gamma_MX; % Unitless
            else
                r0_MX = Param.r0_MX; % nm
                epsilon_MX = Param.epsilon_MX; % kJ/mol
                gamma_MM = Param.gamma_MM; % Unitless
                gamma_XX = Param.gamma_XX; % Unitless
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
            
            % Convert to scaling w.r.t. TF
            Settings.S.A.MM = alpha_MM/defTFMM.alpha;
            Settings.S.A.XX = alpha_XX/defTFXX.alpha;
            Settings.S.A.MX = alpha_MX/defTFMX.alpha;
            
            Settings.S.R.MM = B_MM/defTFMM.B;
            Settings.S.R.XX = B_XX/defTFXX.B;
            Settings.S.R.MX = B_MX/defTFMX.B;
            
            Settings.S.D6D.MM = C_MM/defTFMM.C;
            Settings.S.D6D.XX = C_XX/defTFXX.C;
            Settings.S.D6D.MX = C_MX/defTFMX.C;
            
            Settings.S.D8D.MM = D_MM/defTFMM.D;
            Settings.S.D8D.XX = D_XX/defTFXX.D;
            Settings.S.D8D.MX = D_MX/defTFMX.D;
        
            % Scaling Coulombic Charge
            if Settings.Fix_Charge
                Settings.S.Q = Settings.Q_value;
            else
                Settings.S.Q = Param.SQ;
            end
            
        % Tight form of exp-C6-C8 model
        else
            
            % 1/R6 Dispersion (TF only)
            Settings.S.D6D.MM = Param.SD6MM;
            Settings.S.D6D.XX = Param.SD6XX;
            Settings.S.D6D.MX = Param.SD6MX;

            % 1/R8 Dispersion (TF only)
            if Settings.Fix_C8
                % Calculate value of C8 using recursive relations

                % Default TF params: C in units of kJ/mol nm^6, D in units of kJ/mol nm^8
                [TF_MX,TF_MM,TF_XX] = TF_Potential_Parameters(Settings.Salt,Init_Scaling_Object);

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

                % Update the scaling
                Settings.S.D8D.MM = C8.MM/TF_MM.D;
                Settings.S.D8D.XX = C8.XX/TF_XX.D;
                Settings.S.D8D.MX = C8.MX/TF_MX.D;
            else
                Settings.S.D8D.MM = Param.SD8MM;
                Settings.S.D8D.XX = Param.SD8XX;
                Settings.S.D8D.MX = Param.SD8MX;
            end

            % Alpha (TF exponential steepness repulsive parameter)
            if ~Settings.Fix_Alpha
                Settings.S.A.MM = Param.SAMM;
                Settings.S.A.XX = Param.SAXX;
                Settings.S.A.MX = Param.SAMX;
            end

            % Repulsive wall prefactor
            Settings.S.R.MM = Param.SRMM;
            Settings.S.R.XX = Param.SRXX;
            Settings.S.R.MX = Param.SRMX;

            % Scaling Coulombic Charge
            if Settings.Fix_Charge
                Settings.S.Q = Settings.Q_value;
            else
                Settings.S.Q = Param.SQ;
            end
        end
        
    elseif strcmp(Settings.Theory,'BH')
        
        % Loose form of exp-C6 model
        if Settings.SigmaEpsilon
        
            % Default model parameters: all length-scale units are nm, energy scale units are kJ/mol
            PotSettings = Initialize_MD_Settings;
            PotSettings.Salt = Settings.Salt;
            PotSettings.S = Init_Scaling_Object;
            [defBHMX,defBHMM,defBHXX] = BH_Potential_Parameters(PotSettings);
            
            % Input parameters
            r0_MM = Param.r0_MM; % nm
            r0_XX = Param.r0_XX; % nm
            
            epsilon_MM = Param.epsilon_MM; % kJ/mol
            epsilon_XX = Param.epsilon_XX; % kJ/mol
            
            gamma_MX = Param.gamma_MX; % Unitless
            
            if Settings.Additivity
                r0_MX = (r0_MM + r0_XX)/2; % nm
                epsilon_MX = sqrt(epsilon_MM*epsilon_XX); % kJ/mol
                gamma_MM = gamma_MX; % Unitless
                gamma_XX = gamma_MX; % Unitless
            else
                r0_MX = Param.r0_MX; % nm
                epsilon_MX = Param.epsilon_MX; % kJ/mol
                gamma_MM = Param.gamma_MM; % Unitless
                gamma_XX = Param.gamma_XX; % Unitless
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
            Settings.S.A.MM = alpha_MM/defBHMM.alpha;
            Settings.S.A.XX = alpha_XX/defBHXX.alpha;
            Settings.S.A.MX = alpha_MX/defBHMX.alpha;
            
            Settings.S.R.MM = B_MM/defBHMM.B;
            Settings.S.R.XX = B_XX/defBHXX.B;
            Settings.S.R.MX = B_MX/defBHMX.B;
            
            Settings.S.D.MM = C_MM/defBHMM.C;
            Settings.S.D.XX = C_XX/defBHXX.C;
            Settings.S.D.MX = C_MX/defBHMX.C;
        
            % Scaling Coulombic Charge
            if Settings.Fix_Charge
                Settings.S.Q = Settings.Q_value;
            else
                Settings.S.Q = Param.SQ;
            end
            
        % Tight form of the exp-C6 model
        else 
            % Dispersion
            Settings.S.D.MM = Param.SDMM;
            Settings.S.D.XX = Param.SDXX;

            % Repulsion
            Settings.S.R.MM = Param.SRMM;
            Settings.S.R.XX = Param.SRXX;

            if Settings.Additivity
                Settings.S.D.MX = sqrt(Settings.S.D.MM*Settings.S.D.XX);
                Settings.S.R.MX = sqrt(Settings.S.R.MM*Settings.S.R.XX);

                if Settings.Additional_MM_Disp
                    Settings.S.D.MM = Settings.S.D.MM + Param.SDMM2;
                end
            else
                Settings.S.D.MX = Param.SDMX;
                Settings.S.R.MX = Param.SRMX;
            end
            % Alpha (exponential steepness repulsive parameter)
            if ~Settings.Fix_Alpha
                Settings.S.A.MM = Param.SAMM;
                Settings.S.A.XX = Param.SAXX;
                Settings.S.A.MX = Param.SAMX;
            end

            % Scaling Coulombic Charge
            if Settings.Fix_Charge
                Settings.S.Q = Settings.Q_value;
            else
                Settings.S.Q = Param.SQ;
            end
        end
        
    elseif strcmp(Settings.Theory,'JC') % JC models
        
        % sigma/epsilon form (cast in terms of sigma/epsilon scaling internally)
        if Settings.SigmaEpsilon
            PotSettings = Initialize_MD_Settings;
            PotSettings.Salt = Settings.Salt;
            PotSettings.S = Init_Scaling_Object;
            [MXParams,MMParams,XXParams] = JC_Potential_Parameters(PotSettings);
            
            % Sigma scaling
            Settings.S.S.MM = Param.Sigma_MM/MMParams.sigma;
            Settings.S.S.XX = Param.Sigma_XX/XXParams.sigma;
            
            % Epsilon scaling
            Settings.S.E.MM = Param.Epsilon_MM/MMParams.epsilon;
            Settings.S.E.XX = Param.Epsilon_XX/XXParams.epsilon;
            
            % Default MX params
            def_S_MX = MXParams.sigma;
            def_E_MX = sMXParams.epsilon;
            
            if Settings.Additivity
                Sigma_MX = (Param.Sigma_MM + Param.Sigma_XX)/2;
                Epsilon_MX = sqrt(Param.Epsilon_MM*Param.Epsilon_XX);
                
                Settings.S.S.MX = Sigma_MX/def_S_MX;
                Settings.S.E.MX = Epsilon_MX/def_E_MX;
                
                if Settings.Additional_MM_Disp
                    Full_MM_Epsilon = Param.Epsilon_MM + Param.Epsilon_MM2;
                    Settings.S.E.MM = Full_MM_Epsilon/MMParams.epsilon;
                end
            else
                Settings.S.S.MX = Param.Sigma_MX/def_S_MX;
                Settings.S.E.MX = Param.Epsilon_MX/def_E_MX;
            end
        
        % Scaled dispersion/repulsion form
        else
            % Dispersion
            Settings.S.D.MM = Param.SDMM;
            Settings.S.D.XX = Param.SDXX;

            % Repulsion
            Settings.S.R.MM = Param.SRMM;
            Settings.S.R.XX = Param.SRXX;

            if Settings.Additivity
                PotSettings = Initialize_MD_Settings;
                PotSettings.Salt = Settings.Salt;
                PotSettings.S = Init_Scaling_Object;
                [MXParams,MMParams,XXParams] = JC_Potential_Parameters(PotSettings);

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

                if Settings.Additional_MM_Disp
                    Settings.S.D.MM = Settings.S.D.MM + Param.SDMM2;
                end
            else
                Settings.S.D.MX = Param.SDMX;
                Settings.S.R.MX = Param.SRMX;
            end
        end

        % Scaling Coulombic Charge
        if Settings.Fix_Charge
            Settings.S.Q = Settings.Q_value;
        else
            Settings.S.Q = Param.SQ;
        end
    end
    
    % Perturb the potential with Gaussians
    if ~isempty(Settings.Additional_GAdjust)
        mx = 1;
        mm = 1;
        xx = 1;
        for idx = 1:length(Settings.Additional_GAdjust)
            int = [Settings.Additional_GAdjust{idx} '_' num2str(idx)];
            switch int
                case ['MM' '_' num2str(idx)]
                    Settings.GAdjust_MM(mm,1) = Param.(['GA_' int]);
                    Settings.GAdjust_MM(mm,2) = Param.(['GB_' int]);
                    Settings.GAdjust_MM(mm,3) = Param.(['GC_' int]);
                    mm = mm+1;
                case ['XX' '_' num2str(idx)]
                    Settings.GAdjust_XX(xx,1) = Param.(['GA_' int]);
                    Settings.GAdjust_XX(xx,2) = Param.(['GB_' int]);
                    Settings.GAdjust_XX(xx,3) = Param.(['GC_' int]);
                    xx = xx+1;
                case ['MX' '_' num2str(idx)]
                    Settings.GAdjust_MX(mx,1) = Param.(['GA_' int]);
                    Settings.GAdjust_MX(mx,2) = Param.(['GB_' int]);
                    Settings.GAdjust_MX(mx,3) = Param.(['GC_' int]);
                    mx = mx+1;
            end
        end
    end
    
    % Additional functions
    if Settings.Additional_Function.MM.N >= 0 && ~isempty(Settings.Additional_Function.MM.Range)
        Settings.S.N.MM.Scale = Param.SNMM;
        Settings.S.N.MM.Value = Settings.Additional_Function.MM.N;
    end
    if Settings.Additional_Function.XX.N >= 0 && ~isempty(Settings.Additional_Function.XX.Range)
        Settings.S.N.XX.Scale = Param.SNXX;
        Settings.S.N.XX.Value = Settings.Additional_Function.XX.N;
    end
    if Settings.Additional_Function.MX.N >= 0 && ~isempty(Settings.Additional_Function.MX.Range)
        Settings.S.N.MX.Scale = Param.SNMX;
        Settings.S.N.MX.Value = Settings.Additional_Function.MX.N;
    end
    
% When not in table form
else
    if strcmp(Settings.Theory,'TF')
        
        % Loose form of exp-C6-C8 model
        if Settings.SigmaEpsilon
            
            % Default model parameters: all length-scale units are nm,
            % energy scale units are kJ/mol
            PotSettings = Initialize_MD_Settings;
            PotSettings.Salt = Settings.Salt;
            PotSettings.S = Init_Scaling_Object;
            [defTFMX,defTFMM,defTFXX] = TF_Potential_Parameters(PotSettings);
            
            % Input parameters
            r0_MM = Param(1); % nm
            r0_XX = Param(2); % nm
           
            if Settings.Additivity
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
            
            % Convert to scaling w.r.t. TF
            Settings.S.A.MM = alpha_MM/defTFMM.alpha;
            Settings.S.A.XX = alpha_XX/defTFXX.alpha;
            Settings.S.A.MX = alpha_MX/defTFMX.alpha;
            
            Settings.S.R.MM = B_MM/defTFMM.B;
            Settings.S.R.XX = B_XX/defTFXX.B;
            Settings.S.R.MX = B_MX/defTFMX.B;
            
            Settings.S.D6D.MM = C_MM/defTFMM.C;
            Settings.S.D6D.XX = C_XX/defTFXX.C;
            Settings.S.D6D.MX = C_MX/defTFMX.C;
            
            Settings.S.D8D.MM = D_MM/defTFMM.D;
            Settings.S.D8D.XX = D_XX/defTFXX.D;
            Settings.S.D8D.MX = D_MX/defTFMX.D;
        
            % Scaling Coulombic Charge
            if Settings.Fix_Charge
                Settings.S.Q = Settings.Q_value;
            else
                Settings.S.Q = Param(pidx+1);
                pidx = pidx+1;
            end
            
        % Tight form of exp-C6-C8 model
        else
            % 1/R6 Dispersion (TF only)
            Settings.S.D6D.MM = Param(1);
            Settings.S.D6D.XX = Param(2);
            Settings.S.D6D.MX = Param(3);

            % 1/R8 Dispersion (TF only)
            if Settings.Fix_C8
                % Calculate value of C8 using recursive relations            
                % Default TF params: C in units of kJ/mol nm^6, D in units of kJ/mol nm^8
                PotSettings = Initialize_MD_Settings;
                PotSettings.Salt = Settings.Salt;
                PotSettings.S = Init_Scaling_Object;
                [TF_MX,TF_MM,TF_XX] = TF_Potential_Parameters(PotSettings);

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

                % Update the scaling
                Settings.S.D8D.MM = C8.MM/TF_MM.D;
                Settings.S.D8D.XX = C8.XX/TF_XX.D;
                Settings.S.D8D.MX = C8.MX/TF_MX.D;
                pidx = 3;
            else
                Settings.S.D8D.MM = Param(4);
                Settings.S.D8D.XX = Param(5);
                Settings.S.D8D.MX = Param(6);
                pidx = 6;
            end

            % Alpha (TF exponential steepness repulsive parameter)
            if ~Settings.Fix_Alpha
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
            if Settings.Fix_Charge
                Settings.S.Q = Settings.Q_value;
            else
                Settings.S.Q = Param(pidx+1);
                pidx = pidx+1;
            end
        end
        
    elseif strcmp(Settings.Theory,'BH')
        
        % Loose form of exp-C6 model
        if Settings.SigmaEpsilon
        
            % Default model parameters: all length-scale units are nm, energy scale units are kJ/mol
            PotSettings = Initialize_MD_Settings;
            PotSettings.Salt = Settings.Salt;
            PotSettings.S = Init_Scaling_Object;
            [defBHMX,defBHMM,defBHXX] = BH_Potential_Parameters(PotSettings);
            
            % Input parameters
            r0_MM = Param(1); % nm
            r0_XX = Param(2); % nm
            
            if Settings.Additivity
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
            Settings.S.A.MM = alpha_MM/defBHMM.alpha;
            Settings.S.A.XX = alpha_XX/defBHXX.alpha;
            Settings.S.A.MX = alpha_MX/defBHMX.alpha;
            
            Settings.S.R.MM = B_MM/defBHMM.B;
            Settings.S.R.XX = B_XX/defBHXX.B;
            Settings.S.R.MX = B_MX/defBHMX.B;
            
            Settings.S.D.MM = C_MM/defBHMM.C;
            Settings.S.D.XX = C_XX/defBHXX.C;
            Settings.S.D.MX = C_MX/defBHMX.C;
        
            % Scaling Coulombic Charge
            if Settings.Fix_Charge
                Settings.S.Q = Settings.Q_value;
            else
                Settings.S.Q = Param(pidx+1);
                pidx = pidx+1;
            end
            
        % Tight form of the exp-C6 model
        else 
            
            if Settings.Additivity

                % Dispersion
                Settings.S.D.MM = Param(1);
                Settings.S.D.XX = Param(2);
                Settings.S.D.MX = sqrt(Settings.S.D.MM*Settings.S.D.XX);

                % Repulsion
                Settings.S.R.MM = Param(3);
                Settings.S.R.XX = Param(4);
                Settings.S.R.MX = sqrt(Settings.S.R.MM*Settings.S.R.XX);
                pidx = 4;

                if ~Settings.Fix_Alpha
                    Settings.S.A.MM = Param(pidx+1);
                    Settings.S.A.XX = Param(pidx+2);
                    Settings.S.A.MX = Param(pidx+3);
                    pidx = pidx+3;
                end

                % Scaling Coulombic Charge
                if Settings.Fix_Charge
                    Settings.S.Q = Settings.Q_value;

                    if Settings.Additional_MM_Disp
                        Settings.S.D.MM = Settings.S.D.MM + Param(pidx+1);
                        pidx = pidx+1;
                    end
                else
                    Settings.S.Q = Param(pidx+1);
                    pidx = pidx+1;

                    if Settings.Additional_MM_Disp
                        Settings.S.D.MM = Settings.S.D.MM + Param(pidx+1);
                        pidx = pidx+1;
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
                pidx = 6;

                if ~Settings.Fix_Alpha
                    Settings.S.A.MM = Param(pidx+1);
                    Settings.S.A.XX = Param(pidx+2);
                    Settings.S.A.MX = Param(pidx+3);
                    pidx = pidx+3;
                end            

                % Scaling Coulombic Charge
                if Settings.Fix_Charge
                    Settings.S.Q = Settings.Q_value;
                else
                    Settings.S.Q = Param(pidx+1);
                    pidx = pidx+1;
                end
            end
        end
        
    else % JC models
        
        % sigma/epsilon form (cast in terms of sigma/epsilon scaling internally)
        if Settings.SigmaEpsilon
            
            % Default JC params
            PotSettings = Initialize_MD_Settings;
            PotSettings.Salt = Settings.Salt;
            PotSettings.S = Init_Scaling_Object;
            [MXParams,MMParams,XXParams] = JC_Potential_Parameters(PotSettings);
            def_S_MX = MXParams.sigma;
            def_E_MX = MXParams.epsilon;
            
            % D <-> sigma
            % R <-> epsilon
            if Settings.Additivity
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
                if Settings.Fix_Charge
                    Settings.S.Q = Settings.Q_value;

                    if Settings.Additional_MM_Disp
                        Full_MM_Epsilon = Param(3) + Param(5);
                        Settings.S.E.MM = Full_MM_Epsilon/MMParams.epsilon;
                        pidx = 5;
                    else
                        pidx = 4;
                    end
                else
                    Settings.S.Q = Param(5);

                    if Settings.Additional_MM_Disp
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
                if Settings.Fix_Charge
                    Settings.S.Q = Settings.Q_value;
                    pidx = 6;
                else
                    Settings.S.Q = Param(7);
                    pidx = 7;
                end
            end
        % Scaled dispersion/repulsion form
        else
            if Settings.Additivity

                % Dispersion
                Settings.S.D.MM = Param(1);
                Settings.S.D.XX = Param(2);

                % Repulsion
                Settings.S.R.MM = Param(3);
                Settings.S.R.XX = Param(4);
                
                PotSettings = Initialize_MD_Settings;
                PotSettings.Salt = Settings.Salt;
                PotSettings.S = Init_Scaling_Object;
                [MXParams,MMParams,XXParams] = JC_Potential_Parameters(PotSettings);
                
                % Unscaled
                MX_Sigma = MXParams.sigma;
                MX_Epsilon = MXParams.epsilon;

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
                if Settings.Fix_Charge
                    Settings.S.Q = Settings.Q_value;

                    if Settings.Additional_MM_Disp
                        Settings.S.D.MM = Settings.S.D.MM + Param(5);
                        pidx = 5;
                    else
                        pidx = 4;
                    end
                else
                    Settings.S.Q = Param(5);

                    if Settings.Additional_MM_Disp
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
                if Settings.Fix_Charge
                    Settings.S.Q = Settings.Q_value;
                    pidx = 6;
                else
                    Settings.S.Q = Param(7);
                    pidx = 7;
                end
            end
        end
    end
    
    % Perturb the potential with Gaussians
    if isempty(Settings.Additional_GAdjust)
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
        
        for idx = 1:length(Settings.Additional_GAdjust)
            int = [Settings.Additional_GAdjust{idx} '_' num2str(idx)];
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
    if Settings.Additional_Function.MM.N >= 0 && ~isempty(Settings.Additional_Function.MM.Range)
        pidx = pidx+1;
        Settings.S.N.MM.Scale = Param(pidx);
        Settings.S.N.MM.Value = Settings.Additional_Function.MM.N;
    end
    if Settings.Additional_Function.XX.N >= 0 && ~isempty(Settings.Additional_Function.XX.Range)
        pidx = pidx+1;
        Settings.S.N.XX.Scale = Param(pidx);
        Settings.S.N.XX.Value = Settings.Additional_Function.XX.N;
    end
    if Settings.Additional_Function.MX.N >= 0 && ~isempty(Settings.Additional_Function.MX.Range)
        pidx = pidx+1;
        Settings.S.N.MX.Scale = Param(pidx);
        Settings.S.N.MX.Value = Settings.Additional_Function.MX.N;
    end
    
end

end