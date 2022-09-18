function Settings = Get_Scaling_Params(Settings,Param)

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

PotSettings = Initialize_MD_Settings;
PotSettings.Salt = Settings.Salt;
[JC_MX,JC_MM,JC_XX] = JC_Potential_Parameters(PotSettings);
[BH_MX,BH_MM,BH_XX] = BH_Potential_Parameters(PotSettings);
[TF_MX,TF_MM,TF_XX] = TF_Potential_Parameters(PotSettings);

% Potential Scaling
switch Settings.Theory
case 'TF'

    % Loose form of exp-C6-C8 model
    if Settings.SigmaEpsilon

        % Default model parameters: all length-scale units are nm,
        % energy scale units are kJ/mol

        % Input parameters
        r0_MM = Param.r0_MM; % nm
        r0_XX = Param.r0_XX; % nm

        epsilon_MM = Param.epsilon_MM; % kJ/mol
        epsilon_XX = Param.epsilon_XX; % kJ/mol

        gamma_MX = Param.gamma_MX; % Unitless

        if Settings.Additivity
            r0_MX = (r0_MM + r0_XX)./2; % nm
            epsilon_MX = sqrt(epsilon_MM.*epsilon_XX); % kJ/mol
            gamma_MM = gamma_MX; % Unitless
            gamma_XX = gamma_MX; % Unitless
        else
            r0_MX = Param.r0_MX; % nm
            epsilon_MX = Param.epsilon_MX; % kJ/mol
            gamma_MM = Param.gamma_MM; % Unitless
            gamma_XX = Param.gamma_XX; % Unitless
        end
        
        epsilon_MX(gamma_MX > 48/7) = -epsilon_MX(gamma_MX > 48/7);
        epsilon_MM(gamma_MM > 48/7) = -epsilon_MM(gamma_MM > 48/7);
        epsilon_XX(gamma_XX > 48/7) = -epsilon_XX(gamma_XX > 48/7);
        
        % Convert to Condensed form
        alpha_MM = gamma_MM./r0_MM;
        alpha_XX = gamma_XX./r0_XX;
        alpha_MX = gamma_MX./r0_MX;

        B_MM = 48.*epsilon_MM.*exp(gamma_MM)./(48 - 7*gamma_MM);
        B_XX = 48.*epsilon_XX.*exp(gamma_XX)./(48 - 7*gamma_XX);
        B_MX = 48.*epsilon_MX.*exp(gamma_MX)./(48 - 7*gamma_MX);

        C_MM = 4.*epsilon_MM.*gamma_MM.*(r0_MM.^6)./(48 - 7.*gamma_MM);
        C_XX = 4.*epsilon_XX.*gamma_XX.*(r0_XX.^6)./(48 - 7.*gamma_XX);
        C_MX = 4.*epsilon_MX.*gamma_MX.*(r0_MX.^6)./(48 - 7.*gamma_MX);

        D_MM = 3.*epsilon_MM.*gamma_MM.*(r0_MM.^8)./(48 - 7.*gamma_MM);
        D_XX = 3.*epsilon_XX.*gamma_XX.*(r0_XX.^8)./(48 - 7.*gamma_XX);
        D_MX = 3.*epsilon_MX.*gamma_MX.*(r0_MX.^8)./(48 - 7.*gamma_MX);

        % Convert to scaling w.r.t. TF
        Settings.S.A.MM = alpha_MM./TF_MM.alpha;
        Settings.S.A.XX = alpha_XX./TF_XX.alpha;
        Settings.S.A.MX = alpha_MX./TF_MX.alpha;

        Settings.S.R.MM = B_MM./TF_MM.B;
        Settings.S.R.XX = B_XX./TF_XX.B;
        Settings.S.R.MX = B_MX./TF_MX.B;

        Settings.S.D6D.MM = C_MM./TF_MM.C;
        Settings.S.D6D.XX = C_XX./TF_XX.C;
        Settings.S.D6D.MX = C_MX./TF_MX.C;

        Settings.S.D8D.MM = D_MM./TF_MM.D;
        Settings.S.D8D.XX = D_XX./TF_XX.D;
        Settings.S.D8D.MX = D_MX./TF_MX.D;

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

            % Calculate Scaled C8 using recursion relation from D3 paper
            C8_MM = 3.0.*(Settings.S.D6D.MM.*TF_MM.C./c6units).*sqrt_Q.(Settings.Metal).*sqrt_Q.(Settings.Metal).*c8units; % in kJ/mol nm^8
            C8_XX = 3.0.*(Settings.S.D6D.XX.*TF_XX.C./c6units).*sqrt_Q.(Settings.Halide).*sqrt_Q.(Settings.Halide).*c8units; % in kJ/mol nm^8
            C8_MX = 3.0.*(Settings.S.D6D.MX.*TF_MX.C./c6units).*sqrt_Q.(Settings.Metal).*sqrt_Q.(Settings.Halide).*c8units; % in kJ/mol nm^8

            % Update the scaling
            Settings.S.D8D.MM = C8_MM./TF_MM.D;
            Settings.S.D8D.XX = C8_XX./TF_XX.D;
            Settings.S.D8D.MX = C8_MX./TF_MX.D;
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

case {'BH' 'BD'}

    % Loose form of exp-C6 model
    if Settings.SigmaEpsilon
        
        % Input parameters
        r0_MM = Param.r0_MM; % nm
        r0_XX = Param.r0_XX; % nm
        
        epsilon_MM = Param.epsilon_MM; % kJ/mol
        epsilon_XX = Param.epsilon_XX; % kJ/mol
        
        gamma_MX = Param.gamma_MX; % Unitless
        
        if Settings.Additivity
            r0_MX = (r0_MM + r0_XX)./2; % nm
            epsilon_MX = sqrt(epsilon_MM.*epsilon_XX); % kJ/mol
            gamma_MM = gamma_MX; % Unitless
            gamma_XX = gamma_MX; % Unitless
        else
            r0_MX = Param.r0_MX; % nm
            epsilon_MX = Param.epsilon_MX; % kJ/mol
            gamma_MM = Param.gamma_MM; % Unitless
            gamma_XX = Param.gamma_XX; % Unitless
        end
        
        epsilon_MX(gamma_MX < 6) = -epsilon_MX(gamma_MX < 6);
        epsilon_MM(gamma_MM < 6) = -epsilon_MM(gamma_MM < 6);
        epsilon_XX(gamma_XX < 6) = -epsilon_XX(gamma_XX < 6);
        
        % Convert to Condensed form
        alpha_MM = gamma_MM./r0_MM;
        alpha_XX = gamma_XX./r0_XX;
        alpha_MX = gamma_MX./r0_MX;
        
        B_MM = 6.*epsilon_MM.*exp(gamma_MM)./(gamma_MM - 6);
        B_XX = 6.*epsilon_XX.*exp(gamma_XX)./(gamma_XX - 6);
        B_MX = 6.*epsilon_MX.*exp(gamma_MX)./(gamma_MX - 6);
        
        C_MM = epsilon_MM.*gamma_MM.*(r0_MM.^6)./(gamma_MM - 6);
        C_XX = epsilon_XX.*gamma_XX.*(r0_XX.^6)./(gamma_XX - 6);
        C_MX = epsilon_MX.*gamma_MX.*(r0_MX.^6)./(gamma_MX - 6);
        
        % Convert to scaling w.r.t. default BH
        Settings.S.A.MM = alpha_MM./BH_MM.alpha;
        Settings.S.A.XX = alpha_XX./BH_XX.alpha;
        Settings.S.A.MX = alpha_MX./BH_MX.alpha;
        
        Settings.S.R.MM = B_MM./BH_MM.B;
        Settings.S.R.XX = B_XX./BH_XX.B;
        Settings.S.R.MX = B_MX./BH_MX.B;
        
        Settings.S.D.MM = C_MM./BH_MM.C;
        Settings.S.D.XX = C_XX./BH_XX.C;
        Settings.S.D.MX = C_MX./BH_MX.C;
        
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
            Settings.S.D.MX = sqrt(Settings.S.D.MM.*Settings.S.D.XX);
            Settings.S.R.MX = sqrt(Settings.S.R.MM.*Settings.S.R.XX);

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

case 'JC' % JC models

    % sigma/epsilon form (cast in terms of sigma/epsilon scaling internally)
    if Settings.SigmaEpsilon

        % Sigma scaling
        Settings.S.S.MM = Param.sigma_MM./JC_MM.sigma;
        Settings.S.S.XX = Param.sigma_XX./JC_XX.sigma;

        % Epsilon scaling
        Settings.S.E.MM = Param.epsilon_MM./JC_MM.epsilon;
        Settings.S.E.XX = Param.epsilon_XX./JC_XX.epsilon;

        % Default MX params
        def_S_MX = JC_MX.sigma;
        def_E_MX = JC_MX.epsilon;

        if Settings.Additivity
            Sigma_MX = (Param.sigma_MM + Param.sigma_XX)./2;
            Epsilon_MX = sqrt(Param.epsilon_MM.*Param.epsilon_XX);

            Settings.S.S.MX = Sigma_MX./def_S_MX;
            Settings.S.E.MX = Epsilon_MX./def_E_MX;

            if Settings.Additional_MM_Disp
                Full_MM_Epsilon = Param.epsilon_MM + Param.epsilon_MM2;
                Settings.S.E.MM = Full_MM_Epsilon./JC_MM.epsilon;
            end
        else
            Settings.S.S.MX = Param.sigma_MX./def_S_MX;
            Settings.S.E.MX = Param.epsilon_MX./def_E_MX;
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

            % Unscaled
            MX_Epsilon = JC_MX.epsilon;
            MX_Sigma   = JC_MX.sigma;

            MX_R = 4.*MX_Epsilon.*MX_Sigma.^12;
            MX_D = 4.*MX_Epsilon.*MX_Sigma.^6;

            % Scaled
            MM_Epsilon = JC_MM.epsilon.*(Settings.S.D.MM.^2).*(1./Settings.S.R.MM);
            MM_Sigma = JC_MM.sigma.*(1./(Settings.S.D.MM.^(1/6))).*(Settings.S.R.MM.^(1/6));

            XX_Epsilon = JC_XX.epsilon.*(Settings.S.D.XX.^2).*(1./Settings.S.R.XX);
            XX_Sigma = JC_XX.sigma.*(1./(Settings.S.D.XX.^(1/6))).*(Settings.S.R.XX.^(1/6));

            MX_Epsilon = sqrt(MM_Epsilon.*XX_Epsilon);
            MX_Sigma   = (MM_Sigma + XX_Sigma)./2;

            MX_R_scaled = 4.*MX_Epsilon.*MX_Sigma.^12;
            MX_D_scaled = 4.*MX_Epsilon.*MX_Sigma.^6;

            Settings.S.D.MX = MX_D_scaled./MX_D;
            Settings.S.R.MX = MX_R_scaled./MX_R;

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
    
case 'Mie'
    % sigma/epsilon form (cast in terms of sigma/epsilon scaling internally)
    if Settings.SigmaEpsilon

        % Sigma scaling
        Settings.S.S.MM = Param.sigma_MM./JC_MM.sigma;
        Settings.S.S.XX = Param.sigma_XX./JC_XX.sigma;

        % Epsilon scaling
        Settings.S.E.MM = Param.epsilon_MM./JC_MM.epsilon;
        Settings.S.E.XX = Param.epsilon_XX./JC_XX.epsilon;
        
        Settings.S.n.MX = Param.n_MX;
        
        % Default MX params
        def_S_MX = JC_MX.sigma;
        def_E_MX = JC_MX.epsilon;
        
        if Settings.Additivity
            
            Settings.S.n.MM = Param.n_MX;
            Settings.S.n.XX = Param.n_MX;
            
            Sigma_MX = (Param.sigma_MM + Param.sigma_XX)./2;
            Epsilon_MX = sqrt(Param.epsilon_MM.*Param.epsilon_XX);

            Settings.S.S.MX = Sigma_MX./def_S_MX;
            Settings.S.E.MX = Epsilon_MX./def_E_MX;

            if Settings.Additional_MM_Disp
                Full_MM_Epsilon = Param.epsilon_MM + Param.epsilon_MM2;
                Settings.S.E.MM = Full_MM_Epsilon./JC_MM.epsilon;
            end
        else
            Settings.S.S.MX = Param.sigma_MX./def_S_MX;
            Settings.S.E.MX = Param.epsilon_MX./def_E_MX;
            
            Settings.S.n.MM = Param.n_MM;
            Settings.S.n.XX = Param.n_XX;
        end

    % Scaled dispersion/repulsion form
    else
        error('Mie potential only available in sigma-epsilon form')
    end

    % Scaling Coulombic Charge
    if Settings.Fix_Charge
        Settings.S.Q = Settings.Q_value;
    else
        Settings.S.Q = Param.SQ;
    end
end

end