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
        
        % Outputs
        Settings.S.S.MM = r0_MM;
        Settings.S.S.XX = r0_XX;
        Settings.S.S.MX = r0_MX;
        
        Settings.S.E.MM = epsilon_MM;
        Settings.S.E.XX = epsilon_XX;
        Settings.S.E.MX = epsilon_MX;
        
        Settings.S.G.MM = gamma_MM;
        Settings.S.G.XX = gamma_XX;
        Settings.S.G.MX = gamma_MX;

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
            C8_MM = 3.0.*(Settings.S.D6D.MM./c6units).*sqrt_Q.(Settings.Metal).*sqrt_Q.(Settings.Metal).*c8units; % in kJ/mol nm^8
            C8_XX = 3.0.*(Settings.S.D6D.XX./c6units).*sqrt_Q.(Settings.Halide).*sqrt_Q.(Settings.Halide).*c8units; % in kJ/mol nm^8
            C8_MX = 3.0.*(Settings.S.D6D.MX./c6units).*sqrt_Q.(Settings.Metal).*sqrt_Q.(Settings.Halide).*c8units; % in kJ/mol nm^8

            % Update the scaling
            Settings.S.D8D.MM = C8_MM;
            Settings.S.D8D.XX = C8_XX;
            Settings.S.D8D.MX = C8_MX;
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
case {'BH' 'BD' 'BE'}

    % Loose form of exp-C6 model
    if Settings.SigmaEpsilon
        
        % Input parameters
        r0_MM = Param.r0_MM; % nm
        r0_XX = Param.r0_XX; % nm
        
        epsilon_MM = Param.epsilon_MM; % kJ/mol
        epsilon_XX = Param.epsilon_XX; % kJ/mol
        
        if Settings.Additivity
            switch lower(Settings.Comb_rule)
                case 'lorentz-berthelot'
                    r0_MX = (r0_MM + r0_XX)./2; % nm
                    epsilon_MX = sqrt(epsilon_MM.*epsilon_XX); % kJ/mol
                    gamma_MX = Param.gamma_MX; % Unitless
                    gamma_MM = gamma_MX; % Unitless
                    gamma_XX = gamma_MX; % Unitless
                    epsilon_MM(gamma_MM < 6) = -abs(epsilon_MM(gamma_MM < 6));
                    epsilon_XX(gamma_XX < 6) = -abs(epsilon_XX(gamma_XX < 6));
                    epsilon_MX(gamma_MX < 6) = -abs(epsilon_MX(gamma_MX < 6));
                case 'hogervorst'
                    gamma_MM = Param.gamma_MM; % Unitless
                    gamma_XX = Param.gamma_XX; % Unitless
                    gamma_MX = (gamma_MM + gamma_XX)./2;

                    epsilon_MX = 2.*epsilon_MM.*epsilon_XX./(epsilon_MM + epsilon_XX);

                    epsilon_MM(gamma_MM < 6) = -abs(epsilon_MM(gamma_MM < 6));
                    epsilon_XX(gamma_XX < 6) = -abs(epsilon_XX(gamma_XX < 6));
                    epsilon_MX(gamma_MX < 6) = -abs(epsilon_MX(gamma_MX < 6));

                    r0_MX = ( sqrt( ( epsilon_MM.*epsilon_XX.*gamma_MM.*gamma_XX.*(r0_MM.*r0_XX).^6 )...
                        ./((gamma_MM - 6).*(gamma_XX - 6)) ).*(gamma_MX - 6)./(epsilon_MX.*gamma_MX) ).^(1/6);

                case {'kong' 'gromacs'}
                    gamma_MM = Param.gamma_MM; % Unitless
                    gamma_XX = Param.gamma_XX; % Unitless

                    epsilon_MM(gamma_MM < 6) = -abs(epsilon_MM(gamma_MM < 6));
                    epsilon_XX(gamma_XX < 6) = -abs(epsilon_XX(gamma_XX < 6));

                    k_MM = 1 ./ ( gamma_MM - 6 );
                    k_XX = 1 ./ ( gamma_XX - 6 );

                    A_MM = 6.*epsilon_MM.*k_MM.*exp(gamma_MM); % prefactor
                    A_XX = 6.*epsilon_XX.*k_XX.*exp(gamma_XX);

                    B_MM = gamma_MM./r0_MM; % exponent
                    B_XX = gamma_XX./r0_XX;

                    C_MM = epsilon_MM.*gamma_MM.*k_MM.*(r0_MM.^6); % dispersion
                    C_XX = epsilon_XX.*gamma_XX.*k_XX.*(r0_XX.^6);

                    switch lower(Settings.Comb_rule)
                        case 'kong'
                            A_MX = (1/2).*( A_MM.*(A_MM.*B_MM./(A_XX.*B_XX)).^(-B_MM./(B_MM + B_XX)) + ...
                                            A_XX.*(A_XX.*B_XX./(A_MM.*B_MM)).^(-B_XX./(B_MM + B_XX)) );
                            B_MX = 2.*B_MM.*B_XX./(B_MM + B_XX);
                            C_MX = sqrt(C_MM.*C_XX);
                        case 'gromacs'
                            A_MX = sqrt(A_MM.*A_XX);
                            B_MX = 2./( (1./B_MM) + (1./B_XX) );
                            C_MX = sqrt(C_MM.*C_XX);
                    end

                    % Convert back to gamma/epsilon/r0
                    gamma_MX   = -7.*lambertw((-1/7).*(6.*C_MX.*(B_MX.^6)./A_MX).^(1./7));
                    r0_MX      = gamma_MX./B_MX;
                    epsilon_MX = C_MX.*(gamma_MX - 6)./(gamma_MX.*(r0_MX.^6));
                    epsilon_MX(gamma_MX < 6) = -abs(epsilon_MX(gamma_MX < 6));
            end
        else
            r0_MX = Param.r0_MX; % nm
            epsilon_MX = Param.epsilon_MX; % kJ/mol
            gamma_MM = Param.gamma_MM; % Unitless
            gamma_XX = Param.gamma_XX; % Unitless
            epsilon_MX(gamma_MX < 6) = -epsilon_MX(gamma_MX < 6);
            epsilon_MM(gamma_MM < 6) = -epsilon_MM(gamma_MM < 6);
            epsilon_XX(gamma_XX < 6) = -epsilon_XX(gamma_XX < 6);
        end
        
        % Outputs
        Settings.S.S.MM = r0_MM;
        Settings.S.S.XX = r0_XX;
        Settings.S.S.MX = r0_MX;
        
        Settings.S.E.MM = epsilon_MM;
        Settings.S.E.XX = epsilon_XX;
        Settings.S.E.MX = epsilon_MX;
        
        Settings.S.G.MM = gamma_MM;
        Settings.S.G.XX = gamma_XX;
        Settings.S.G.MX = gamma_MX;
        
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
case 'LJ'
    
    % sigma/epsilon form (cast in terms of sigma/epsilon scaling internally)
    if Settings.SigmaEpsilon

        % Sigma scaling
        Settings.S.S.MM = Param.sigma_MM;
        Settings.S.S.XX = Param.sigma_XX;

        % Epsilon scaling
        Settings.S.E.MM = Param.epsilon_MM;
        Settings.S.E.XX = Param.epsilon_XX;
        
        if Settings.Additivity
            Settings.S.S.MX = (Param.sigma_MM + Param.sigma_XX)./2;
            Settings.S.E.MX = sqrt(Param.epsilon_MM.*Param.epsilon_XX);

            if Settings.Additional_MM_Disp
                Full_MM_Epsilon = Param.epsilon_MM + Param.epsilon_MM2;
                Settings.S.E.MM = Full_MM_Epsilon;
            end
        else
            Settings.S.S.MX = Param.sigma_MX;
            Settings.S.E.MX = Param.epsilon_MX;
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
            
            MM_Epsilon = (Settings.S.D.MM.^2).*(1./Settings.S.R.MM);
            MM_Sigma = (1./(Settings.S.D.MM.^(1/6))).*(Settings.S.R.MM.^(1/6));

            XX_Epsilon = (Settings.S.D.XX.^2).*(1./Settings.S.R.XX);
            XX_Sigma = (1./(Settings.S.D.XX.^(1/6))).*(Settings.S.R.XX.^(1/6));

            MX_Epsilon = sqrt(MM_Epsilon.*XX_Epsilon);
            MX_Sigma   = (MM_Sigma + XX_Sigma)./2;
            
            Settings.S.D.MX = 4.*MX_Epsilon.*MX_Sigma.^6;
            Settings.S.R.MX = 4.*MX_Epsilon.*MX_Sigma.^12;

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
case 'BF'
    % Input parameters
    Settings.S.S.MM = Param.sigma_MM; % nm
    Settings.S.S.XX = Param.sigma_XX; % nm
    
    Settings.S.E.MM = Param.epsilon_MM; % kJ/mol
    Settings.S.E.XX = Param.epsilon_XX; % kJ/mol
        
    if Settings.Additivity       
        switch lower(Settings.Comb_rule)
            case 'lorentz-berthelot'
                Settings.S.G.MX = Param.gamma_MX; % Unitless
                Settings.S.S.MX = (Settings.S.S.MM + Settings.S.S.XX)./2; % nm
                Settings.S.E.MX = sqrt(Settings.S.E.MM.*Settings.S.E.XX); % kJ/mol
                Settings.S.G.MM = Settings.S.G.MX; % Unitless
                Settings.S.G.XX = Settings.S.G.MX; % Unitless
            case 'hogervorst'
                Settings.S.G.MM = Param.gamma_MM; % Unitless
                Settings.S.G.XX = Param.gamma_XX; % Unitless
                Settings.S.G.MX = (Settings.S.G.MM + Settings.S.G.XX)./2;

                Settings.S.E.MX = 2.*Settings.S.E.MM.*Settings.S.E.XX./(Settings.S.E.MM + Settings.S.E.XX);

                Settings.S.S.MX = ( sqrt( ( Settings.S.E.MM.*Settings.S.E.XX.*Settings.S.G.MM.*Settings.S.G.XX.*(Settings.S.S.MM.*Settings.S.S.XX).^6 )...
                    ./((Settings.S.G.MM - 6).*(Settings.S.G.XX - 6)) ).*(Settings.S.G.MX - 6)./(Settings.S.E.MX.*Settings.S.G.MX) ).^(1/6);
            case 'gromacs'
                Settings.S.G.MM = Param.gamma_MM; % Unitless
                Settings.S.G.XX = Param.gamma_XX; % Unitless
                Settings.S.G.MX = sqrt(Settings.S.G.MM.*Settings.S.G.XX);

                Settings.S.E.MX = sqrt(Settings.S.E.MM.*Settings.S.E.XX);
                Settings.S.S.MX = sqrt(Settings.S.S.MM.*Settings.S.S.XX);
        end
        
    else
        Settings.S.S.MX = Param.sigma_MX; % nm
        Settings.S.E.MX = Param.epsilon_MX; % kJ/mol
        Settings.S.G.MM = Param.gamma_MM; % Unitless
        Settings.S.G.XX = Param.gamma_XX; % Unitless
        Settings.S.G.MX = Param.gamma_MX;
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
        Settings.S.S.MM = Param.sigma_MM;
        Settings.S.S.XX = Param.sigma_XX;

        % Epsilon scaling
        Settings.S.E.MM = Param.epsilon_MM;
        Settings.S.E.XX = Param.epsilon_XX;
        
        if Settings.Additivity
            
            if ~Settings.Fix_Mie_n
                Settings.S.n.MX = Param.n_MX;
                Settings.S.n.MM = Param.n_MX;
                Settings.S.n.XX = Param.n_MX;
            end
            
            Sigma_MX = (Param.sigma_MM + Param.sigma_XX)./2;
            Epsilon_MX = sqrt(Param.epsilon_MM.*Param.epsilon_XX);

            Settings.S.S.MX = Sigma_MX;
            Settings.S.E.MX = Epsilon_MX;

            if Settings.Additional_MM_Disp
                Settings.S.E.MM = Param.epsilon_MM + Param.epsilon_MM2;
            end
        else
            Settings.S.S.MX = Param.sigma_MX;
            Settings.S.E.MX = Param.epsilon_MX;
            
            if ~Settings.Fix_Mie_n
                Settings.S.n.MX = Param.n_MX;
                Settings.S.n.MM = Param.n_MM;
                Settings.S.n.XX = Param.n_XX;
            end
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