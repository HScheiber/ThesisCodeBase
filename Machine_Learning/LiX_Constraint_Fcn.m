function tf = LiX_Constraint_Fcn(Settings,Param,varargin)

% Optional inputs
p = inputParser;
p.FunctionName = 'LiX_Constraint_Fcn';
addOptional(p,'Output_Loss',false,@(x)validateattributes(x,{'logical'},{'nonempty'}))

parse(p,varargin{:});
Output_Loss = p.Results.Output_Loss;

Settings.Table_Length = 10; % nm
Settings.Table_StepSize = 0.01;

if ~isfolder(Settings.home)
    Settings.home = find_home;
end

% Conversion factors
Bohr_nm = 0.0529177; % a_0 - > Angstrom
c6conv = 1e-3/2625.4999/((0.052917726)^6); % J/mol nm^6 - > au (from sourcecode)
J_kJ = 1e-3; % J - > kJ
Ha_kJmol = 2625.4999; % Ha - > kJ/mol
c6units = (1/c6conv)*J_kJ; % au - > kJ/mol nm^6
c8units = (Ha_kJmol)*(Bohr_nm^8); % au - > kJ/mol nm^8

% Factors used to calculate C8
sqrt_Q.Li = 5.019869340000000;
sqrt_Q.Na = 6.585855360000000;
sqrt_Q.K  = 7.977627530000000;
sqrt_Q.Rb = 9.554616980000000;
sqrt_Q.Cs = 11.02204549000000;
sqrt_Q.F  = 2.388252500000000;
sqrt_Q.Cl = 3.729323560000000;
sqrt_Q.Br = 4.590896470000000;
sqrt_Q.I  = 5.533218150000000;

% Param will always be a table.
N_par = size(Param,1);
Loss = zeros(N_par,1);
ParNames = Param.Properties.VariableNames;

% Coupled or each salt independent?
[Salts,MultiMetal,MultiHalide] = Select_MultiSalt(Settings);
for salt_idx = 1:numel(Salts) % Loop through coupled salts
    Settings.Salt = Salts{salt_idx};
    [Settings.Metal,Settings.Halide] = Separate_Metal_Halide(Settings.Salt);
    
    if Settings.GaussianCharge
        Settings = Alexandria_Potential_Parameters(Settings,'Coulomb_Only',true); % Loads Gaussian charge parameters
    end
    
    if MultiMetal
        for idx = 1:numel(ParNames)
            if contains(ParNames{idx},[Settings.Metal Settings.Metal])
                pname = strrep(ParNames{idx},[Settings.Metal Settings.Metal],'MM');
                Param.(pname) = Param.(ParNames{idx});
            end
        end
    end
    if MultiHalide
        for idx = 1:numel(ParNames)
            if contains(ParNames{idx},[Settings.Halide Settings.Halide])
                pname = strrep(ParNames{idx},[Settings.Halide Settings.Halide],'XX');
                Param.(pname) = Param.(ParNames{idx});
            end
        end
    end
    
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
                        gamma_MX = Param.gamma_MX; % Unitless
                        r0_MX = (r0_MM + r0_XX)./2; % nm
                        epsilon_MX = sqrt(epsilon_MM.*epsilon_XX); % kJ/mol
                        gamma_MM = gamma_MX; % Unitless
                        gamma_XX = gamma_MX; % Unitless
                        epsilon_MM(gamma_MM < 6) = -abs(epsilon_MM(gamma_MM < 6));
                        epsilon_XX(gamma_XX < 6) = -abs(epsilon_XX(gamma_XX < 6));
                        epsilon_MX(gamma_MX < 6) = -abs(epsilon_MX(gamma_MX < 6));
                    case {'hogervorst' 'hogervorst-wbk'}
                        gamma_MM = Param.gamma_MM; % Unitless
                        gamma_XX = Param.gamma_XX; % Unitless
                        gamma_MX = (gamma_MM + gamma_XX)./2;

                        epsilon_MX = 2.*epsilon_MM.*epsilon_XX./(epsilon_MM + epsilon_XX);

                        epsilon_MM(gamma_MM < 6) = -abs(epsilon_MM(gamma_MM < 6));
                        epsilon_XX(gamma_XX < 6) = -abs(epsilon_XX(gamma_XX < 6));
                        epsilon_MX(gamma_MX < 6) = -abs(epsilon_MX(gamma_MX < 6));

                        r0_MX = ( sqrt( ( epsilon_MM.*epsilon_XX.*gamma_MM.*gamma_XX.*(r0_MM.*r0_XX).^6 )...
                            ./((gamma_MM - 6).*(gamma_XX - 6)) ).*(gamma_MX - 6)./(epsilon_MX.*gamma_MX) ).^(1/6);

                    case {'kong' 'gromacs' 'kong-sr' 'gromacs-sr' 'kong-mid' 'gromacs-mid' 'kong-lr' 'gromacs-lr'}
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
                            case {'kong' 'kong-sr' 'kong-mid' 'kong-lr'}
                                A_MX = (1/2).*( A_MM.*(A_MM.*B_MM./(A_XX.*B_XX)).^(-B_MM./(B_MM + B_XX)) + ...
                                                A_XX.*(A_XX.*B_XX./(A_MM.*B_MM)).^(-B_XX./(B_MM + B_XX)) );
                                B_MX = 2.*B_MM.*B_XX./(B_MM + B_XX);
                                C_MX = sqrt(C_MM.*C_XX);
                            case {'gromacs' 'gromacs-sr' 'gromacs-mid' 'gromacs-lr'}
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
                gamma_MX = Param.gamma_MX; % Unitless
                epsilon_MM(gamma_MM < 6) = -abs(epsilon_MM(gamma_MM < 6));
                epsilon_XX(gamma_XX < 6) = -abs(epsilon_XX(gamma_XX < 6));
                epsilon_MX(gamma_MX < 6) = -abs(epsilon_MX(gamma_MX < 6));
            end

            % Remove any nonsense
            r0_MX(abs(imag(r0_MX)) > sqrt(eps)) = nan;
            gamma_MX(abs(imag(gamma_MX)) > sqrt(eps)) = nan;
            epsilon_MX(abs(imag(epsilon_MX)) > sqrt(eps)) = nan;
            r0_MX = real(r0_MX);
            gamma_MX = real(gamma_MX);
            epsilon_MX = real(epsilon_MX);

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
    case 'BF'
        % Input parameters
        sigma_MM = Param.sigma_MM; % nm
        sigma_XX = Param.sigma_XX; % nm

        epsilon_MM = Param.epsilon_MM; % kJ/mol
        epsilon_XX = Param.epsilon_XX; % kJ/mol

        if Settings.Additivity
                switch lower(Settings.Comb_rule)
                    case 'lorentz-berthelot'
                        gamma_MX = Param.gamma_MX; % Unitless
                        sigma_MX = (sigma_MM + sigma_XX)./2; % nm
                        epsilon_MX = sqrt(epsilon_MM.*epsilon_XX); % kJ/mol
                        gamma_MM = gamma_MX; % Unitless
                        gamma_XX = gamma_MX; % Unitless
                    case 'hogervorst'
                        gamma_MM = Param.gamma_MM; % Unitless
                        gamma_XX = Param.gamma_XX; % Unitless
                        gamma_MX = (gamma_MM + gamma_XX)./2;

                        epsilon_MX = 2.*epsilon_MM.*epsilon_XX./(epsilon_MM + epsilon_XX);

                        sigma_MX = ( sqrt( ( epsilon_MM.*epsilon_XX.*gamma_MM.*gamma_XX.*(sigma_MM.*sigma_XX).^6 )...
                            ./((gamma_MM - 6).*(gamma_XX - 6)) ).*(gamma_MX - 6)./(epsilon_MX.*gamma_MX) ).^(1/6);
                    case {'hogervorst-wbk'}
                        gamma_MM = Param.gamma_MM; % Unitless
                        gamma_XX = Param.gamma_XX; % Unitless
                        gamma_MX = (gamma_MM + gamma_XX)./2;
                        
                        epsilon_MX = 2.*epsilon_MM.*epsilon_XX./(epsilon_MM + epsilon_XX);
                        
                        sigma_MX = ( sqrt( ( epsilon_MM.*epsilon_XX.*(gamma_MM + 3).*(gamma_XX + 3).*(sigma_MM.*sigma_XX).^6 )...
                            ./(gamma_MM.*gamma_XX) ).*gamma_MX./(epsilon_MX.*(gamma_MX + 3)) ).^(1/6);
                    case {'kong' 'gromacs' 'kong-sr' 'gromacs-sr' 'kong-mid' 'gromacs-mid' 'kong-lr' 'gromacs-lr'}
                        gamma_MM = Param.gamma_MM; % Unitless
                        gamma_XX = Param.gamma_XX; % Unitless
                        
                        % These apply exactly in all cases
                        B_MM = gamma_MM./sigma_MM; % exponent
                        B_XX = gamma_XX./sigma_XX;
                        
                        switch lower(Settings.Comb_rule)
                            case {'kong' 'kong-sr' 'gromacs' 'gromacs-sr'}
                                % Define approximate tight form constants (valid as r->0)
                                A_MM = 6.*epsilon_MM.*exp(gamma_MM)./gamma_MM; % prefactor
                                A_XX = 6.*epsilon_XX.*exp(gamma_XX)./gamma_XX;
                                
                                C_MM = 2.*epsilon_MM.*(gamma_MM + 3)./gamma_MM; % dispersion
                                C_XX = 2.*epsilon_XX.*(gamma_XX + 3)./gamma_XX;
                                
                                % Apply mixing rules
                                switch lower(Settings.Comb_rule)
                                    case {'kong' 'kong-sr'}
                                        A_MX = (1/2).*( A_MM.*(A_MM.*B_MM./(A_XX.*B_XX)).^(-B_MM./(B_MM + B_XX)) + ...
                                                        A_XX.*(A_XX.*B_XX./(A_MM.*B_MM)).^(-B_XX./(B_MM + B_XX)) );
                                        B_MX = 2.*B_MM.*B_XX./(B_MM + B_XX);
                                        C_MX = sqrt(C_MM.*C_XX);
                                    case {'gromacs' 'gromacs-sr'}
                                        A_MX = sqrt(A_MM.*A_XX);
                                        B_MX = 2./( (1./B_MM) + (1./B_XX) );
                                        C_MX = sqrt(C_MM.*C_XX);
                                end
                                
                                % Convert MX parameters back to gamma/epsilon/sigma
                                gamma_MX = -lambertw(-1,-3*C_MX./(A_MX.*exp(3))) - 3;
                                sigma_MX = gamma_MX./B_MX;
                                epsilon_MX = A_MX.*gamma_MX.*exp(-gamma_MX)./6;
                                
                            case {'kong-mid' 'gromacs-mid'}
                                % Define approximate tight form constants (valid as r->0)
                                A_MM = 3.*epsilon_MM.*exp(gamma_MM)./gamma_MM; % prefactor
                                A_XX = 3.*epsilon_XX.*exp(gamma_XX)./gamma_XX;
                                
                                C_MM = epsilon_MM.*(gamma_MM + 3)./gamma_MM; % dispersion
                                C_XX = epsilon_XX.*(gamma_XX + 3)./gamma_XX;
                                
                                % Apply mixing rules
                                switch lower(Settings.Comb_rule)
                                    case 'kong-mid'
                                        A_MX = (1/2).*( A_MM.*(A_MM.*B_MM./(A_XX.*B_XX)).^(-B_MM./(B_MM + B_XX)) + ...
                                                        A_XX.*(A_XX.*B_XX./(A_MM.*B_MM)).^(-B_XX./(B_MM + B_XX)) );
                                        B_MX = 2.*B_MM.*B_XX./(B_MM + B_XX);
                                        C_MX = sqrt(C_MM.*C_XX);
                                    case 'gromacs-mid'
                                        A_MX = sqrt(A_MM.*A_XX);
                                        B_MX = 2./( (1./B_MM) + (1./B_XX) );
                                        C_MX = sqrt(C_MM.*C_XX);
                                end
                                
                                % Convert MX parameters back to gamma/epsilon/sigma
                                gamma_MX = -lambertw(-1,-3*C_MX./(A_MX.*exp(3))) - 3;
                                sigma_MX = gamma_MX./B_MX;
                                epsilon_MX = A_MX.*gamma_MX.*exp(-gamma_MX)./3;
                                
                            case {'kong-lr' 'gromacs-lr'}
                                % Define approximate tight form constants (valid as r->0)
                                A_MM = 6.*epsilon_MM.*exp(gamma_MM).*(sigma_MM.^6)./gamma_MM; % prefactor
                                A_XX = 6.*epsilon_XX.*exp(gamma_XX).*(sigma_XX.^6)./gamma_XX;
                                
                                C_MM = 2.*epsilon_MM.*(gamma_MM + 3).*(sigma_MM.^6)./gamma_MM; % dispersion
                                C_XX = 2.*epsilon_XX.*(gamma_XX + 3).*(sigma_XX.^6)./gamma_XX;
                                
                                % Apply mixing rules
                                switch lower(Settings.Comb_rule)
                                    case 'kong-lr'
                                        A_MX = (1/2).*( A_MM.*(A_MM.*B_MM./(A_XX.*B_XX)).^(-B_MM./(B_MM + B_XX)) + ...
                                                        A_XX.*(A_XX.*B_XX./(A_MM.*B_MM)).^(-B_XX./(B_MM + B_XX)) );
                                        B_MX = 2.*B_MM.*B_XX./(B_MM + B_XX);
                                        C_MX = sqrt(C_MM.*C_XX);
                                    case 'gromacs-lr'
                                        A_MX = sqrt(A_MM.*A_XX);
                                        B_MX = 2./( (1./B_MM) + (1./B_XX) );
                                        C_MX = sqrt(C_MM.*C_XX);
                                end
                                
                                % Convert MX parameters back to gamma/epsilon/sigma
                                gamma_MX = -lambertw(-1,-3*C_MX./(A_MX.*exp(3))) - 3;
                                sigma_MX = gamma_MX./B_MX;
                                epsilon_MX = A_MX.*(B_MX.^6).*(gamma_MX.^(-5)).*exp(-gamma_MX)./6;
                        end
                end
        else
            sigma_MX = Param.sigma_MX; % nm
            epsilon_MX = Param.epsilon_MX; % kJ/mol
            gamma_MM = Param.gamma_MM; % Unitless
            gamma_XX = Param.gamma_XX; % Unitless
            gamma_MX = Param.gamma_MX; % Unitless
        end

        % Remove any complex numbers
        sigma_MX(abs(imag(sigma_MX)) > sqrt(eps)) = nan;
        gamma_MX(abs(imag(gamma_MX)) > sqrt(eps)) = nan;
        epsilon_MX(abs(imag(epsilon_MX)) > sqrt(eps)) = nan;
        sigma_MX = real(sigma_MX);
        gamma_MX = real(gamma_MX);
        epsilon_MX = real(epsilon_MX);

        % Outputs
        Settings.S.S.MM = sigma_MM;
        Settings.S.S.XX = sigma_XX;
        Settings.S.S.MX = sigma_MX;

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
    case 'JC'
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

            % Sigma
            Settings.S.S.MM = Param.sigma_MM;
            Settings.S.S.XX = Param.sigma_XX;

            % Epsilon
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

        % Dispersion/repulsion form
        else
            % Dispersion
            Settings.S.D.MM = Param.SDMM;
            Settings.S.D.XX = Param.SDXX;

            % Repulsion
            Settings.S.R.MM = Param.SRMM;
            Settings.S.R.XX = Param.SRXX;

            if Settings.Additivity

                % Scaled
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

                Settings.S.S.MX = (Param.sigma_MM + Param.sigma_XX)./2;
                Settings.S.E.MX = sqrt(Param.epsilon_MM.*Param.epsilon_XX);

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

    % Calculate loss due to infeasible models with no well minima (only works reliably for BH/TF models in sigma-epsilon form)
    if strcmp(Settings.Theory,'BH')
        U = BH_Potential_Generator_vec(Settings);
    elseif strcmp(Settings.Theory,'BD')
        U = BD_Potential_Generator_vec(Settings);
    elseif strcmp(Settings.Theory,'BE')
        U = BE_Potential_Generator_vec(Settings);
    elseif strcmp(Settings.Theory,'BF')
        U = BF_Potential_Generator_vec(Settings);
    elseif strcmp(Settings.Theory,'TF')
        U = TF_Potential_Generator_vec(Settings);
    elseif strcmp(Settings.Theory,'JC')
        U = JC_Potential_Generator_vec(Settings);
    elseif strcmp(Settings.Theory,'LJ')
        U = LJ_Potential_Generator_vec(Settings);
    elseif strcmp(Settings.Theory,'Mie')
        U = Mie_Potential_Generator_vec(Settings);
    end

    %% Check MX attractive potential
    % Check for NaN functions
    badfcn_idx = any(isnan(U.MX),2);
    Loss(badfcn_idx) = Loss(badfcn_idx) + Settings.BadFcnLossPenalty;

    % Grab the peaks and valleys of the MX attractive potential
    peaks_idx = islocalmax(U.MX,2,'MinProminence',1e-8);
    valleys_idx = islocalmin(U.MX,2,'MinProminence',1e-8);
    r = repmat(U.r,size(U.MX,1),1);

    Num_peaks = sum(peaks_idx,2);
    Num_valleys = sum(valleys_idx,2);

    nv_idx = Num_valleys == 0;                        % Potentials that contain no valley
    ovnp_idx = (Num_valleys == 1) & (Num_peaks == 0); % Potentials with 1 valley and no peak
    ovop_idx = (Num_valleys == 1) & (Num_peaks == 1); % Potentials with 1 peak and 1 valley

    % If no well minimum exists in MX interaction, add a loss penalty
    Loss(nv_idx) = Loss(nv_idx) + Settings.BadFcnLossPenalty;

    % If no peak exists in MX interaction, but valley does, check well depth. 
    % This is normal for JC potential and some BH/TF potentials
    U_MX_ovnp_set = U.MX(ovnp_idx,:);
    r_ovnp_set = r(ovnp_idx,:);
    ov_idx = sum( cumprod(valleys_idx(ovnp_idx,:) == 0, 2), 2) + 1;
    ovl_idx = sub2ind(size(U_MX_ovnp_set),(1:numel(ov_idx)).',ov_idx);
    U_min = U_MX_ovnp_set(ovl_idx);
    Loss(ovnp_idx) = Loss(ovnp_idx) + max(Settings.MaxAttWellDepth - U_min,0).*Settings.BadFcnLossPenalty; % valley too deep
    % Also check well location
    r_min = r_ovnp_set(ovl_idx);
    Loss(ovnp_idx) = Loss(ovnp_idx) + max(r_min - (Settings.MaxMXWellR/10),0).*Settings.BadFcnLossPenalty; % valley too far out
    Loss(ovnp_idx) = Loss(ovnp_idx) + max((Settings.MinMXWellR/10) - r_min,0).*Settings.BadFcnLossPenalty; % valley too close in

    % If both a peak and minimum exist, find the well depth
    U_MX_ovop_set = U.MX(ovop_idx,:);
    r_ovop_set = r(ovop_idx,:);
    op_idx = sum( cumprod(peaks_idx(ovop_idx,:) == 0, 2), 2)  + 1;
    opl_idx = sub2ind(size(U_MX_ovop_set),(1:numel(op_idx)).',op_idx);
    ov_idx = sum( cumprod(valleys_idx(ovop_idx,:) == 0, 2), 2) + 1;
    ovl_idx = sub2ind(size(U_MX_ovop_set),(1:numel(ov_idx)).',ov_idx);
    dU = U_MX_ovop_set(opl_idx) - U_MX_ovop_set(ovl_idx);
    U_min = U_MX_ovop_set(ovl_idx);
    % Penalize walls that are too small and valleys that are too deep
    Loss(ovop_idx) = Loss(ovop_idx) + max(Settings.MinExpWallHeight - dU,0).*Settings.BadFcnLossPenalty; % repulsive wall too low
    Loss(ovop_idx) = Loss(ovop_idx) + max(Settings.MaxAttWellDepth - U_min,0).*Settings.BadFcnLossPenalty; % valley too deep
    % Also check valley location
    r_min = r_ovop_set(ovl_idx);
    Loss(ovop_idx) = Loss(ovop_idx) + max(r_min - (Settings.MaxMXWellR/10),0).*Settings.BadFcnLossPenalty; % valley too far
    Loss(ovop_idx) = Loss(ovop_idx) + max((Settings.MinMXWellR/10) - r_min,0).*Settings.BadFcnLossPenalty; % valley too close

%     % Loss = zeros(N_par,1);
%     idxes = find(Loss ~= 0);
%     for jdx = 1:min(numel(idxes),100)
%         idx = idxes(jdx);
%         plot(U.r,U.MX(idx,:))
%         hold on
%         scatter(U.r(peaks_idx(idx,:)),U.MX(idx,peaks_idx(idx,:)))
%         scatter(U.r(valleys_idx(idx,:)),U.MX(idx,valleys_idx(idx,:)))
%     end
%     ylim([-1000 1000])
%     xlim([0 1])

    %% Grab the peaks and valleys of the MM/XX potentials
    YY = {'MM' 'XX'};
    r_wall_sv = cell(2,1);
    for j = 1:2
        jj = YY{j};

        peaks_idx = islocalmax(U.(jj),2,'MinProminence',1e-8);
        valleys_idx = islocalmin(U.(jj),2,'MinProminence',1e-8);

        Num_peaks = sum(peaks_idx,2);
        Num_valleys = sum(valleys_idx,2);

        % 5 possible cases
        %npnv_idx = (Num_peaks == 0) & (Num_valleys == 0); % Potentials that contain no peaks and no valleys
        opnv_idx = (Num_peaks == 1) & (Num_valleys == 0); % Potentials that contain one peak and no valleys
        tpov_idx = (Num_peaks == 2);                      % Potentials that countain 2 peaks and one valley
        opov_idx = (Num_peaks == 1) & (Num_valleys == 1); % Potentials that contain one peak and one valley
        npov_idx = (Num_peaks == 0) & (Num_valleys == 1); % Potentials that contain no peaks and one valley

        % No peak and no valley exists: Do nothing, this is normal for JC and sometimes BH/TF
        % Since peak height is presumably infinite, do not need to check it.
        %npnv_idx;

        % One peak, no valley: This is a normal case for TF and BH
        % repulsive potentials, check to ensure the peak height is at
        % least Settings.MinExpWallHeight
        opov_set = peaks_idx(opnv_idx,:);
        U_jj_opnv_set = U.(jj)(opnv_idx,:);
        op_idx = sum(cumprod(opov_set == 0, 2), 2) + 1;
        opl_idx = sub2ind(size(U_jj_opnv_set),(1:numel(op_idx)).',op_idx);
        U_peak = U_jj_opnv_set(opl_idx);
        Loss(opnv_idx) = Loss(opnv_idx) + max(Settings.MinExpWallHeight - U_peak,0).*Settings.BadFcnLossPenalty; % peak too low

        % two peaks are visible implies one valley in between them: 
        % Penalize any model with a non-zero well depth greater than Settings.MaxRepWellDepth
        U_jj_tpov_set = U.(jj)(tpov_idx,:);
        sp_idx = size(U.r,2) - sum( cumprod(fliplr(peaks_idx(tpov_idx,:)) == 0, 2), 2); % second peak
        spl_idx = sub2ind(size(U_jj_tpov_set),(1:numel(sp_idx)).',sp_idx);
        fp_idx = sum( cumprod(peaks_idx(tpov_idx,:) == 0, 2), 2) + 1; % first peak
        fpl_idx = sub2ind(size(U_jj_tpov_set),(1:numel(fp_idx)).',fp_idx);

        ov_idx = sum( cumprod(valleys_idx(tpov_idx,:) == 0, 2), 2) + 1;
        ovl_idx = sub2ind(size(U_jj_tpov_set),(1:numel(ov_idx)).',ov_idx);
        dU_sp = U_jj_tpov_set(spl_idx) - U_jj_tpov_set(ovl_idx);
        Loss(tpov_idx) = Loss(tpov_idx) + max(dU_sp - Settings.MaxRepWellDepth,0).*Settings.BadFcnLossPenalty; % valley exists
        % Also check the peak height
        dU_fp = U_jj_tpov_set(fpl_idx) - U_jj_tpov_set(ovl_idx);
        Loss(tpov_idx) = Loss(tpov_idx) + max(Settings.MinExpWallHeight - dU_fp,0).*Settings.BadFcnLossPenalty; % peak too low

        % One peak visible + one valley, and (possibly) one hidden peak to the right
        % In this case, check if the valley is closer-in than the peak
        % If the valley is closer-in, there is no hidden peak, use dU = peak_U - valley_U
        % If the peak is closer-in, there is a second hidden peak to the right, use dU = U(end) - valley_U
        opov_set = peaks_idx(opov_idx,:);
        U_jj_opov_set = U.(jj)(opov_idx,:);
        op_idx = sum(cumprod(opov_set == 0, 2), 2) + 1;
        ov_idx = sum(cumprod(valleys_idx(opov_idx,:) == 0, 2), 2) + 1;

        % split into two subsets
        ovcc_idx = (ov_idx < op_idx); % valley is closer in
        Lossn = zeros(size(ovcc_idx));

        % Valley closer set: in this subset, there is no inner peak. Hence no repulsive wall peak to worry about 
        % use peak_U - valley_U for the valley depth
        U_jj_opovc_set = U_jj_opov_set(ovcc_idx,:);
        opc_idx = op_idx(ovcc_idx,:);
        opcl_idx = sub2ind(size(U_jj_opovc_set),(1:numel(opc_idx)).',opc_idx);
        ovc_idx = ov_idx(ovcc_idx,:);
        ovcl_idx = sub2ind(size(U_jj_opovc_set),(1:numel(ovc_idx)).',ovc_idx);
        dUc = U_jj_opovc_set(opcl_idx) - U_jj_opovc_set(ovcl_idx);
        Lossn(ovcc_idx) =  max(dUc - Settings.MaxRepWellDepth,0).*Settings.BadFcnLossPenalty; % valley exists

        % Valley further set: In this subset, there is an inner peak but no outer peak.
        % U(end) - valley_U for valley depth
        U_jj_opovf_set = U_jj_opov_set(~ovcc_idx,:);
        fpf_idx = op_idx(~ovcc_idx,:); % set the first peak
        fpfl_idx = sub2ind(size(U_jj_opovf_set),(1:numel(fpf_idx)).',fpf_idx);
        spf_idx = repmat(length(U.r),sum(~ovcc_idx),1); % set the second "peak" to the furthest possible place
        spfl_idx = sub2ind(size(U_jj_opovf_set),(1:numel(spf_idx)).',spf_idx);
        ovf_idx = ov_idx(~ovcc_idx,:);
        ovfl_idx = sub2ind(size(U_jj_opovf_set),(1:numel(ovf_idx)).',ovf_idx);

        U_ov = U_jj_opovf_set(ovfl_idx); % valley energy
        U_fp = U_jj_opovf_set(fpfl_idx); % first peak energy
        U_sp = U_jj_opovf_set(spfl_idx); % second "peak" energy
        dUf = U_sp - U_ov; % Well depth relative to the furthest out place
        dU_fp = U_fp - U_ov; % 
        Lossn(~ovcc_idx) =  max(dUf - Settings.MaxRepWellDepth,0).*Settings.BadFcnLossPenalty; % valley exists
        % Also check the repulsive peak height
        Lossn(~ovcc_idx) = Lossn(~ovcc_idx) + max(Settings.MinExpWallHeight - dU_fp,0).*Settings.BadFcnLossPenalty; % wall too low
        % Add everything back
        Loss(opov_idx) = Loss(opov_idx) + Lossn;

        % No peak + one valley: there must be a hidden peak to the right
        % use dU = U(end) - valley_U
        U_jj_npov_set = U.(jj)(npov_idx,:);
        ov_idx = sum(cumprod(valleys_idx(npov_idx,:) == 0, 2), 2) + 1;
        ovl_idx = sub2ind(size(U_jj_npov_set),(1:numel(ov_idx)).',ov_idx);
        op_idx = repmat(length(U.r),sum(npov_idx),1); % set the "peak" to the furthest possible place
        opl_idx = sub2ind(size(U_jj_npov_set),(1:numel(op_idx)).',op_idx);
        dU = U_jj_npov_set(opl_idx) - U_jj_npov_set(ovl_idx); % Peak height
        Loss(npov_idx) = Loss(npov_idx) +  max(dU - Settings.MaxRepWellDepth,0).*Settings.BadFcnLossPenalty; % valley exists

        % For all functions, ensure the repulsive wall is not too far out
        walls_idx = any(U.(jj) >= Settings.MinExpWallHeight,2);
        U_jj_rw_set = U.(jj)(walls_idx,:);
        r_jj_rw_set = r(walls_idx,:);
        walls_set = U_jj_rw_set >= Settings.MinExpWallHeight;
        rw_idx = size(U_jj_rw_set,2) - sum( cumprod(fliplr(walls_set(:,:)) == 0, 2), 2); % where the wall reaches min height
        rwl_idx = sub2ind(size(U_jj_rw_set),(1:numel(rw_idx)).',rw_idx);
        r_wall = r_jj_rw_set(rwl_idx);
        Loss(walls_idx) = Loss(walls_idx) + max(r_wall - (Settings.MaxMXWellR/10),0).*Settings.BadFcnLossPenalty; % wall too far
        Loss(walls_idx) = Loss(walls_idx) + max((Settings.MinMXWellR/10) - r_wall,0).*Settings.BadFcnLossPenalty; % wall too close

        r_wall_sv{j} = nan(size(U.MX,1),1);
        r_wall_sv{j}(walls_idx) = r_wall;

    %     % Loss = zeros(N_par,1);
    %     idxes = find(Loss == 0);
    %     for jdx = 1:min(numel(idxes),100)
    %         idx=idxes(jdx);
    %         hold on
    %         plot(U.r,U.(jj)(idx,:))
    %         scatter(U.r(peaks_idx(idx,:)),U.(jj)(idx,peaks_idx(idx,:)))
    %         scatter(U.r(valleys_idx(idx,:)),U.(jj)(idx,valleys_idx(idx,:)))
    %     end
    %     ylim([-1000 1000])
    %     xlim([0 2])
    end

    if Settings.EnforceRR && Settings.SigmaEpsilon && Settings.Additivity
        switch Settings.Theory
            case {'JC' 'Mie' 'BF' 'LJ'}
                RR_Walls = Param.sigma_MM./Param.sigma_XX; % Ratio of M/X size, should be < 1
            case {'BH' 'BD' 'BE' 'TF'}
                RR_Walls = Param.r0_MM./Param.r0_XX; % Ratio of M/X size, should be < 1
        end
        Loss = Loss + max(RR_Walls - 1,0).*Settings.BadFcnLossPenalty; % Radius ratios are incorrect
    elseif Settings.EnforceRR
        % Check that the repulsive wall of the XX interaction is further out than the MM interaction
        RR_Walls = r_wall_sv{1}./r_wall_sv{2}; % Ratio of M/X size, should be < 1
        Loss = Loss + max(RR_Walls - 1,0).*Settings.BadFcnLossPenalty; % Radius ratios are incorrect
    end
end

if Output_Loss
    tf = Loss;
else
    tf = log1p(Loss) < sqrt(eps);
end

end