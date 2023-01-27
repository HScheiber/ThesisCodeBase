% This function assumes sigma-epsilon form
function Output = Gen_Param_Table(Settings,Param,varargin)

DefParams = bayesopt_params(Settings);
ParNames = {DefParams.Name};
Output = struct;
for idx = 1:numel(Param)
    Output.(ParNames{idx}) = Param(idx);
end

if nargin < 3
    Negative_epsilon = true;
else
    Negative_epsilon = varargin{1};
end


switch Settings.Theory
case 'TF'
    % Loose form of exp-C6-C8 model
    if Settings.SigmaEpsilon
        if Settings.Additivity
            Output.r0_MX = (Output.r0_MM + Output.r0_XX)/2; % nm
            Output.epsilon_MX = sqrt(Output.epsilon_MM*Output.epsilon_XX); % kJ/mol
            Output.gamma_MM = Output.gamma_MX; % Unitless
            Output.gamma_XX = Output.gamma_MX; % Unitless
        end

        if Output.gamma_MX > 48/7
            Output.epsilon_MX = -Output.epsilon_MX;
        end
        if Output.gamma_MM > 48/7
            Output.epsilon_MM = -Output.epsilon_MM;
        end
        if Output.gamma_XX > 48/7
            Output.epsilon_XX = -Output.epsilon_XX;
        end

        % Convert to Condensed form
        Output.alpha_MM = Output.gamma_MM/Output.r0_MM; % nm^(-1)
        Output.alpha_XX = Output.gamma_XX/Output.r0_XX; % nm^(-1)
        Output.alpha_MX = Output.gamma_MX/Output.r0_MX; % nm^(-1)

        Output.B_MM = 48*Output.epsilon_MM*exp(Output.gamma_MM)/(48 - 7*Output.gamma_MM); % kJ/mol
        Output.B_XX = 48*Output.epsilon_XX*exp(Output.gamma_XX)/(48 - 7*Output.gamma_XX); % kJ/mol
        Output.B_MX = 48*Output.epsilon_MX*exp(Output.gamma_MX)/(48 - 7*Output.gamma_MX); % kJ/mol

        Output.C_MM = 4*Output.epsilon_MM*Output.gamma_MM*(Output.r0_MM^6)/(48 - 7*Output.gamma_MM); % kJ/mol nm^6
        Output.C_XX = 4*Output.epsilon_XX*Output.gamma_XX*(Output.r0_XX^6)/(48 - 7*Output.gamma_XX); % kJ/mol nm^6
        Output.C_MX = 4*Output.epsilon_MX*Output.gamma_MX*(Output.r0_MX^6)/(48 - 7*Output.gamma_MX); % kJ/mol nm^6

        Output.D_MM = 3*Output.epsilon_MM*Output.gamma_MM*(Output.r0_MM^8)/(48 - 7*Output.gamma_MM); % kJ/mol nm^8
        Output.D_XX = 3*Output.epsilon_XX*Output.gamma_XX*(Output.r0_XX^8)/(48 - 7*Output.gamma_XX); % kJ/mol nm^8
        Output.D_MX = 3*Output.epsilon_MX*Output.gamma_MX*(Output.r0_MX^8)/(48 - 7*Output.gamma_MX); % kJ/mol nm^8

        % Scaling Coulombic Charge
        if Settings.Fix_Charge
            Output.Q = Settings.Q_value;
        end
        
    else
        
        % 1/R8 Dispersion (TF only)
        if Settings.Fix_C8
            % Calculate value of C8 using recursive relations
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
            
            % Conversion factors
            Bohr_nm = 0.0529177; % a_0 - > Angstrom
            c6conv = 1e-3/2625.4999/((0.052917726)^6); % J/mol nm^6 - > au (from sourcecode)
            J_kJ = 1e-3; % J - > kJ
            Ha_kJmol = 2625.4999; % Ha - > kJ/mol
            c6units = (1/c6conv)*J_kJ; % au - > kJ/mol nm^6
            c8units = (Ha_kJmol)*(Bohr_nm^8); % au - > kJ/mol nm^8

            % Calculate Scaled C8 using recursion relation from D3 paper
            Output.SD8MM = 3.0.*(Output.SD6MM./c6units).*sqrt_Q.(Settings.Metal).*sqrt_Q.(Settings.Metal).*c8units; % in kJ/mol nm^8
            Output.SD8XX = 3.0.*(Output.SD6XX./c6units).*sqrt_Q.(Settings.Halide).*sqrt_Q.(Settings.Halide).*c8units; % in kJ/mol nm^8
            Output.SD8MX = 3.0.*(Output.SD6MX./c6units).*sqrt_Q.(Settings.Metal).*sqrt_Q.(Settings.Halide).*c8units; % in kJ/mol nm^8
        end

        % Alpha (TF exponential steepness repulsive parameter)
        if Settings.Fix_Alpha
            Output.SAMM = Output.SAMX;
            Output.SAXX = Output.SAMX;
        end

        % Scaling Coulombic Charge
        if Settings.Fix_Charge
            Output.Q = Settings.Q_value;
        end
    end
    
case {'BH' 'BD' 'BE'} % Buckingham model
    if Settings.SigmaEpsilon
        if Settings.Additivity
        	switch lower(Settings.Comb_rule)
                case 'lorentz-berthelot'
                    Output.r0_MX      = (Output.r0_MM + Output.r0_XX)/2; % nm
                    Output.epsilon_MX = sqrt(Output.epsilon_MM*Output.epsilon_XX); % kJ/mol
                    Output.gamma_MM   = Output.gamma_MX; % Unitless
                    Output.gamma_XX   = Output.gamma_MX; % Unitless
                    
                    if Negative_epsilon
                        if Output.gamma_MX < 6
                            Output.epsilon_MX = -Output.epsilon_MX;
                        end
                        if Output.gamma_MM < 6
                            Output.epsilon_MM = -Output.epsilon_MM;
                        end
                        if Output.gamma_XX < 6
                            Output.epsilon_XX = -Output.epsilon_XX;
                        end
                    end
                    
                case {'hogervorst' 'hogervorst-wbk'}
                    Output.gamma_MX = (Output.gamma_MM + Output.gamma_XX)/2;
                    
                    Output.epsilon_MX = 2*Output.epsilon_MM.*Output.epsilon_XX/(Output.epsilon_MM + Output.epsilon_XX);
                    
                    if Output.gamma_MX < 6
                        Output.epsilon_MX = -Output.epsilon_MX;
                    end
                    if Output.gamma_MM < 6
                        Output.epsilon_MM = -Output.epsilon_MM;
                    end
                    if Output.gamma_XX < 6
                        Output.epsilon_XX = -Output.epsilon_XX;
                    end
                    
                    Output.r0_MX = ( sqrt( ( Output.epsilon_MM*Output.epsilon_XX*Output.gamma_MM*Output.gamma_XX*(Output.r0_MM*Output.r0_XX)^6 )...
                        /((Output.gamma_MM - 6)*(Output.gamma_XX - 6)) )*(Output.gamma_MX - 6)/(Output.epsilon_MX*Output.gamma_MX) )^(1/6);
                    
                    if ~Negative_epsilon
                        if Output.gamma_MX < 6
                            Output.epsilon_MX = -Output.epsilon_MX;
                        end
                        if Output.gamma_MM < 6
                            Output.epsilon_MM = -Output.epsilon_MM;
                        end
                        if Output.gamma_XX < 6
                            Output.epsilon_XX = -Output.epsilon_XX;
                        end
                    end
                    
                case {'kong' 'gromacs'}
                    
                    if Output.gamma_MM < 6
                        Output.epsilon_MM = -Output.epsilon_MM;
                    end
                    if Output.gamma_XX < 6
                        Output.epsilon_XX = -Output.epsilon_XX;
                    end
                    
                    k_MM = 1 / ( Output.gamma_MM - 6 );
                    k_XX = 1 / ( Output.gamma_XX - 6 );
                    
                    A_MM = 6*Output.epsilon_MM*k_MM*exp(Output.gamma_MM); % prefactor
                    A_XX = 6*Output.epsilon_XX*k_XX*exp(Output.gamma_XX);
                    
                    B_MM = Output.gamma_MM/Output.r0_MM; % exponent
                    B_XX = Output.gamma_XX/Output.r0_XX;
                    
                    C_MM = Output.epsilon_MM*Output.gamma_MM*k_MM*(Output.r0_MM^6); % dispersion
                    C_XX = Output.epsilon_XX*Output.gamma_XX*k_XX*(Output.r0_XX^6);
                    
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
                    Output.gamma_MX   = -7*lambertw((-1/7)*(6*C_MX*(B_MX^6)/A_MX)^(1/7));
                    Output.r0_MX      = Output.gamma_MX/B_MX;
                    Output.epsilon_MX = C_MX*(Output.gamma_MX - 6)/(Output.gamma_MX*(Output.r0_MX^6));
                    if Output.gamma_MX < 6
                        Output.epsilon_MX = -abs(Output.epsilon_MX);
                    end
                    if ~Negative_epsilon
                        if Output.gamma_MX < 6
                            Output.epsilon_MX = -Output.epsilon_MX;
                        end
                        if Output.gamma_MM < 6
                            Output.epsilon_MM = -Output.epsilon_MM;
                        end
                        if Output.gamma_XX < 6
                            Output.epsilon_XX = -Output.epsilon_XX;
                        end
                    end
        	end
        end
        
        % Convert to Condensed form
        Output.alpha_MM = Output.gamma_MM/Output.r0_MM; % nm^(-1)
        Output.alpha_XX = Output.gamma_XX/Output.r0_XX; % nm^(-1)
        Output.alpha_MX = Output.gamma_MX/Output.r0_MX; % nm^(-1)

        Output.B_MM = 6*Output.epsilon_MM*exp(Output.gamma_MM)/(Output.gamma_MM - 6); % kJ/mol
        Output.B_XX = 6*Output.epsilon_XX*exp(Output.gamma_XX)/(Output.gamma_XX - 6); % kJ/mol
        Output.B_MX = 6*Output.epsilon_MX*exp(Output.gamma_MX)/(Output.gamma_MX - 6); % kJ/mol

        Output.C_MM = Output.epsilon_MM*Output.gamma_MM*(Output.r0_MM^6)/(Output.gamma_MM - 6); % kJ/mol nm^6
        Output.C_XX = Output.epsilon_XX*Output.gamma_XX*(Output.r0_XX^6)/(Output.gamma_XX - 6); % kJ/mol nm^6
        Output.C_MX = Output.epsilon_MX*Output.gamma_MX*(Output.r0_MX^6)/(Output.gamma_MX - 6); % kJ/mol nm^6

        % Scaling Coulombic Charge
        if Settings.Fix_Charge
            Output.Q = Settings.Q_value;
        end
    else % Tight form of the exp-C6 model

        if Settings.Additivity
            Output.C_MX = sqrt(Output.C_MM.*Output.C_XX);
            Output.B_MX = sqrt(Output.B_MM.*Output.B_XX);

            if Settings.Additional_MM_Disp
                Output.SDMM = Output.SDMM + Output.SDMM2;
            end
        end
        
        % Alpha (exponential steepness repulsive parameter)
        if Settings.Fix_Alpha
            Output.SAMM = Output.SAMX;
            Output.SAXX = Output.SAMX;
        end

        % Scaling Coulombic Charge
        if Settings.Fix_Charge
            Output.Q = Settings.Q_value;
        end
    end
case 'BF' % Wang-Buckingham model
    if Settings.SigmaEpsilon
        if Settings.Additivity
            switch lower(Settings.Comb_rule)
                case 'lorentz-berthelot'
                    Output.sigma_MX = (Output.sigma_MM + Output.sigma_XX)/2; % nm
                    Output.epsilon_MX = sqrt(Output.epsilon_MM*Output.epsilon_XX); % kJ/mol
                    Output.gamma_MM = Output.gamma_MX; % Unitless
                    Output.gamma_XX = Output.gamma_MX; % Unitless
                case 'hogervorst'
                    Output.gamma_MX = (Output.gamma_MM + Output.gamma_XX)/2;
                    Output.epsilon_MX = 2*Output.epsilon_MM*Output.epsilon_XX/(Output.epsilon_MM + Output.epsilon_XX);
                    
                    Output.sigma_MX = ( sqrt( ( Output.epsilon_MM*Output.epsilon_XX*Output.gamma_MM*Output.gamma_XX*(Output.sigma_MM*Output.sigma_XX)^6 )...
                        /((Output.gamma_MM - 6)*(Output.gamma_XX - 6)) )*(Output.gamma_MX - 6)/(Output.epsilon_MX.*Output.gamma_MX) )^(1/6);
                case 'hogervorst-wbk'
                    Output.gamma_MX = (Output.gamma_MM + Output.gamma_XX)/2;
                    Output.epsilon_MX = 2*Output.epsilon_MM*Output.epsilon_XX/(Output.epsilon_MM + Output.epsilon_XX);

                    Output.sigma_MX = ( sqrt( ( Output.epsilon_MM.*Output.epsilon_XX.*(Output.gamma_MM + 3).*(Output.gamma_XX + 3).*(Output.sigma_MM.*Output.sigma_XX).^6 )...
                        ./(Output.gamma_MM.*Output.gamma_XX) ).*Output.gamma_MX./(Output.epsilon_MX.*(Output.gamma_MX + 3)) ).^(1/6);
                case {'kong' 'gromacs'}
                    
                    % Define approximate tight form constants (valid as r->0)
                    B_MM = 6.*Output.epsilon_MM.*exp(Output.gamma_MM)/Output.gamma_MM; % prefactor
                    B_XX = 6.*Output.epsilon_XX.*exp(Output.gamma_XX)/Output.gamma_XX;

                    alpha_MM = Output.gamma_MM./Output.sigma_MM; % exponent
                    alpha_XX = Output.gamma_XX./Output.sigma_XX;

                    C_MM = 2.*Output.epsilon_MM.*(Output.gamma_MM + 3)./Output.gamma_MM; % dispersion
                    C_XX = 2.*Output.epsilon_XX.*(Output.gamma_XX + 3)./Output.gamma_XX;

                    switch lower(Settings.Comb_rule)
                        case 'kong'
                            B_MX = (1/2).*( B_MM.*(B_MM.*alpha_MM./(B_XX.*alpha_XX)).^(-alpha_MM./(alpha_MM + alpha_XX)) + ...
                                            B_XX.*(B_XX.*alpha_XX./(B_MM.*alpha_MM)).^(-alpha_XX./(alpha_MM + alpha_XX)) );
                            alpha_MX = 2.*alpha_MM.*alpha_XX./(alpha_MM + alpha_XX);
                            C_MX = sqrt(C_MM.*C_XX);
                        case 'gromacs'
                            B_MX = sqrt(B_MM*B_XX);
                            alpha_MX = 2/( (1/alpha_MM) + (1/alpha_XX) );
                            C_MX = sqrt(C_MM*C_XX);
                    end

                    % Convert back to gamma/epsilon/r0
                    Output.gamma_MX = -lambertw(-1,-3*C_MX/(B_MX*exp(3))) - 3;
                    Output.sigma_MX = Output.gamma_MX/alpha_MX;
                    Output.epsilon_MX = B_MX*Output.gamma_MX*exp(-Output.gamma_MX)/6; % = gamma_MX*C_MX/(2*(gamma_MX + 3));
            end
        end

        % Convert to Condensed form
        Output.alpha_MM = Output.gamma_MM/Output.sigma_MM; % nm^(-1)
        Output.alpha_XX = Output.gamma_XX/Output.sigma_XX; % nm^(-1)
        Output.alpha_MX = Output.gamma_MX/Output.sigma_MX; % nm^(-1)

        Output.B_MM = 6*Output.epsilon_MM*exp(Output.gamma_MM)/(Output.gamma_MM - 6); % kJ/mol
        Output.B_XX = 6*Output.epsilon_XX*exp(Output.gamma_XX)/(Output.gamma_XX - 6); % kJ/mol
        Output.B_MX = 6*Output.epsilon_MX*exp(Output.gamma_MX)/(Output.gamma_MX - 6); % kJ/mol

        Output.C_MM = Output.epsilon_MM*Output.gamma_MM*(Output.sigma_MM^6)/(Output.gamma_MM - 6); % kJ/mol nm^6
        Output.C_XX = Output.epsilon_XX*Output.gamma_XX*(Output.sigma_XX^6)/(Output.gamma_XX - 6); % kJ/mol nm^6
        Output.C_MX = Output.epsilon_MX*Output.gamma_MX*(Output.sigma_MX^6)/(Output.gamma_MX - 6); % kJ/mol nm^6

        % Scaling Coulombic Charge
        if Settings.Fix_Charge
            Output.Q = Settings.Q_value;
        end
    else
        error('SigmaEpsilon form required for use with Wang-Buckingham potential form.')
    end
case 'JC'
    if Settings.SigmaEpsilon
        % D <-> sigma
        % R <-> epsilon
        if Settings.Additivity
            Output.sigma_MX = (Output.sigma_MM + Output.sigma_XX)/2;
            Output.epsilon_MX = sqrt(Output.epsilon_MM*Output.epsilon_XX); % kJ/mol
        end
        if Settings.Fix_Charge
            Output.Q = Settings.Q_value;
        end

        % Convert to Condensed form
        Output.B_MM = 4*Output.epsilon_MM*(Output.sigma_MM^12); % kJ/mol nm^(12)
        Output.B_XX = 4*Output.epsilon_XX*(Output.sigma_XX^12); % kJ/mol nm^(12)
        Output.B_MX = 4*Output.epsilon_MX*(Output.sigma_MX^12); % kJ/mol nm^(12)
        Output.C_MM = 4*Output.epsilon_MM*(Output.sigma_MM^6); % kJ/mol nm^(6)
        Output.C_XX = 4*Output.epsilon_XX*(Output.sigma_XX^6); % kJ/mol nm^(6)
        Output.C_MX = 4*Output.epsilon_MX*(Output.sigma_MX^6); % kJ/mol nm^(6)
    else
        error('SigmaEpsilon form required for use with JC potential form.')
    end
case 'Mie'
    if Settings.SigmaEpsilon
        if Settings.Additivity
            Output.sigma_MX = (Output.sigma_MM + Output.sigma_XX)/2;
            Output.epsilon_MX = sqrt(Output.epsilon_MM*Output.epsilon_XX); % kJ/mol
            Output.n_MM = Output.n_MX;
            Output.n_XX = Output.n_MX;
        end
        if Settings.Fix_Charge
            Output.Q = Settings.Q_value;
        end

        % Calculate prefactor
        P_MM = (Output.n_MM/(Output.n_MM - 6))*((Output.n_MM/6)^(6/(Output.n_MM - 6)));
        P_XX = (Output.n_XX/(Output.n_XX - 6))*((Output.n_XX/6)^(6/(Output.n_XX - 6)));
        P_MX = (Output.n_MX/(Output.n_MX - 6))*((Output.n_MX/6)^(6/(Output.n_MX - 6)));

        % Convert to Condensed form
        Output.B_MM = P_MM*Output.epsilon_MM*(Output.sigma_MM^Output.n_MM); % kJ/mol nm^(n)
        Output.B_XX = P_XX*Output.epsilon_XX*(Output.sigma_XX^Output.n_XX); % kJ/mol nm^(n)
        Output.B_MX = P_MX*Output.epsilon_MX*(Output.sigma_MX^Output.n_MX); % kJ/mol nm^(n)
        Output.C_MM = P_MM*Output.epsilon_MM*(Output.sigma_MM^6); % kJ/mol nm^(6)
        Output.C_XX = P_XX*Output.epsilon_XX*(Output.sigma_XX^6); % kJ/mol nm^(6)
        Output.C_MX = P_MX*Output.epsilon_MX*(Output.sigma_MX^6); % kJ/mol nm^(6)
    else
        error('SigmaEpsilon form required for use with Mie potential form.')
    end
end


end