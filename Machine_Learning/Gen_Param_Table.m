% This function assumes sigma-epsilon form
function Output = Gen_Param_Table(Settings,Param)

DefParams = bayesopt_params(Settings);
ParNames = {DefParams.Name};
Output = struct;
for idx = 1:numel(Param)
    Output.(ParNames{idx}) = Param(idx);
end


switch Settings.Theory
case 'TF'
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
    
case {'BH' 'BD'} % Buckingham model
    if Settings.Additivity
        Output.r0_MX = (Output.r0_MM + Output.r0_XX)/2; % nm
        Output.epsilon_MX = sqrt(Output.epsilon_MM*Output.epsilon_XX); % kJ/mol
        Output.gamma_MM = Output.gamma_MX; % Unitless
        Output.gamma_XX = Output.gamma_MX; % Unitless
    end

    if Output.gamma_MX < 6
        Output.epsilon_MX = -Output.epsilon_MX;
    end
    if Output.gamma_MM < 6
        Output.epsilon_MM = -Output.epsilon_MM;
    end
    if Output.gamma_XX < 6
        Output.epsilon_XX = -Output.epsilon_XX;
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
case 'JC'

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
case 'Mie'
    
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
end


end