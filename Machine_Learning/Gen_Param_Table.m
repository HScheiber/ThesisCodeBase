% This function assumes sigma-epsilon form
function Output = Gen_Param_Table(Settings,Param)

Output = struct;
if strcmp(Settings.Theory,'TF')
    
    % Input parameters
    Output.r0_MM = Param(1); % nm
    Output.r0_XX = Param(2); % nm

    if Settings.Additivity
        Output.r0_MX = (Output.r0_MM + Output.r0_XX)/2; % nm

        Output.epsilon_MM = Param(3); % kJ/mol
        Output.epsilon_XX = Param(4); % kJ/mol
        Output.epsilon_MX = sqrt(Output.epsilon_MM*Output.epsilon_XX); % kJ/mol

        Output.gamma_MX = Param(5); % Unitless
        Output.gamma_MM = Output.gamma_MX; % Unitless
        Output.gamma_XX = Output.gamma_MX; % Unitless
        pidx = 5;
    else
        Output.r0_MX = Param(3); % nm

        Output.epsilon_MM = Param(4); % kJ/mol
        Output.epsilon_XX = Param(5); % kJ/mol
        Output.epsilon_MX = Param(6); % kJ/mol

        Output.gamma_MM = Param(7); % Unitless
        Output.gamma_XX = Param(8); % Unitless
        Output.gamma_MX = Param(9); % Unitless
        pidx = 9;
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
    else
        Output.Q = Param(pidx+1);
    end
    
elseif strcmp(Settings.Theory,'BH') % Buckingham model
    % Input parameters
    Output.r0_MM = Param(1); % nm
    Output.r0_XX = Param(2); % nm

    if Settings.Additivity
        Output.r0_MX = (Output.r0_MM + Output.r0_XX)/2; % nm

        Output.epsilon_MM = Param(3); % kJ/mol
        Output.epsilon_XX = Param(4); % kJ/mol
        Output.epsilon_MX = sqrt(Output.epsilon_MM*Output.epsilon_XX); % kJ/mol

        Output.gamma_MX = Param(5); % Unitless
        Output.gamma_MM = Output.gamma_MX; % Unitless
        Output.gamma_XX = Output.gamma_MX; % Unitless
        pidx = 5;
    else
        Output.r0_MX = Param(3); % nm

        Output.epsilon_MM = Param(4); % kJ/mol
        Output.epsilon_XX = Param(5); % kJ/mol
        Output.epsilon_MX = Param(6); % kJ/mol

        Output.gamma_MM = Param(7); % Unitless
        Output.gamma_XX = Param(8); % Unitless
        Output.gamma_MX = Param(9); % Unitless
        pidx = 9;
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
    else
        Output.Q = Param(pidx+1);
    end
else % JC models

    % D <-> sigma
    % R <-> epsilon
    if Settings.Additivity
        % Sigma scaling
        Output.sigma_MM = Param(1); % nm
        Output.sigma_XX = Param(2); % nm

        % Epsilon scaling
        Output.epsilon_MM = Param(3); % kJ/mol
        Output.epsilon_XX= Param(4); % kJ/mol

        Output.sigma_MX = (Param(1) + Param(2))/2; % nm
        Output.epsilon_MX = sqrt(Param(3)*Param(4)); % kJ/mol

        % Scaling Coulombic Charge
        if Settings.Fix_Charge
            Output.Q = Settings.Q_value;
        else
            Output.Q = Param(5);
        end
    else
        % Sigma scaling
        Output.sigma_MM = Param(1); % nm
        Output.sigma_XX = Param(2); % nm
        Output.sigma_MX = Param(3); % nm

        % Epsilon scaling
        Output.epsilon_MM = Param(4); % kJ/mol
        Output.epsilon_XX= Param(5);  % kJ/mol
        Output.epsilon_MX = Param(6); % kJ/mol

        % Scaling Coulombic Charge
        if Settings.Fix_Charge
            Output.Q = Settings.Q_value;
        else
            Output.Q = Param(7);
        end
    end
    
    % Convert to Condensed form
    Output.B_MM = 4*Output.epsilon_MM*(Output.sigma_MM^12); % kJ/mol nm^(12)
    Output.B_XX = 4*Output.epsilon_XX*(Output.sigma_XX^12); % kJ/mol nm^(12)
    Output.B_MX = 4*Output.epsilon_MX*(Output.sigma_MX^12); % kJ/mol nm^(12)
    Output.C_MM = 4*Output.epsilon_MM*(Output.sigma_MM^6); % kJ/mol nm^(6)
    Output.C_XX = 4*Output.epsilon_XX*(Output.sigma_XX^6); % kJ/mol nm^(6)
    Output.C_MX = 4*Output.epsilon_MX*(Output.sigma_MX^6); % kJ/mol nm^(6)
end


end