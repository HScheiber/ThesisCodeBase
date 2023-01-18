% Damp_Types:
% 0 = no (default) damping. This is default of JC model.
% 1 = BJ/rational damping (same as in D3(BJ))
% 2 = Tang Damping (Essentially removes dispersion in JC)
% 3 = MMDRE Damping function (Very weak damping, damps mainly at mid range)
% 4 = PAMoC Damping function (fairly weak damping, damps mainly at mid range)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)

% TF Parameter sets for C6/C8 coefficients
% 0 = default TF Parameters
% 1 = D3 values
% 2 = Best literature values available
% 3 = D4 with C6 and C8 generated on the fly

% GAdjust are N x 3 arrays of gaussian parameters
% (i , 1) is the Gaussian height of the ith adjustment (may be negative or
% positive)
% (i , 2) is the center point of the ith Gaussian (should be positive)
% (i , 3) is the standard deviation or width (negative and positive values
% are the same)
% When Model.Parallel_Struct_Min = true, this uses the parallel version of the subroutine Structure_Minimization (Note: each instance of gromacs is single-core in either mode)
%% Inputs
function [Loss,coupledconstraints,UserData] = LiX_Minimizer(Settings,Param,varargin)

% The coupled constraint is for finite temperature simulations that
% give NaN outputs, or for Structure_Minimization calculations that are not feasible. 
% These are caused when 
%   (1) The lattice params of a structure are too large or too small
%   (2) The lattice energy of a structure is lower than Settings.MinMDP.E_Unphys
%   (3) A finite-T simulation fails 
%   (4) A melting point cannot be found for the structure of interest 
%   (5) The liquid is amorphous at the experimental MP 
%   (6) The liquid or solid converts to another structure at the experimental MP

if Settings.UseCoupledConstraint
    coupledconstraints = -1;
else
    coupledconstraints = [];
end

% Optional inputs
p = inputParser;
p.FunctionName = 'LiX_Minimizer';
addOptional(p,'Extra_Properties',Settings.Extra_Properties,@(x)validateattributes(x,{'logical'},{'nonempty'}))
addOptional(p,'Therm_Prop_Override',false,@(x)validateattributes(x,{'logical'},{'nonempty'}))

parse(p,varargin{:});
Extra_Properties = p.Results.Extra_Properties;
Settings.Therm_Prop_Override = p.Results.Therm_Prop_Override;

Settings.MinMDP.Verbose = Settings.Verbose;
if Settings.Parallel_LiX_Minimizer && Settings.Parallel_Struct_Min
    Settings.Parallel_Struct_Min = false;
end

% Convert Params to table form
if ~istable(Param) && ~isstruct(Param)
    List_par = bayesopt_params(Settings);
    ParNames = {List_par.Name}; % gives the order of names
    ParamNums = Param;
    Param = table;
    for idx = 1:numel(ParNames)
        Param.(ParNames{idx}) = ParamNums(idx);
    end
else
    ParNames = Param.Properties.VariableNames;
end

% Coupled or each salt independent?
[Salts,MultiMetal,MultiHalide] = Select_MultiSalt(Settings);

Loss_Salt = 0;
for salt_idx = 1:numel(Salts) % Loop through coupled salts
    Settings.Salt = Salts{salt_idx};
    [Settings.Metal,Settings.Halide] = Separate_Metal_Halide(Settings.Salt);
    Settings.Minimization_Data = Initialize_Minimization_Data(Settings);
    Settings.Finite_T_Data = Initialize_Finite_T_Data(Settings);
    
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

    % Potential Scaling
    switch Settings.Theory
    case 'TF'
        % Loose form of exp-C6-C8 model
        if Settings.SigmaEpsilon

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

            if gamma_MX > 48/7
                epsilon_MX = -epsilon_MX;
            end
            if gamma_MM > 48/7
                epsilon_MM = -epsilon_MM;
            end
            if gamma_XX > 48/7
                epsilon_XX = -epsilon_XX;
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

        % Tight form of exp-C6-C8 model
        else

            % 1/R6 Dispersion (TF only)
            Settings.S.D6D.MM = Param.SD6MM;
            Settings.S.D6D.XX = Param.SD6XX;
            Settings.S.D6D.MX = Param.SD6MX;

            % 1/R8 Dispersion (TF only)
            if Settings.Fix_C8
                % Calculate value of C8 using recursive relations

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
                C8.MM = 3.0*(Settings.S.D6D.MM/c6units)*sqrt_Q.(Settings.Metal)*sqrt_Q.(Settings.Metal)*c8units; % in kJ/mol nm^8
                C8.XX = 3.0*(Settings.S.D6D.XX/c6units)*sqrt_Q.(Settings.Halide)*sqrt_Q.(Settings.Halide)*c8units; % in kJ/mol nm^8
                C8.MX = 3.0*(Settings.S.D6D.MX/c6units)*sqrt_Q.(Settings.Metal)*sqrt_Q.(Settings.Halide)*c8units; % in kJ/mol nm^8

                % Update the scaling
                Settings.S.D8D.MM = C8.MM;
                Settings.S.D8D.XX = C8.XX;
                Settings.S.D8D.MX = C8.MX;
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
                        r0_MX = (r0_MM + r0_XX)/2; % nm
                        epsilon_MX = sqrt(epsilon_MM*epsilon_XX); % kJ/mol
                        gamma_MM = gamma_MX; % Unitless
                        gamma_XX = gamma_MX; % Unitless
                        if gamma_MX < 6 && epsilon_MX > 0
                            epsilon_MX = -epsilon_MX;
                        end
                        if gamma_MM < 6 && epsilon_MM > 0
                            epsilon_MM = -epsilon_MM;
                        end
                        if gamma_XX < 6 && epsilon_XX > 0
                            epsilon_XX = -epsilon_XX;
                        end

                    case 'hogervorst'
                        gamma_MM = Param.gamma_MM; % Unitless
                        gamma_XX = Param.gamma_XX; % Unitless
                        gamma_MX = (gamma_MM + gamma_XX)/2;

                        if gamma_MM < 6 && epsilon_MM > 0
                            epsilon_MM = -epsilon_MM;
                        end
                        if gamma_XX < 6 && epsilon_XX > 0
                            epsilon_XX = -epsilon_XX;
                        end

                        epsilon_MX = 2*epsilon_MM*epsilon_XX/(epsilon_MM + epsilon_XX);
                        if gamma_MX < 6 && epsilon_MX > 0
                            epsilon_MX = -epsilon_MX;
                        end

                        r0_MX = ( sqrt( ( epsilon_MM*epsilon_XX*gamma_MM*gamma_XX*(r0_MM*r0_XX)^6 )...
                            /((gamma_MM - 6)*(gamma_XX - 6)) )*(gamma_MX - 6)/(epsilon_MX*gamma_MX) )^(1/6);

                    case {'kong' 'gromacs'}
                        gamma_MM = Param.gamma_MM; % Unitless
                        gamma_XX = Param.gamma_XX; % Unitless

                        if gamma_MM < 6 && epsilon_MM > 0
                            epsilon_MM = -epsilon_MM;
                        end
                        if gamma_XX < 6 && epsilon_XX > 0
                            epsilon_XX = -epsilon_XX;
                        end

                        k_MM = 1 / ( gamma_MM - 6 );
                        k_XX = 1 / ( gamma_XX - 6 );

                        A_MM = 6*epsilon_MM*k_MM*exp(gamma_MM); % prefactor
                        A_XX = 6*epsilon_XX*k_XX*exp(gamma_XX);

                        B_MM = gamma_MM/r0_MM; % exponent
                        B_XX = gamma_XX/r0_XX;

                        C_MM = epsilon_MM*gamma_MM*k_MM*(r0_MM^6); % dispersion
                        C_XX = epsilon_XX*gamma_XX*k_XX*(r0_XX^6);

                        switch lower(Settings.Comb_rule)
                            case 'kong'
                                A_MX = (1/2)*( A_MM*(A_MM*B_MM/(A_XX*B_XX))^(-B_MM/(B_MM + B_XX)) + ...
                                               A_XX*(A_XX*B_XX/(A_MM*B_MM))^(-B_XX/(B_MM + B_XX)) );
                                B_MX = 2*B_MM*B_XX/(B_MM + B_XX);
                                C_MX = sqrt(C_MM*C_XX);
                            case 'gromacs'
                                A_MX = sqrt(A_MM*A_XX);
                                B_MX = 2/( (1/B_MM) + (1/B_XX) );
                                C_MX = sqrt(C_MM*C_XX);
                        end

                        % Convert back to gamma/epsilon/r0
                        gamma_MX   = -7*lambertw((-1/7)*(6*C_MX*(B_MX^6)/A_MX)^(1/7));
                        r0_MX      = gamma_MX/B_MX;
                        epsilon_MX = C_MX*(gamma_MX - 6)/(gamma_MX*(r0_MX^6));
                        if gamma_MX < 6 && epsilon_MX > 0
                            epsilon_MX = -epsilon_MX;
                        end
                end
            else
                r0_MX = Param.r0_MX; % nm
                epsilon_MX = Param.epsilon_MX; % kJ/mol
                gamma_MX = Param.gamma_MX; % Unitless
                gamma_MM = Param.gamma_MM; % Unitless
                gamma_XX = Param.gamma_XX; % Unitless
                if gamma_MX < 6 && epsilon_MX > 0
                    epsilon_MX = -epsilon_MX;
                end
                if gamma_MM < 6 && epsilon_MM > 0
                    epsilon_MM = -epsilon_MM;
                end
                if gamma_XX < 6 && epsilon_XX > 0
                    epsilon_XX = -epsilon_XX;
                end
            end

            % Check for nonsense
            if abs(imag(r0_MX)) > sqrt(eps)
                r0_MX = nan;
            end
            if abs(imag(gamma_MX)) > sqrt(eps)
                gamma_MX = nan;
            end
            if abs(imag(epsilon_MX)) > sqrt(eps)
                epsilon_MX = nan;
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
                        sigma_MX = (sigma_MM + sigma_XX)/2; % nm
                        epsilon_MX = sqrt(epsilon_MM*epsilon_XX); % kJ/mol
                        gamma_MM = gamma_MX; % Unitless
                        gamma_XX = gamma_MX; % Unitless
                    case 'hogervorst'
                        gamma_MM = Param.gamma_MM; % Unitless
                        gamma_XX = Param.gamma_XX; % Unitless
                        gamma_MX = (gamma_MM + gamma_XX)/2;

                        epsilon_MX = 2*epsilon_MM*epsilon_XX/(epsilon_MM + epsilon_XX);

                        sigma_MX = ( sqrt( ( epsilon_MM*epsilon_XX*gamma_MM*gamma_XX*(sigma_MM*sigma_XX)^6 )...
                            /((gamma_MM - 6)*(gamma_XX - 6)) )*(gamma_MX - 6)/(epsilon_MX*gamma_MX) )^(1/6);
                    case 'gromacs'
                        gamma_MM = Param.gamma_MM; % Unitless
                        gamma_XX = Param.gamma_XX; % Unitless
                        gamma_MX = sqrt(gamma_MM.*gamma_XX);

                        epsilon_MX = sqrt(epsilon_MM.*epsilon_XX);
                        sigma_MX = sqrt(sigma_MM.*sigma_XX);
                end
        else
            sigma_MX = Param.sigma_MX; % nm
            epsilon_MX = Param.epsilon_MX; % kJ/mol
            gamma_MM = Param.gamma_MM; % Unitless
            gamma_XX = Param.gamma_XX; % Unitless
            gamma_MX = Param.gamma_MX; % Unitless
        end

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
    case 'JC' % JC models

        % sigma/epsilon form (cast in terms of sigma/epsilon scaling internally)
        if Settings.SigmaEpsilon
            PotSettings = Initialize_MD_Settings;
            PotSettings.Salt = Settings.Salt;
            [MXParams,MMParams,XXParams] = JC_Potential_Parameters(PotSettings);

            % Sigma scaling
            Settings.S.S.MM = Param.sigma_MM/MMParams.sigma;
            Settings.S.S.XX = Param.sigma_XX/XXParams.sigma;

            % Epsilon scaling
            Settings.S.E.MM = Param.epsilon_MM/MMParams.epsilon;
            Settings.S.E.XX = Param.epsilon_XX/XXParams.epsilon;

            % Default MX params
            def_S_MX = MXParams.sigma;
            def_E_MX = MXParams.epsilon;

            if Settings.Additivity
                Sigma_MX = (Param.sigma_MM + Param.sigma_XX)/2;
                Epsilon_MX = sqrt(Param.epsilon_MM*Param.epsilon_XX);

                Settings.S.S.MX = Sigma_MX/def_S_MX;
                Settings.S.E.MX = Epsilon_MX/def_E_MX;

                if Settings.Additional_MM_Disp
                    Full_MM_Epsilon = Param.epsilon_MM + Param.epsilon_MM2;
                    Settings.S.E.MM = Full_MM_Epsilon/MMParams.epsilon;
                end
            else
                Settings.S.S.MX = Param.sigma_MX/def_S_MX;
                Settings.S.E.MX = Param.epsilon_MX/def_E_MX;
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
    case 'LJ' % LJ models

        % sigma/epsilon form (cast in terms of sigma/epsilon scaling internally)
        if Settings.SigmaEpsilon

            % Sigma scaling
            Settings.S.S.MM = Param.sigma_MM;
            Settings.S.S.XX = Param.sigma_XX;

            % Epsilon scaling
            Settings.S.E.MM = Param.epsilon_MM;
            Settings.S.E.XX = Param.epsilon_XX;
            
            if Settings.Additivity
                
                Settings.S.S.MX = (Param.sigma_MM + Param.sigma_XX)/2;
                Settings.S.E.MX = sqrt(Param.epsilon_MM*Param.epsilon_XX);

                if Settings.Additional_MM_Disp
                    Settings.S.E.MM = Param.epsilon_MM + Param.epsilon_MM2;
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
                % Scaled
                MM_Epsilon = (Settings.S.D.MM^2)*(1/Settings.S.R.MM);
                MM_Sigma = (1/(Settings.S.D.MM^(1/6)))*(Settings.S.R.MM^(1/6));

                XX_Epsilon = (Settings.S.D.XX^2)*(1/Settings.S.R.XX);
                XX_Sigma = (1/(Settings.S.D.XX^(1/6)))*(Settings.S.R.XX^(1/6));

                MX_Epsilon = sqrt(MM_Epsilon*XX_Epsilon);
                MX_Sigma   = (MM_Sigma + XX_Sigma)/2;

                Settings.S.D.MX = 4*MX_Epsilon*MX_Sigma^6;
                Settings.S.R.MX = 4*MX_Epsilon*MX_Sigma^12;

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

    % Calculate loss due to infeasible models
    Loss_add = 0;
    if Settings.CheckBadFcn  
        tl = Settings.Table_Length;
        ss = Settings.Table_StepSize;
        Settings.Table_Length = 10; % nm
        Settings.Table_StepSize = 0.01;
        if strcmp(Settings.Theory,'BH')
            [U,~] = BH_Potential_Generator(Settings,'Include_Dispersion_Scale',true,...
                'Startpoint',0.01);
        elseif strcmp(Settings.Theory,'TF')
            [U,~] = TF_Potential_Generator(Settings,'Include_Dispersion_Scale',true,...
                'Startpoint',0.01);
        elseif strcmp(Settings.Theory,'JC')
            [U,~] = JC_Potential_Generator(Settings,'Include_Dispersion_Scale',true,...
                'Startpoint',0.01);
        elseif strcmp(Settings.Theory,'LJ')
            [U,~] = LJ_Potential_Generator(Settings,'Include_Dispersion_Scale',true,...
                'Startpoint',0.01);
        elseif strcmp(Settings.Theory,'Mie')
            [U,~] = Mie_Potential_Generator(Settings,'Include_Dispersion_Scale',true,...
                'Startpoint',0.01);
        elseif strcmp(Settings.Theory,'BD')
            [U,~] = BD_Potential_Generator(Settings,'Include_Dispersion_Scale',true,...
                'Startpoint',0.01);
        elseif strcmp(Settings.Theory,'BE')
            [U,~] = BE_Potential_Generator(Settings,'Include_Dispersion_Scale',true,...
                'Startpoint',0.01);
        elseif strcmp(Settings.Theory,'BF')
            [U,~] = BF_Potential_Generator(Settings,'Include_Dispersion_Scale',true,...
                'Startpoint',0.01);
        end
        Settings.Table_Length = tl; % nm
        Settings.Table_StepSize = ss;

        %% Grab the peaks and valleys of the MX attractive potential

        if any(isnan(U.MX.Total),2)
            Loss_add = Loss_add + Settings.BadFcnLossPenalty;
        end

        peaks_idx = islocalmax(U.MX.Total,'MinProminence',1e-8);
        valleys_idx = islocalmin(U.MX.Total,'MinProminence',1e-8);

        U_peak = U.MX.Total(peaks_idx);
        U_valley = U.MX.Total(valleys_idx);
        r_peak = U.r(peaks_idx);
        if numel(r_peak) > 1
            r_peak = r_peak(end);
        end
        r_valley = U.r(valleys_idx);
        U_valley(r_valley<=r_peak) = [];
        r_valley(r_valley<=r_peak) = [];

        if isempty(U_valley) % If no well minimum exists in MX interaction
            Loss_add = Loss_add + Settings.BadFcnLossPenalty;
        elseif isempty(U_peak) % If no peak exists in MX interaction, but valley does, check well depth. 
            % This is normal for JC potential and some BH/TF potentials
            Loss_add = Loss_add + max(Settings.MaxAttWellDepth - U_valley,0)*Settings.BadFcnLossPenalty;
            % Also check valley location
            Loss_add = Loss_add + max(r_valley - (Settings.MaxMXWellR/10),0).*Settings.BadFcnLossPenalty;
            Loss_add = Loss_add + max((Settings.MinMXWellR/10) - r_valley,0).*Settings.BadFcnLossPenalty;

            % Also check well location
        else % Otherwise, a well minimum exists and at least one peak exists
            % Penalize wells that are too shallow and wells that are too deep
            dU = U_peak - U_valley;
            Loss_add = Loss_add + max(Settings.MinExpWallHeight - dU,0)*Settings.BadFcnLossPenalty;
            Loss_add = Loss_add + max(Settings.MaxAttWellDepth - U_valley,0)*Settings.BadFcnLossPenalty;
            % Also check valley location
            Loss_add = Loss_add + max(r_valley - (Settings.MaxMXWellR/10),0).*Settings.BadFcnLossPenalty;
            Loss_add = Loss_add + max((Settings.MinMXWellR/10) - r_valley,0).*Settings.BadFcnLossPenalty;
        end
    %     plot(U.r,U.MX.Total)
    %     hold on
    %     scatter(U.r(peaks_idx),U.MX.Total(peaks_idx))
    %     scatter(U.r(valleys_idx),U.MX.Total(valleys_idx))
    %     ylim([-1000 1000])

        %% Grab the peaks and valleys of the MM/XX potentials
        r_wall_sv = nan(2,1);
        ints = {'MM' 'XX'};
        for idx = 1:numel(ints)
            int = ints{idx};
            peaks_idx = islocalmax(U.(int).Total,'MinProminence',1e-8);
            valleys_idx = islocalmin(U.(int).Total,'MinProminence',1e-8);

            U_peak = U.(int).Total(peaks_idx);
            U_valley = U.(int).Total(valleys_idx);

            r_peak = U.r(peaks_idx);
            r_valley = U.r(valleys_idx);

    %         plot(U.r,U.(int).Total)
    %         hold on
    %         scatter(U.r(peaks_idx),U.(int).Total(peaks_idx))
    %         scatter(U.r(valleys_idx),U.(int).Total(valleys_idx))
    %         ylim([-1000 1000])

            if isempty(U_peak) % No peak exists
                % Do nothing, this is normal for JC and sometimes BH/TF
            elseif length(U_peak) == 1 && isempty(U_valley) % One peak, no valley 
                % this is a normal case for TF and BH repulsive potentials
                % check to ensure the peak height is at least Settings.MinExpWallHeight
                Loss_add = Loss_add + max(Settings.MinExpWallHeight - U_peak,0).*Settings.BadFcnLossPenalty;

            elseif length(U_peak) > 1 % two peaks are visible + one valley in between them
                % Penalize any model with a non-zero well depth between
                % like-like interactions
                dU = U_peak(2) - U_valley;
                Loss_add = Loss_add + max(dU - Settings.MaxRepWellDepth,0)*Settings.BadFcnLossPenalty; % valley exists
                dU_fp = U_peak(1) - U_valley;
                Loss_add = Loss_add + max(Settings.MinExpWallHeight - dU_fp,0).*Settings.BadFcnLossPenalty; % peak too low

            elseif length(U_peak) == 1 && length(U_valley) == 1 % One peak visible + one valley in between, and (possibly) one hidden peak to the right

                if r_valley < r_peak % valley closer set
                    dU = U_peak - U_valley; % valley exists
                    Loss_add = Loss_add + max(dU - Settings.MaxRepWellDepth,0)*Settings.BadFcnLossPenalty;
                else % peak closer set
                    % Case of hidden peak to the right
                    % Ensure valley depth is not greater than the threshold
                    dU = U.(int).Total(end) - U_valley;
                    Loss_add = Loss_add + max(dU - Settings.MaxRepWellDepth,0)*Settings.BadFcnLossPenalty;

                    % Also check the repulsive peak height
                    dUfp = U_peak - U_valley;
                    Loss_add = Loss_add + max(Settings.MinExpWallHeight - dUfp,0).*Settings.BadFcnLossPenalty; % wall too low
                end

            elseif length(U_valley) == 1 % valley exists but no peaks are visible, there must be a hidden peak to the right
                dU = U.(int).Total(end) - U_valley;
                Loss_add = Loss_add + max(dU - Settings.MaxRepWellDepth,0)*Settings.BadFcnLossPenalty;
            else
                % This should never be reached...
                warning('Possible issue with the potential!')
            end

            % Check that the repulsive wall is not too far out or close in
            if any(U.(int).Total >= Settings.MinExpWallHeight)
                r_above = U.r(U.(int).Total >= Settings.MinExpWallHeight);
                r_wall = r_above(end);
                r_wall_sv(idx) = r_wall;
                Loss_add = Loss_add + max(r_wall - (Settings.MaxMXWellR/10),0).*Settings.BadFcnLossPenalty; % wall too far
                Loss_add = Loss_add + max((Settings.MinMXWellR/10) - r_wall,0).*Settings.BadFcnLossPenalty; % wall too close
            end
        end

        % Check that the repulsive wall of the XX interaction is further out than the MM interaction
        if Settings.EnforceRR && Settings.SigmaEpsilon && Settings.Additivity
            switch Settings.Theory
                case {'JC' 'Mie' 'BF' 'LJ'}
                    RR_Walls = Param.sigma_MM./Param.sigma_XX; % Ratio of M/X size, should be < 1
                case {'BH' 'BD' 'BE'}
                    RR_Walls = Param.r0_MM./Param.r0_XX; % Ratio of M/X size, should be < 1
            end
            Loss_add = Loss_add + max(RR_Walls - 1,0).*Settings.BadFcnLossPenalty; % Radius ratios are incorrect
        elseif Settings.EnforceRR
            RR_Walls = r_wall_sv(1)./r_wall_sv(2); % Ratio of M/X size, should be < 1
            Loss_add = Loss_add + max(RR_Walls - 1,0).*Settings.BadFcnLossPenalty; % Radius ratios are incorrect
        end
    end

    %% Parallel Setup
    N = length(Settings.Structures);
    Structure_Min_Calc_Fail = false;
    if Settings.Parallel_LiX_Minimizer
        % Set up matlab parallel features
        Parcores = feature('numcores');
        PrefCores = min(Parcores,N);
        if ~isempty(gcp('nocreate'))
            Cur_Pool = gcp;
            Cur_Workers = Cur_Pool.NumWorkers;

            % Start the parallel pool
            if Cur_Workers < PrefCores
                delete(Cur_Pool);
                % Create a "local" cluster object
                local_cluster = parcluster('local');

                % Modify the JobStorageLocation to a temporary directory
                [~,~,computer] = find_home;
                switch computer
                    case {'cedar' 'graham' 'narval'}
                        tmp = fullfile(getenv('SLURM_TMPDIR'),'local_cluster_jobs');
                    case 'sockeye'
                        tmp = fullfile(getenv('TMPDIR'),'local_cluster_jobs');
                    otherwise
                        tmp = fullfile(tempname,'local_cluster_jobs');
                end
                if ~isfolder(tmp)
                    mkdir(tmp)
                end
                local_cluster.JobStorageLocation = tmp;
                ppool = parpool(local_cluster,min(Parcores,N));
            else
                ppool = Cur_Pool;
            end
        else
            ppool = parpool(min(Parcores,N)); 
        end

        f(1:N) = parallel.FevalFuture;

        % Run in parallel with parfavel
        for idx = 1:N
            Settings.Structure = Settings.Structures{idx};
            f(idx) = parfeval(ppool,@Structure_Minimization,1,Settings,...
                'Extra_Properties',Extra_Properties);
        end

        wait(f); % Wait for parallel jobs to finish

        % Collect outputs into cell array
        for idx = 1:N
            Settings.Minimization_Data{idx} = f(idx).fetchOutputs;
            if Settings.Minimization_Data{idx}.CalcFail
                Structure_Min_Calc_Fail = true;
            end
        end
    %% Serial mode    
    else
        for idx = 1:N
            Settings.Structure = Settings.Structures{idx};
            Settings.Minimization_Data{idx} = Structure_Minimization(Settings,...
                'Extra_Properties',Extra_Properties);
            if Settings.Minimization_Data{idx}.CalcFail
                Structure_Min_Calc_Fail = true;
            end
        end
    end

    % This catches Structure_Minimization calculations that produce a result
    % outside of the allowed energy or volume bounds
    if Structure_Min_Calc_Fail && ~Settings.Therm_Prop_Override
        if Settings.UseCoupledConstraint
            Loss = nan;
        else
            Loss = real(log1p(Loss_Salt + Loss_add + Settings.BadFcnLossPenalty));
        end
        if MultiMetal || MultiHalide
            UserData.Finite_T_Data.(Settings.Salt) = Settings.Finite_T_Data;
            UserData.Minimization_Data.(Settings.Salt) = Settings.Minimization_Data;
        else
            UserData.Finite_T_Data = Settings.Finite_T_Data;
            UserData.Minimization_Data = Settings.Minimization_Data;
        end
        return
    elseif Structure_Min_Calc_Fail && Settings.Therm_Prop_Override
        if Settings.UseCoupledConstraint
            Loss_add = nan;
        else
            Loss_add = Loss_add + Settings.BadFcnLossPenalty;
        end
    end

    % Initialize Finite T Data structure and update Settings
    tol = sqrt(eps);
    Settings.skip_finite_T = false;
    if any([Settings.Loss_Options.Fusion_Enthalpy ...
            Settings.Loss_Options.MP_Volume_Change ...
            Settings.Loss_Options.Liquid_MP_Volume ...
            Settings.Loss_Options.Solid_MP_Volume ...
            Settings.Loss_Options.MP] > tol) || Settings.Therm_Prop_Override

        Settings.T0 = Settings.Finite_T_Data.Exp_MP; % K, Initial temperature
        Settings.Target_T = Settings.Finite_T_Data.Exp_MP; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings.MDP.Initial_T = Settings.Finite_T_Data.Exp_MP; % Initial termpature at which to generate velocities

        Settings.Structure = Settings.Finite_T_Data.Structure;
        Settings.Geometry = Default_Crystal(Settings,'Center_Coordinates',true);

        % Get an estimate for the density in molecules/nm^3
        N = length(Settings.Minimization_Data);
        Structures = cell(1,N);
        for idx = 1:N
            Structures{idx} = Settings.Minimization_Data{idx}.Structure;
        end
        strmatch = strcmp(Settings.Finite_T_Data.Structure,Structures);
        V0_model = Settings.Minimization_Data{strmatch}.V; % Volume of model in A^3/molecule

        if ~isfield(Settings,'MaxModelVolume')
            defSettings = Initialize_LiX_BO_Settings;
            Settings.MaxModelVolume = defSettings.MaxModelVolume;
        end

        if V0_model > Settings.MaxModelVolume
            Model_Mismatch = (V0_model - Settings.MaxModelVolume)/Settings.MaxModelVolume;
            Loss_add_Vol = Model_Mismatch*Settings.BadFcnLossPenalty;
            if Settings.UseCoupledConstraint
                coupledconstraints = real(log1p(Loss_add_Vol));
            end
        elseif V0_model < Settings.MinModelVolume
            Model_Mismatch = (Settings.MinModelVolume - V0_model)/Settings.MinModelVolume;
            Loss_add_Vol = Model_Mismatch*Settings.BadFcnLossPenalty;
            if Settings.UseCoupledConstraint
                coupledconstraints = real(log1p(Loss_add_Vol));
            end
        else
            Loss_add_Vol = 0;
        end

        Loss_add = Loss_add + Loss_add_Vol;
        if real(log1p(Loss_add)) >= Settings.MinSkipLoss && ~Settings.Therm_Prop_Override
            Settings.skip_finite_T = true;
        else
            Settings.Ref_Density = 1/(Settings.Minimization_Data{strmatch}.V*(0.1^3)); % molecules / nm^3
            Settings.Geometry.a = Settings.Minimization_Data{strmatch}.a;
            Settings.Geometry.b = Settings.Minimization_Data{strmatch}.b;
            Settings.Geometry.c = Settings.Minimization_Data{strmatch}.c;
        end
    end

    % Melting point
    if ( Settings.Loss_Options.MP > tol && ~Settings.skip_finite_T ) || Settings.Therm_Prop_Override

        if Settings.Parallel_Bayesopt % Run serially
            env.OMP_NUM_THREADS = getenv('OMP_NUM_THREADS');
            env.GMX_PME_NUM_THREADS = getenv('GMX_PME_NUM_THREADS');
            env.GMX_PME_NTHREADS = getenv('GMX_PME_NTHREADS');
            env.GMX_OPENMP_MAX_THREADS = getenv('GMX_OPENMP_MAX_THREADS');
            env.KMP_AFFINITY = getenv('KMP_AFFINITY');
            setenv('OMP_NUM_THREADS','1');
            setenv('GMX_PME_NUM_THREADS','1');
            setenv('GMX_PME_NTHREADS','1');
            setenv('GMX_OPENMP_MAX_THREADS','1');
            setenv('KMP_AFFINITY','disabled');
            Settings.mdrun_opts = ' -pin on -ntmpi 1 -ntomp 1';
            Settings.gmx = Settings.gmx_loc;
        elseif ~isempty(gcp('nocreate')) % close the current ppool, it will likely close on its own anyway, causing issues
            delete(gcp);
        end

        [WorkDir,Settings.JobName,Settings.Full_Model_Name] = GetMDWorkdir(Settings);
        Settings.WorkDir = [WorkDir '_MP'];
        ThermFolder = fullfile(Settings.OuterDir,'BestPoint_Thermal','Melting_Point');
        if Settings.Therm_Prop_Override && isfolder(ThermFolder)
            Settings.WorkDir = ThermFolder;
        end

        Settings.BatchMode = false;
        Settings.Submit_Jobs = false;
        Settings.Skip_Minimization = true; % Skip the automatic geometry minimization
        Settings.RefStructure = Settings.Finite_T_Data.Structure;
        Verbose = Settings.Verbose;
        Settings.Verbose = true;
        [Tm_estimate,~,Aborted,T_dat] = Find_Melting_Point(Settings);
        Settings.Verbose = Verbose;

        if Settings.Therm_Prop_Override && ~isfolder(ThermFolder)
            copyfile(Settings.WorkDir,ThermFolder)
        end

        Settings.Finite_T_Data.T_dat = T_dat;
        if Aborted
            Settings.Finite_T_Data.MP = nan;
            Loss_add = nan;
        else
            Settings.Finite_T_Data.MP = Tm_estimate;
            if Settings.Delete_Equil && isfolder(Settings.WorkDir)
                try
                    rmdir(Settings.WorkDir,'s')
                catch
                    if Settings.Verbose
                        warning(['Unable to remove directory: ' Settings.WorkDir])
                    end
                end
            end
        end
        if Settings.Parallel_Bayesopt
            setenv('OMP_NUM_THREADS',env.OMP_NUM_THREADS);
            setenv('GMX_PME_NUM_THREADS',env.GMX_PME_NUM_THREADS);
            setenv('GMX_PME_NTHREADS',env.GMX_PME_NTHREADS);
            setenv('GMX_OPENMP_MAX_THREADS',env.GMX_OPENMP_MAX_THREADS);
            setenv('KMP_AFFINITY',env.KMP_AFFINITY);
        end
    end

    % High T liquid properties
    if ( any([Settings.Loss_Options.Fusion_Enthalpy ...
            Settings.Loss_Options.MP_Volume_Change ...
            Settings.Loss_Options.Liquid_MP_Volume ...
            Settings.Loss_Options.Liquid_DM_MP] > tol) ...
            && ~Settings.skip_finite_T ) || Settings.Therm_Prop_Override

        [WorkDir,Settings.JobName,Settings.Full_Model_Name] = GetMDWorkdir(Settings);
        Settings.WorkDir = [WorkDir '_LP'];

        ThermFolder = fullfile(Settings.OuterDir,'BestPoint_Thermal','Liq_Properties_at_MP');
        if Settings.Therm_Prop_Override && isfolder(ThermFolder)
            Settings.WorkDir = ThermFolder;
        end

        if Settings.Parallel_Bayesopt % Run serially
            env.OMP_NUM_THREADS = getenv('OMP_NUM_THREADS');
            env.GMX_PME_NUM_THREADS = getenv('GMX_PME_NUM_THREADS');
            env.GMX_PME_NTHREADS = getenv('GMX_PME_NTHREADS');
            env.GMX_OPENMP_MAX_THREADS = getenv('GMX_OPENMP_MAX_THREADS');
            env.KMP_AFFINITY = getenv('KMP_AFFINITY');
            setenv('OMP_NUM_THREADS','1');
            setenv('GMX_PME_NUM_THREADS','1');
            setenv('GMX_PME_NTHREADS','1');
            setenv('GMX_OPENMP_MAX_THREADS','1');
            setenv('KMP_AFFINITY','disabled');
            Settings.mdrun_opts = ' -pin on -ntmpi 1 -ntomp 1';
            Settings.gmx = Settings.gmx_loc;
            Verbose = Settings.Verbose;
            Settings.Verbose = true;
            Liq_Output = Calc_Liquid_Properties_at_MP(Settings); % Output is nan if liquid converts to >0.85 solid
            Settings.Verbose = Verbose;
        else
            dd = Settings.dd;
            npme = Settings.npme;
            Settings.dd = [];
            Settings.npme = [];
            [~,Settings] = MD_Batch_Template(Settings);
            Verbose = Settings.Verbose;
            Settings.Verbose = true;
            Liq_Output = Calc_Liquid_Properties_at_MP(Settings); % Output is nan if liquid converts to >0.85 solid
            Settings.Verbose = Verbose;
            Settings.dd = dd;
            Settings.npme = npme;
            [~,Settings] = MD_Batch_Template(Settings);
        end

        if Settings.Therm_Prop_Override && ~isfolder(ThermFolder)
            copyfile(Settings.WorkDir,ThermFolder)
        end

        Settings.Finite_T_Data.Liquid_V_MP = Liq_Output.Liquid_V_MP;
        Settings.Finite_T_Data.Liquid_H_MP = Liq_Output.Liquid_H_MP;
        Settings.Finite_T_Data.Liquid_DM_MP = Liq_Output.Liquid_DM_MP; % cm^2 / s
        Settings.Finite_T_Data.Liquid_DX_MP = Liq_Output.Liquid_DX_MP;

        if isnan(Liq_Output.Liquid_H_MP)
            Loss_add = nan;
            if ~Settings.Therm_Prop_Override
                Loss = nan;
                if MultiMetal || MultiHalide
                    UserData.Finite_T_Data.(Settings.Salt) = Settings.Finite_T_Data;
                    UserData.Minimization_Data.(Settings.Salt) = Settings.Minimization_Data;
                else
                    UserData.Finite_T_Data = Settings.Finite_T_Data;
                    UserData.Minimization_Data = Settings.Minimization_Data;
                end
                return
            end
        end

        if Settings.Parallel_Bayesopt
            setenv('OMP_NUM_THREADS',env.OMP_NUM_THREADS);
            setenv('GMX_PME_NUM_THREADS',env.GMX_PME_NUM_THREADS);
            setenv('GMX_PME_NTHREADS',env.GMX_PME_NTHREADS);
            setenv('GMX_OPENMP_MAX_THREADS',env.GMX_OPENMP_MAX_THREADS);
            setenv('KMP_AFFINITY',env.KMP_AFFINITY);
        end
    end

    % High T solid properties
    if ( any([Settings.Loss_Options.Fusion_Enthalpy ...
            Settings.Loss_Options.MP_Volume_Change ...
            Settings.Loss_Options.Solid_MP_Volume] > tol) ...
            && ~Settings.skip_finite_T ) || Settings.Therm_Prop_Override

        [WorkDir,Settings.JobName,Settings.Full_Model_Name] = GetMDWorkdir(Settings);
        Settings.WorkDir = [WorkDir '_SP'];

        ThermFolder = fullfile(Settings.OuterDir,'BestPoint_Thermal','Sol_Properties_at_MP');
        if Settings.Therm_Prop_Override && isfolder(ThermFolder)
            Settings.WorkDir = ThermFolder;
        end

        if Settings.Parallel_Bayesopt % Run serially
            env.OMP_NUM_THREADS = getenv('OMP_NUM_THREADS');
            env.GMX_PME_NUM_THREADS = getenv('GMX_PME_NUM_THREADS');
            env.GMX_PME_NTHREADS = getenv('GMX_PME_NTHREADS');
            env.GMX_OPENMP_MAX_THREADS = getenv('GMX_OPENMP_MAX_THREADS');
            env.KMP_AFFINITY = getenv('KMP_AFFINITY');
            setenv('OMP_NUM_THREADS','1');
            setenv('GMX_PME_NUM_THREADS','1');
            setenv('GMX_PME_NTHREADS','1');
            setenv('GMX_OPENMP_MAX_THREADS','1');
            setenv('KMP_AFFINITY','disabled');
            Settings.mdrun_opts = ' -pin on -ntmpi 1 -ntomp 1';
            Settings.gmx = Settings.gmx_loc;
            Verbose = Settings.Verbose;
            Settings.Verbose = true;
            Sol_Output = Calc_Solid_Properties_at_MP(Settings);
            Settings.Verbose = Verbose;
        else
            dd = Settings.dd;
            npme = Settings.npme;
            Settings.dd = [];
            Settings.npme = [];
            [~,Settings] = MD_Batch_Template(Settings);
            Verbose = Settings.Verbose;
            Settings.Verbose = true;
            Sol_Output = Calc_Solid_Properties_at_MP(Settings);
            Settings.Verbose = Verbose;
            Settings.dd = dd;
            Settings.npme = npme;
            [~,Settings] = MD_Batch_Template(Settings);
        end

        if Settings.Therm_Prop_Override && ~isfolder(ThermFolder)
            copyfile(Settings.WorkDir,ThermFolder)
        end
        Settings.Finite_T_Data.Solid_V_MP = Sol_Output.Solid_V_MP;
        Settings.Finite_T_Data.Solid_H_MP = Sol_Output.Solid_H_MP;

        Settings.Finite_T_Data.Fusion_dH = Settings.Finite_T_Data.Liquid_H_MP - ...
            Settings.Finite_T_Data.Solid_H_MP;

        Settings.Finite_T_Data.Fusion_dV = Settings.Finite_T_Data.Liquid_V_MP - ...
            Settings.Finite_T_Data.Solid_V_MP;

        if isnan(Sol_Output.Solid_H_MP)
            Loss_add = nan;
            if ~Settings.Therm_Prop_Override
                Loss = nan;
                if MultiMetal || MultiHalide
                    UserData.Finite_T_Data.(Settings.Salt) = Settings.Finite_T_Data;
                    UserData.Minimization_Data.(Settings.Salt) = Settings.Minimization_Data;
                else
                    UserData.Finite_T_Data = Settings.Finite_T_Data;
                    UserData.Minimization_Data = Settings.Minimization_Data;
                end
                return
            end
        end
        if Settings.Parallel_Bayesopt
            setenv('OMP_NUM_THREADS',env.OMP_NUM_THREADS);
            setenv('GMX_PME_NUM_THREADS',env.GMX_PME_NUM_THREADS);
            setenv('GMX_PME_NTHREADS',env.GMX_PME_NTHREADS);
            setenv('GMX_OPENMP_MAX_THREADS',env.GMX_OPENMP_MAX_THREADS);
            setenv('KMP_AFFINITY',env.KMP_AFFINITY);
        end
    end

    % Delete previous calculations that did not complete
    if Settings.Therm_Prop_Override
        files = dir(Settings.scratch_dir);
        dirFlags = [files.isdir];
        subFolders = files(dirFlags);
        subFolderNames = {subFolders(3:end).name};
        prev_calcs = subFolderNames(cellfun(@(x) ~isempty(x),regexp(subFolderNames,'.+?_[S|M|L|O]P','once')));
        for idx = 1:length(prev_calcs)
            try
                rmdir(fullfile(Settings.scratch_dir,prev_calcs{idx}),'s')
            catch
                if Settings.Verbose
                    disp(['Unable to remove failed calculation directory: ' fullfile(Settings.scratch_dir,prev_calcs{idx})])
                end
            end
        end
    end
    Loss_Salt = Loss_Salt + LiX_Loss(Settings) + Loss_add;
    
    if MultiMetal || MultiHalide
        UserData.Finite_T_Data.(Settings.Salt) = Settings.Finite_T_Data;
        UserData.Minimization_Data.(Settings.Salt) = Settings.Minimization_Data;
    else
        UserData.Finite_T_Data = Settings.Finite_T_Data;
        UserData.Minimization_Data = Settings.Minimization_Data;
    end
end

% Remember: If skipping the finite_T data, then loss comes through the Loss_add
% variable, so we don't need to double count. Use Settings.skip_finite_T as
% a switch to disable re-calculating loss from finite T data
% Calculate Loss function from data
Loss = real(log1p(Loss_Salt));

end