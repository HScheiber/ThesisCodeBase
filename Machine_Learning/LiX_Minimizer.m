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

% Initialize coupled constraints to -1, this indicates the constraints are satisfied by default
% The first coupled constraint is for model volumes at 0 K that are either
% larger than Settings.MaxModelVolume or smaller than Settings.MinModelVolume.
%
% The second coupled constraint is for finite temperature simulations that
% give NaN outputs, or for Structure_Minimization calculations that are not feasible. 
% These are caused when 
%   (1) The lattice params of a structure are too large or too small
%   (2) The lattice energy of a structure is lower than Settings.MinMDP.E_Unphys
%   (3) A finite-T simulation fails 
%   (4) A melting point cannot be found for the structure of interest 
%   (5) The liquid is amorphous at the experimental MP 
%   (6) The liquid or solid converts to another structure at the experimental MP
coupledconstraints = -ones(1,2);

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

% Maintain backwards compatibility, add defaults that may be missing to C6Damp
% C6D = Init_C6Damping_Object;
% C6D_n = fieldnames(C6D);
% for idx = 1:length(C6D_n)
%     if ~isfield(Model.C6Damp,C6D_n{idx})
%         Model.C6Damp.(C6D_n{idx}) = C6D.(C6D_n{idx});
%     end
% end

[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);

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
            def_E_MX = MXParams.epsilon;
            
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

Settings.Minimization_Data = Initialize_Minimization_Data(Settings);
Settings.Finite_T_Data = Initialize_Finite_T_Data(Settings);

% Calculate loss due to infeasible TF / BH models with no well minima (only works reliably in sigma-epsilon form)
Loss_add = 0;
if Settings.CheckBadFcn
    tl = Settings.Table_Length;
    ss = Settings.Table_StepSize;
    Settings.Table_Length = 10; % nm
    Settings.Table_StepSize = 0.002;
    if strcmp(Settings.Theory,'BH')
        [U_MX, U_MM, U_XX] = BH_Potential_Generator(Settings,...
            'Startpoint',0.01,'ReturnAsStructure',true);
    elseif strcmp(Settings.Theory,'TF')
        [U_MX, U_MM, U_XX] = TF_Potential_Generator(Settings,...
            'Startpoint',0.01,'ReturnAsStructure',true);
    elseif strcmp(Settings.Theory,'JC')
        [U_MX, U_MM, U_XX] = JC_Potential_Generator(Settings,...
            'Startpoint',0.01,'ReturnAsStructure',true);
    end
    Settings.Table_Length = tl; % nm
    Settings.Table_StepSize = ss;
    
    %% Grab the peaks and valleys of the MX attractive potential
    peaks_idx = islocalmax(U_MX.Total);
    valleys_idx = islocalmin(U_MX.Total);
    
    maxima_U = U_MX.Total(peaks_idx);
    minima_U = U_MX.Total(valleys_idx);
    
    if isempty(minima_U) % If no well minimum exists in MX interaction
        Loss_add = Loss_add + Settings.BadFcnLossPenalty;
    elseif isempty(maxima_U) % If no peak exists in MX interaction
        % Do nothing, this is normal for JC potential and some BH/TF
        % potentials
    else % Otherwise, a well minimum exists and at least one peak exists
        % ensure peak - well height is greater than specified threshold
        Threshold = Settings.MinExpWallHeight; % kJ/mol
        dU = maxima_U - minima_U;
        Loss_add = Loss_add + max(Threshold - dU,0)*Settings.BadFcnLossPenalty/Threshold;
    end
%     plot(U_MX.r,U_MX.Total)
%     hold on
%     scatter(U_MX.r(peaks_idx),U_MX.Total(peaks_idx))
%     scatter(U_MX.r(valleys_idx),U_MX.Total(valleys_idx))
%     ylim([-1000 1000])
    
    %% Grab the peaks and valleys of the MM/XX potentials
    for U = [U_MM,U_XX]
        peaks_idx = islocalmax(U.Total);
        valleys_idx = islocalmin(U.Total);
        
        maxima_U = U.Total(peaks_idx);
        minima_U = U.Total(valleys_idx);
        
        maxima_r = U.r(peaks_idx);
        minima_r = U.r(valleys_idx);

%         plot(U.r,U.Total)
%         hold on
%         scatter(U.r(peaks_idx),U.Total(peaks_idx))
%         scatter(U.r(valleys_idx),U.Total(valleys_idx))
%         ylim([-1000 1000])

        if isempty(maxima_U) % No peak exists
            % Do nothing, this is normal for JC and sometimes BH/TF
        elseif length(maxima_U) > 1 % two peaks are visible + one valley in between them
            % Penalize any model with a non-zero well depth between
            % like-like interactions
            Threshold = Settings.MaxRepWellDepth; % kJ/mol
            dU = maxima_U(2) - minima_U;
            Loss_add = Loss_add + max(dU - Threshold,0)*Settings.BadFcnLossPenalty;
        elseif length(maxima_U) == 1 && isempty(minima_U) % One peak, no valley (this is a normal case for TF and BH repulsive potentials)
            % Do nothing
        elseif length(maxima_U) == 1 && length(minima_U) == 1 % One peak visible + one valley in between, and (possibly) one hidden peak to the right
            Threshold = Settings.MaxRepWellDepth; % kJ/mol
            if minima_r < maxima_r
                dU = maxima_U - minima_U;
            else
                % Case of hidden peak to the right
                % Ensure valley depth is not greater than the threshold
                dU = U.Total(end) - minima_U;
            end
            Loss_add = Loss_add + max(dU - Threshold,0)*Settings.BadFcnLossPenalty;
        elseif length(minima_U) == 1 % well minima is available but no peaks are visible, there must be a hidden peak to the right
            Threshold = Settings.MaxRepWellDepth; % kJ/mol
            dU = U.Total(end) - minima_U;
            Loss_add = Loss_add + max(dU - Threshold,0)*Settings.BadFcnLossPenalty;
        else
            % This should never be reached...
            warning('Possible issue with the potential!')
        end
    end
end

if ~isfield(Settings,'MinSkipLoss')
    defSettings = Initialize_LiX_BO_Settings;
    Settings.MinSkipLoss = defSettings.MinSkipLoss;
end

% if Loss_add >= Settings.MinSkipLoss
%     Loss = Loss_add;
%     UserData.Minimization_Data = Settings.Minimization_Data;
%     UserData.Finite_T_Data = Settings.Finite_T_Data;
%     return
% end

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
    coupledconstraints(2) = 1;
    Loss = real(log1p(Loss_add + Settings.BadFcnLossPenalty));
    UserData.Minimization_Data = Settings.Minimization_Data;
    UserData.Finite_T_Data = Settings.Finite_T_Data;
    return
elseif Structure_Min_Calc_Fail && Settings.Therm_Prop_Override
    coupledconstraints(2) = 1;
    Loss_add = Loss_add + Settings.BadFcnLossPenalty;
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
        coupledconstraints(1) = Model_Mismatch;
        Loss_add_Vol = Model_Mismatch*Settings.BadFcnLossPenalty;
    elseif V0_model < Settings.MinModelVolume
        Model_Mismatch = (Settings.MinModelVolume - V0_model)/Settings.MinModelVolume;
        coupledconstraints(1) = Model_Mismatch;
        Loss_add_Vol = Model_Mismatch*Settings.BadFcnLossPenalty;
    else
        Loss_add_Vol = 0;
    end

    Loss_add = Loss_add + Loss_add_Vol;
    if Loss_add_Vol >= Settings.MinSkipLoss && ~Settings.Therm_Prop_Override
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
        coupledconstraints(2) = 1;
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
        dd = Settings.JobSettings.dd;
        npme = Settings.JobSettings.npme;
        Settings.JobSettings.dd = [];
        Settings.JobSettings.npme = [];
        [~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts] = MD_Batch_Template(Settings.JobSettings);
        Verbose = Settings.Verbose;
        Settings.Verbose = true;
        Liq_Output = Calc_Liquid_Properties_at_MP(Settings); % Output is nan if liquid converts to >0.85 solid
        Settings.Verbose = Verbose;
        Settings.JobSettings.dd = dd;
        Settings.JobSettings.npme = npme;
        [~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts] = MD_Batch_Template(Settings.JobSettings);
    end
    
    if Settings.Therm_Prop_Override && ~isfolder(ThermFolder)
        copyfile(Settings.WorkDir,ThermFolder)
    end
    
    if Settings.Delete_Equil && isfolder(Settings.WorkDir)
        try
            rmdir(Settings.WorkDir,'s')
        catch
            if Settings.Verbose
                warning(['Unable to remove directory: ' Settings.WorkDir])
            end
        end
    end
    
    if isnan(Liq_Output.Liquid_H_MP)
        coupledconstraints(2) = 1;
    end
    
    Settings.Finite_T_Data.Liquid_V_MP = Liq_Output.Liquid_V_MP;
    Settings.Finite_T_Data.Liquid_H_MP = Liq_Output.Liquid_H_MP;
    Settings.Finite_T_Data.Liquid_DM_MP = Liq_Output.Liquid_DM_MP; % cm^2 / s
    Settings.Finite_T_Data.Liquid_DX_MP = Liq_Output.Liquid_DX_MP;
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
        dd = Settings.JobSettings.dd;
        npme = Settings.JobSettings.npme;
        Settings.JobSettings.dd = [];
        Settings.JobSettings.npme = [];
        [~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts] = MD_Batch_Template(Settings.JobSettings);
        Verbose = Settings.Verbose;
        Settings.Verbose = true;
        Sol_Output = Calc_Solid_Properties_at_MP(Settings);
        Settings.Verbose = Verbose;
        Settings.JobSettings.dd = dd;
        Settings.JobSettings.npme = npme;
        [~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts] = MD_Batch_Template(Settings.JobSettings);
    end
    
    if Settings.Therm_Prop_Override && ~isfolder(ThermFolder)
        copyfile(Settings.WorkDir,ThermFolder)
    end
    
    if Settings.Delete_Equil && isfolder(Settings.WorkDir)
        try
            rmdir(Settings.WorkDir,'s')
        catch
            if Settings.Verbose
                warning(['Unable to remove directory: ' Settings.WorkDir])
            end
        end
    end
    
    if isnan(Sol_Output.Solid_H_MP)
        coupledconstraints(2) = 1;
    end
    
    Settings.Finite_T_Data.Solid_V_MP = Sol_Output.Solid_V_MP;
    Settings.Finite_T_Data.Solid_H_MP = Sol_Output.Solid_H_MP;
    
    Settings.Finite_T_Data.Fusion_dH = Settings.Finite_T_Data.Liquid_H_MP - ...
        Settings.Finite_T_Data.Solid_H_MP;
    
    Settings.Finite_T_Data.Fusion_dV = Settings.Finite_T_Data.Liquid_V_MP - ...
        Settings.Finite_T_Data.Solid_V_MP;
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
    files = dir(Settings.OuterDir);
    dirFlags = [files.isdir];
    subFolders = files(dirFlags);
    subFolderNames = {subFolders(3:end).name};
    prev_calcs = subFolderNames(cellfun(@(x) ~isempty(x),regexp(subFolderNames,'.+?_[S|M|L]P','once')));
    for idx = 1:length(prev_calcs)
        try
            rmdir(fullfile(Settings.OuterDir,prev_calcs{idx}),'s')
        catch
            if Settings.Verbose
                disp(['Unable to remove failed calculation directory: ' fullfile(Settings.OuterDir,prev_calcs{idx})])
            end
        end
    end
end

% Remember: If skipping the finite_T data, then loss comes through the Loss_add
% variable, so we don't need to double count. Use Settings.skip_finite_T as
% a switch to disable re-calculating loss from finite T data
% Calculate Loss function from data
Loss = real(log1p(LiX_Loss(Settings) + Loss_add));
UserData.Finite_T_Data = Settings.Finite_T_Data;
UserData.Minimization_Data = Settings.Minimization_Data;

end