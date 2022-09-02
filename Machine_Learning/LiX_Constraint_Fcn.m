function tf = LiX_Constraint_Fcn(Settings,Param)

Settings.Table_Length = 10; % nm
Settings.Table_StepSize = 0.01;

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

% Param will always be a table.
N_par = size(Param,1);
Loss = zeros(N_par,1);
% Potential Scaling
if strcmp(Settings.Theory,'TF')

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

elseif strcmp(Settings.Theory,'BH')

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

elseif strcmp(Settings.Theory,'JC') % JC models

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
end

% Calculate loss due to infeasible models with no well minima (only works reliably for BH/TF models in sigma-epsilon form)
if strcmp(Settings.Theory,'BH')
    U = BH_Potential_Generator_vec(Settings);
elseif strcmp(Settings.Theory,'TF')
    U = TF_Potential_Generator_vec(Settings);
elseif strcmp(Settings.Theory,'JC')
    U = JC_Potential_Generator_vec(Settings);
end

%% Grab the peaks and valleys of the MX attractive potential
peaks_idx = islocalmax(U.MX,2);
valleys_idx = islocalmin(U.MX,2);

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
ov_idx = sum( cumprod(valleys_idx(ovnp_idx,:) == 0, 2), 2) + 1;
ovl_idx = sub2ind(size(U_MX_ovnp_set),(1:numel(ov_idx)).',ov_idx);
U_min = U_MX_ovnp_set(ovl_idx);
Loss(ovnp_idx) = Loss(ovnp_idx) + max(Settings.MaxAttWellDepth - U_min,0).*Settings.BadFcnLossPenalty;

% If both a peak and minimum exist, find the well depth
U_MX_ovop_set = U.MX(ovop_idx,:);
op_idx = sum( cumprod(peaks_idx(ovop_idx,:) == 0, 2), 2)  + 1;
opl_idx = sub2ind(size(U_MX_ovop_set),(1:numel(op_idx)).',op_idx);
ov_idx = sum( cumprod(valleys_idx(ovop_idx,:) == 0, 2), 2) + 1;
ovl_idx = sub2ind(size(U_MX_ovop_set),(1:numel(ov_idx)).',ov_idx);
dU = U_MX_ovop_set(opl_idx) - U_MX_ovop_set(ovl_idx);
U_min = U_MX_ovop_set(ovl_idx);
% Penalize wells that are too shallow and wells that are too deep
Loss(ovop_idx) = Loss(ovop_idx) + max(Settings.MinExpWallHeight - dU,0).*Settings.BadFcnLossPenalty; % too shallow
Loss(ovop_idx) = Loss(ovop_idx) + max(Settings.MaxAttWellDepth - U_min,0).*Settings.BadFcnLossPenalty; % too deep

% idx = 9983;
% plot(U.r,U.MX(idx,:))
% hold on
% scatter(U.r(peaks_idx(idx,:)),U.MX(idx,peaks_idx(idx,:)))
% scatter(U.r(valleys_idx(idx,:)),U.MX(idx,valleys_idx(idx,:)))
% ylim([-1000 1000])

%% Grab the peaks and valleys of the MM/XX potentials
YY = {'MM' 'XX'};
for j = 1:2
    jj = YY{j};

    peaks_idx = islocalmax(U.(jj),2);
    valleys_idx = islocalmin(U.(jj),2);

    Num_peaks = sum(peaks_idx,2);
    Num_valleys = sum(valleys_idx,2);

    % 5 possible cases
    %npnv_idx = (Num_peaks == 0) & (Num_valleys == 0); % Potentials that contain no peaks and no valleys
    opnv_idx = (Num_peaks == 1) & (Num_valleys == 0); % Potentials that contain one peak and no valleys
    tpov_idx = (Num_peaks == 2);                      % Potentials that countain 2 peaks and one valley
    opov_idx = (Num_peaks == 1) & (Num_valleys == 1); % Potentials that contain one peak and one valley
    npov_idx = (Num_peaks == 0) & (Num_valleys == 1); % Potentials that contain no peaks and one valley

    % No peak and no valley exists: Do nothing, this is normal for JC and sometimes BH/TF
    %npnv_idx;

    % One peak, no valley: This is a normal case for TF and BH
    % repulsive potentials, check to ensure the peak height is at
    % least Settings.MinExpWallHeight
    U_jj_opnv = U.(jj)(opnv_idx,:)';
    dU = U_jj_opnv(peaks_idx(opnv_idx,:)');
    Loss(opnv_idx) = Loss(opnv_idx) + max(Settings.MinExpWallHeight - dU,0).*Settings.BadFcnLossPenalty;
    
    % two peaks are visible implies one valley in between them: 
    % Penalize any model with a non-zero well depth greater than Settings.MaxRepWellDepth
    U_jj_tpov_set = U.(jj)(tpov_idx,:);
    sp_idx = length(U.r) - sum( cumprod(fliplr(peaks_idx(tpov_idx,:)) == 0, 2), 2);
    spl_idx = sub2ind(size(U_jj_tpov_set),(1:numel(sp_idx)).',sp_idx);
    ov_idx = sum( cumprod(valleys_idx(tpov_idx,:) == 0, 2), 2) + 1;
    ovl_idx = sub2ind(size(U_jj_tpov_set),(1:numel(ov_idx)).',ov_idx);
    dU = U_jj_tpov_set(spl_idx) - U_jj_tpov_set(ovl_idx);
    Loss(tpov_idx) = Loss(tpov_idx) + max(dU - Settings.MaxRepWellDepth,0).*Settings.BadFcnLossPenalty;

    % One peak visible + one valley, and (possibly) one hidden peak to the right
    % In this case, check if the valley is closer-in than the peak
    % If the valley is closer-in, there is no hidden peak, use dU = peak_U - valley_U
    % If the peak is closer-in, there is a second hidden peak to the right, use dU = U(end) - valley_U
    opov_set = peaks_idx(opov_idx,:);
    U_jj_opov_set = U.(jj)(opov_idx,:);

    op_idx = sum(cumprod(opov_set == 0, 2), 2) + 1;
    ov_idx = sum(cumprod(valleys_idx(opov_idx,:) == 0, 2), 2) + 1;

    % split into two subsets
    ovcc_idx = (ov_idx < op_idx);
    Lossn = zeros(size(ovcc_idx));

    % Valley closer set: use peak_U - valley_U
    U_jj_opovc_set = U_jj_opov_set(ovcc_idx,:);
    opc_idx = op_idx(ovcc_idx,:);
    opcl_idx = sub2ind(size(U_jj_opovc_set),(1:numel(opc_idx)).',opc_idx);
    ovc_idx = ov_idx(ovcc_idx,:);
    ovcl_idx = sub2ind(size(U_jj_opovc_set),(1:numel(ovc_idx)).',ovc_idx);
    dUc = U_jj_opovc_set(opcl_idx) - U_jj_opovc_set(ovcl_idx);
    Lossn(ovcc_idx) =  max(dUc - Settings.MaxRepWellDepth,0).*Settings.BadFcnLossPenalty;

    % Valley further set: use U(end) - valley_U
    U_jj_opovf_set = U_jj_opov_set(~ovcc_idx,:);
    opf_idx = repmat(length(U.r),sum(~ovcc_idx),1); % set the "peak" to the furthest possible place
    opfl_idx = sub2ind(size(U_jj_opovf_set),(1:numel(opf_idx)).',opf_idx);
    ovf_idx = ov_idx(~ovcc_idx,:);
    ovfl_idx = sub2ind(size(U_jj_opovf_set),(1:numel(ovf_idx)).',ovf_idx);
    dUf = U_jj_opovf_set(opfl_idx) - U_jj_opovf_set(ovfl_idx);
    Lossn(~ovcc_idx) =  max(dUf - Settings.MaxRepWellDepth,0).*Settings.BadFcnLossPenalty;

    Loss(opov_idx) = Loss(opov_idx) + Lossn;

    % No peak + one valley: there must be a hidden peak to the right
    % use dU = U(end) - valley_U
    U_jj_npov_set = U.(jj)(npov_idx,:);
    ov_idx = sum(cumprod(valleys_idx(npov_idx,:) == 0, 2), 2) + 1;
    ovl_idx = sub2ind(size(U_jj_npov_set),(1:numel(ov_idx)).',ov_idx);
    op_idx = repmat(length(U.r),sum(npov_idx),1); % set the "peak" to the furthest possible place
    opl_idx = sub2ind(size(U_jj_npov_set),(1:numel(op_idx)).',op_idx);
    dU = U_jj_npov_set(opl_idx) - U_jj_npov_set(ovl_idx);
    Loss(npov_idx) = Loss(npov_idx) +  max(dU - Settings.MaxRepWellDepth,0).*Settings.BadFcnLossPenalty;

%             idxes = find(Loss > 0);
%             for jdx = 1:100
%                 idx=idxes(jdx);
%                 hold on
%                 plot(U.r,U.(jj)(idx,:))
%                 scatter(U.r(peaks_idx(idx,:)),U.(jj)(idx,peaks_idx(idx,:)))
%                 scatter(U.r(valleys_idx(idx,:)),U.(jj)(idx,valleys_idx(idx,:)))
%             end
%             ylim([-1000 1000])
end

tf = log1p(Loss) < sqrt(eps);

% % Plot result to visualize
% tf_num = double(tf);
% 
% % 'r0_MM'  'r0_XX'  'epsilon_MM'  'epsilon_XX'  'gamma_MX'
% % 'sigma_MM'  'sigma_XX'  'epsilon_MM'  'epsilon_XX'
% 
% ax1 = 'r0_MM';
% ax2 = 'epsilon_MM';
% ax3 = 'gamma_MX';
% 
% scatter3(Param.(ax1),Param.(ax2),Param.(ax3),50,tf_num,'filled')
% if strcmp(Settings.Theory,'TF')
%     set(gca, 'YScale', 'log')
%     set(gca, 'ZScale', 'log')
% end
% xlabel(ax1);
% ylabel(ax2);
% zlabel(ax3);
% 
% clear;
end