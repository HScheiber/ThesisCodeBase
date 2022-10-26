function Loss = Loss_Recalculate(Settings,Param,UserData,Structures)
UserData = [UserData{:}];
N_Dat = size(UserData,2);
N_Structures = numel(Structures);
Loss = zeros(N_Dat,1);

% Calculate loss due to infeasible models
if Settings.CheckBadFcn
    Loss_add = LiX_Constraint_Fcn(Settings,Param,'Output_Loss',true);
end

% Place minimization data into cell array
MinDat = [UserData.Minimization_Data]';
tol = sqrt(eps);

% Any finite T data
if any([Settings.Loss_Options.Fusion_Enthalpy ...
        Settings.Loss_Options.MP_Volume_Change ...
        Settings.Loss_Options.Liquid_MP_Volume ...
        Settings.Loss_Options.Solid_MP_Volume ...
        Settings.Loss_Options.MP] > tol)
    
    % Penalize for crystal size
    RefStruct = UserData(1).Finite_T_Data.Structure;
    V0_model = struct2table([MinDat{:,strcmp(RefStruct,Structures)}],'AsArray',true).V;
    Model_Mismatch = abs(V0_model - Settings.MaxModelVolume)./Settings.MaxModelVolume;
    Loss_add_Vol = zeros(N_Dat,1);
    Loss_add_Vol(V0_model > Settings.MaxModelVolume) = Model_Mismatch(V0_model > Settings.MaxModelVolume).*Settings.BadFcnLossPenalty;
    Loss_add_Vol(V0_model < Settings.MinModelVolume) = Model_Mismatch(V0_model < Settings.MinModelVolume).*Settings.BadFcnLossPenalty;
    Loss = Loss + Loss_add_Vol;
end

% Check min skip loss
skip_finite_T = real(log1p(Loss + Loss_add)) >= Settings.MinSkipLoss;

% Melting point
FinTDat = struct2table([UserData.Finite_T_Data],'AsArray',true);
if Settings.Loss_Options.MP > tol
    Loss(isnan(FinTDat.MP) & ~skip_finite_T) = nan;
end

% High T Liquid properties
if any([Settings.Loss_Options.Fusion_Enthalpy ...
        Settings.Loss_Options.MP_Volume_Change ...
        Settings.Loss_Options.Liquid_MP_Volume ...
        Settings.Loss_Options.Liquid_DM_MP] > tol)
    Loss(isnan(FinTDat.Liquid_H_MP) & ~skip_finite_T) = nan;
end

% High T solid properties
if any([Settings.Loss_Options.Fusion_Enthalpy ...
        Settings.Loss_Options.MP_Volume_Change ...
        Settings.Loss_Options.Solid_MP_Volume] > tol)
    Loss(isnan(FinTDat.Solid_H_MP) & ~skip_finite_T) = nan;
end

%% Calc loss only on remaining data
Unfinished_Loss_idx = ~isnan(Loss);
N_Remain = sum(Unfinished_Loss_idx);
MinDat_remain = MinDat(Unfinished_Loss_idx,:);
FinTDat_remain = FinTDat(Unfinished_Loss_idx,:);
Loss_remain = zeros(N_Remain,1);
do_finite_T_remain = ~skip_finite_T(Unfinished_Loss_idx);

if N_Remain > 0
    %% Define the regularization function
    if isa(Settings.Loss_Options.regularization,'function_handle')
        reg = Settings.Loss_Options.regularization;
    else
        switch Settings.Loss_Options.regularization
            case 'L1'
                reg = @(x) abs(x);
            case 'L2'
                reg = @(x) (x).^2;
            otherwise
                error(['Unknown regularization function: ' Settings.Loss_Options.regularization]);
        end
    end

    %% Load DFT/experimental data
    DFT = Load_Best_DFT_Data;
    if Settings.Loss_Options.Experimental_LE || Settings.Loss_Options.Experimental_LP
        Experimental = Load_Experimental_Data;

        if Settings.Loss_Options.Experimental_LE
            E_Correction = Experimental.(Settings.Salt).Rocksalt.E - DFT.(Settings.Salt).Rocksalt.Energy;        
            for idx = 1:N_Structures
                try
                    DFT.(Settings.Salt).(Structures{idx}).Energy = DFT.(Settings.Salt).(Structures{idx}).Energy + E_Correction;
                catch
                    DFT.(Settings.Salt).(Structures{idx}).Energy = nan;
                    DFT.(Settings.Salt).(Structures{idx}).a = nan;
                    DFT.(Settings.Salt).(Structures{idx}).b = nan;
                    DFT.(Settings.Salt).(Structures{idx}).c = nan;
                end
            end
        end
        if Settings.Loss_Options.Experimental_LP
            RS_a_Correction = Experimental.(Settings.Salt).Rocksalt.a_zero - DFT.(Settings.Salt).Rocksalt.a;
            RS_V_correction = Experimental.(Settings.Salt).Rocksalt.V_zero - DFT.(Settings.Salt).Rocksalt.V;
            DFT.(Settings.Salt).Rocksalt.a = DFT.(Settings.Salt).Rocksalt.a + RS_a_Correction;
            DFT.(Settings.Salt).Rocksalt.V = DFT.(Settings.Salt).Rocksalt.V + RS_V_correction;
            if isfield(Experimental.(Settings.Salt),'Wurtzite')
                WZ_a_Correction = Experimental.(Settings.Salt).Wurtzite.a_zero - DFT.(Settings.Salt).Wurtzite.a;
                WZ_c_Correction = Experimental.(Settings.Salt).Wurtzite.c_zero - DFT.(Settings.Salt).Wurtzite.c;
                WZ_V_correction = Experimental.(Settings.Salt).Wurtzite.V_zero - DFT.(Settings.Salt).Wurtzite.V;
                DFT.(Settings.Salt).Wurtzite.a = DFT.(Settings.Salt).Wurtzite.a + WZ_a_Correction;
                DFT.(Settings.Salt).Wurtzite.c = DFT.(Settings.Salt).Wurtzite.c + WZ_c_Correction;
                DFT.(Settings.Salt).Wurtzite.V = DFT.(Settings.Salt).Wurtzite.V + WZ_V_correction;
            end
        end
    end

    %% Calculate Losses for remaining points
    try
        Min_RS_E = struct2table([MinDat_remain{:,strcmp('Rocksalt',Structures)}],'AsArray',true).E; % Model-Calculated RS lattice energy
    catch
        Min_RS_E = nan(N_Remain,1);
    end
    DFT_RS_E = DFT.(Settings.Salt).Rocksalt.Energy; % DFT/Experimental-calculated   RS lattice energy

    % Solid State T = 0 properties
    for idx = 1:N_Structures
        Structure = Structures{idx};
        Min_dat = struct2table([MinDat_remain{:,strcmp(Structure,Structures)}],'AsArray',true);

        if ~isfield(Settings.Loss_Options,Structure)
            % Skip structures that are not contained in the loss options
            continue
        end
        
        %% Incorporate the absolute error in lattice energy for each structure.
        if Settings.Loss_Options.(Structure).LE >= tol
            fr_idx = Settings.Loss_Options.(Structure).LE_freedom < abs( Min_dat.E - DFT.(Settings.Salt).(Structure).Energy );

            Loss_remain(fr_idx) = Loss_remain(fr_idx) + Settings.Loss_Options.(Structure).LE.*...
                reg((Min_dat.E - DFT.(Settings.Salt).(Structure).Energy)./DFT.(Settings.Salt).(Structure).Energy);
        end
        
        %% Incorporate the relative error in lattice energy for each structure.
        % Grab the relative lattice energies
        DFT_diff_E = DFT.(Settings.Salt).(Structure).Energy - DFT_RS_E; 
        Min_diff_E = Min_dat.E - Min_RS_E;

        if Settings.Loss_Options.(Structure).RLE >= tol
            fr_idx = Settings.Loss_Options.(Structure).RLE_freedom < abs(DFT_diff_E - Min_diff_E);

            Loss_remain(fr_idx) = Loss_remain(fr_idx) + Settings.Loss_Options.(Structure).RLE.*...
                reg( (Min_diff_E - DFT_diff_E)./DFT_RS_E );
        end
        
        %% Incorporate the relative absolute error in lattice parameter / volume for each structure.
        if Settings.Loss_Options.(Structure).a > tol
            Loss_remain = Loss_remain + Settings.Loss_Options.(Structure).a.*...
                reg( (Min_dat.a - DFT.(Settings.Salt).(Structure).a)./DFT.(Settings.Salt).(Structure).a );
        end
        
        if Settings.Loss_Options.(Structure).b > tol
            Loss_remain = Loss_remain + Settings.Loss_Options.(Structure).b.*...
                reg( (Min_dat.b - DFT.(Settings.Salt).(Structure).b)./DFT.(Settings.Salt).(Structure).b );
        end
        
        if Settings.Loss_Options.(Structure).c > tol
            Loss_remain = Loss_remain + Settings.Loss_Options.(Structure).c.*...
                reg( (Min_dat.c - DFT.(Settings.Salt).(Structure).c)./DFT.(Settings.Salt).(Structure).c );
        end
        
        if isfield(Settings.Loss_Options.(Structure),'V') && Settings.Loss_Options.(Structure).V > tol
            Loss_remain = Loss_remain + Settings.Loss_Options.(Structure).V.*...
                reg( (Min_dat.V - DFT.(Settings.Salt).(Structure).V)./DFT.(Settings.Salt).(Structure).V );
        end
        
        %% Incorporate energy gaps for each structure
        % Comparison function: less than, greater than, equal to, etc.
        % Less than:    the Actual_Gap must be less than    the Target_Gap or the loss is non-zero.
        % Greater than: the Actual_Gap must be greater than the Target_Gap or the loss is non-zero.
        % Equal to:     the Actual_Gap must be equal to     the Target_Gap or the loss is non-zero.
        compare_fun = Settings.Loss_Options.(Structure).Gap.Type;
        
        % The reference structure for the gap. 
        % Current structure is looped through.
        Gap_ref_structure = Settings.Loss_Options.(Structure).Gap.Ref;
        
        % The target gap between reference and current structure
        % Gap.Value < 0 -> "reference structure" is favoured
        % Gap.Value > 0 -> "current structure" is favoured
        % Gap.Value = 0 -> structures are equal in energy
        Target_Gap = Settings.Loss_Options.(Structure).Gap.Value;
        
        % Difference in energy between the reference structure and the current structure. 
        % Negative values means the "reference structure" is favoured. 
        % Positive values means the "current structure" is favoured.
        try
            Actual_gap = struct2table([MinDat_remain{:,strcmp(Structures,Gap_ref_structure)}],'AsArray',true).E - Min_dat.E;
        catch
            Actual_gap = nan(N_Remain,1);
        end
        
        % Is the target satified?
        Delta_Gap = zeros(N_Remain,1);
        Gap_sat_idx = compare_fun(Actual_gap,Target_Gap);
        Delta_Gap(~Gap_sat_idx) = Target_Gap - Actual_gap(~Gap_sat_idx);
        
        % Add to the total loss function
        Loss_remain = Loss_remain + Settings.Loss_Options.(Structure).Gap.Weight.*reg(Delta_Gap);
    end
    
    % Calculate loss for finite temperature properties
    % Fit the Fusion enthalpy
    if Settings.Loss_Options.Fusion_Enthalpy > tol

        rel_er = Settings.Loss_Options.Fusion_Enthalpy.*...
            reg((FinTDat_remain.Fusion_dH - FinTDat_remain.Exp_Fusion_dH)./...
            FinTDat_remain.Exp_Fusion_dH);
        rel_er(isnan(rel_er)) = Settings.Loss_Options.Fusion_Enthalpy*Settings.BadFcnLossPenalty;
        Loss_remain(do_finite_T_remain) = Loss_remain(do_finite_T_remain) + rel_er(do_finite_T_remain);
    end

    % Fit the change in volume due to melting
    if Settings.Loss_Options.MP_Volume_Change > tol

        rel_er = Settings.Loss_Options.MP_Volume_Change.*...
            reg((FinTDat_remain.Fusion_dV - FinTDat_remain.Exp_Fusion_dV)./...
            FinTDat_remain.Exp_Fusion_dV);

        rel_er(isnan(rel_er)) = Settings.Loss_Options.MP_Volume_Change*Settings.BadFcnLossPenalty;
        Loss_remain(do_finite_T_remain) = Loss_remain(do_finite_T_remain) + rel_er(do_finite_T_remain);
    end

    % Fit the liquid molar volume at the MP
    if Settings.Loss_Options.Liquid_MP_Volume > tol

        rel_er = Settings.Loss_Options.Liquid_MP_Volume.*...
            reg((FinTDat_remain.Liquid_V_MP - FinTDat_remain.Exp_Liquid_V_MP)./...
            FinTDat_remain.Exp_Liquid_V_MP);

        rel_er(isnan(rel_er)) = Settings.Loss_Options.Liquid_MP_Volume*Settings.BadFcnLossPenalty;
        Loss_remain(do_finite_T_remain) = Loss_remain(do_finite_T_remain) + rel_er(do_finite_T_remain);
    end

    % Fit the solid molar volume at the MP
    if Settings.Loss_Options.Solid_MP_Volume > tol

        rel_er = Settings.Loss_Options.Solid_MP_Volume.*...
            reg((FinTDat_remain.Solid_V_MP - FinTDat_remain.Exp_Solid_V_MP)./...
            FinTDat_remain.Exp_Solid_V_MP);

        rel_er(isnan(rel_er)) = Settings.Loss_Options.Solid_MP_Volume*Settings.BadFcnLossPenalty;
        Loss_remain(do_finite_T_remain) = Loss_remain(do_finite_T_remain) + rel_er(do_finite_T_remain);
    end

    % Fit the metal ion diffusivity at the MP
    if Settings.Loss_Options.Liquid_DM_MP > tol

        rel_er = Settings.Loss_Options.Liquid_DM_MP.*...
            reg((FinTDat_remain.Liquid_DM_MP - FinTDat_remain.Exp_DM_MP)./...
            FinTDat_remain.Exp_DM_MP);

        rel_er(isnan(rel_er)) = Settings.Loss_Options.Liquid_DM_MP*Settings.BadFcnLossPenalty;
        Loss_remain(do_finite_T_remain) = Loss_remain(do_finite_T_remain) + rel_er(do_finite_T_remain);
    end

    % Fit the experimental melting point
    if Settings.Loss_Options.MP > tol

        rel_er = Settings.Loss_Options.MP.*...
            reg((FinTDat_remain.MP - FinTDat_remain.Exp_MP)./...
            FinTDat_remain.Exp_MP);

        rel_er(isnan(rel_er)) = Settings.Loss_Options.MP*Settings.BadFcnLossPenalty;
        Loss_remain(do_finite_T_remain) = Loss_remain(do_finite_T_remain) + rel_er(do_finite_T_remain);
    end
else
    Loss_remain = [];
end

% Add up loss
Loss(Unfinished_Loss_idx) = Loss(Unfinished_Loss_idx) + Loss_remain;
Loss = real(log1p(Loss + Loss_add));

% This catches Structure_Minimization calculations that produced a result
% outside of the allowed energy or volume bounds
Failed_Calcs = any(cellfun(@(x) x.CalcFail, MinDat),2);
if Settings.UseCoupledConstraint
    Loss(Failed_Calcs) = nan;
else
    Loss(Failed_Calcs) = real(log1p(Loss_add(Failed_Calcs) + Settings.BadFcnLossPenalty));
end

end