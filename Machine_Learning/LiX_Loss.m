function Loss = LiX_Loss(Settings)
    tol = sqrt(eps);

    % Insure the correct structure order
    N = length(Settings.Minimization_Data);
    Structures = cell(1,N);
    for idx = 1:N
        Structures{idx} = Settings.Minimization_Data{idx}.Structure;
    end
    
    % Define the regularization function
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

    % Load DFT/experimental data
    DFT = Load_Best_DFT_Data;
    if Settings.Loss_Options.Experimental_LE || Settings.Loss_Options.Experimental_LP
        Experimental = Load_Experimental_Data;
        
        if Settings.Loss_Options.Experimental_LE
            E_Correction = Experimental.(Settings.Salt).Rocksalt.E - DFT.(Settings.Salt).Rocksalt.Energy;        
            for idx = 1:N
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
    
    %% Calculate Loss function
    try
        Min_RS_E = Settings.Minimization_Data{strcmpi(Structures,'Rocksalt')}.E; % Model-Calculated RS lattice energy
    catch
        Min_RS_E = nan;
    end
    DFT_RS_E = DFT.(Settings.Salt).Rocksalt.Energy;                          % DFT/Experimental-calculated   RS lattice energy
    
    % Solid State T = 0 properties
    Loss = 0; % Initialize loss
    for idx = 1:N
        Structure = Settings.Minimization_Data{idx}.Structure;
        Min_dat = Settings.Minimization_Data{idx};

        if ~isfield(Settings.Loss_Options,Structure)
            % Skip structures that are not contained in the loss options
            continue
        end
        
        %% Incorporate the absolute error in lattice energy for each structure.
        if Settings.Loss_Options.(Structure).LE < tol
            % Skip this if the weight on LE is below tolerance
        elseif ~isfield(Settings.Loss_Options.(Structure),'LE_freedom') || ...
                Settings.Loss_Options.(Structure).LE_freedom < abs( Min_dat.E - DFT.(Settings.Salt).(Structure).Energy )
            
            Loss = Loss + Settings.Loss_Options.(Structure).LE*...
                reg((Min_dat.E - DFT.(Settings.Salt).(Structure).Energy)/DFT.(Settings.Salt).(Structure).Energy);
        end

        %% Incorporate the relative error in lattice energy for each structure.
        % Grab the relative lattice energies
        DFT_diff_E = DFT.(Settings.Salt).(Structure).Energy - DFT_RS_E; 
        Min_diff_E = Min_dat.E                     - Min_RS_E;
        
        if Settings.Loss_Options.(Structure).RLE < tol
            % Skip this if the weight on RLE is below tolerance
        elseif ~isfield(Settings.Loss_Options.(Structure),'RLE_freedom') || ...
                Settings.Loss_Options.(Structure).RLE_freedom < abs(DFT_diff_E - Min_diff_E)
            
            Loss = Loss + Settings.Loss_Options.(Structure).RLE*...
                reg( (Min_diff_E - DFT_diff_E)/DFT_RS_E );
        end

        %% Incorporate the relative absolute error in lattice parameter / volume for each structure.
        if Settings.Loss_Options.(Structure).a > tol
            Loss = Loss + Settings.Loss_Options.(Structure).a*...
                reg( (Min_dat.a - DFT.(Settings.Salt).(Structure).a)/DFT.(Settings.Salt).(Structure).a );
        end
        
        if Settings.Loss_Options.(Structure).b > tol
            Loss = Loss + Settings.Loss_Options.(Structure).b*...
                reg( (Min_dat.b - DFT.(Settings.Salt).(Structure).b)/DFT.(Settings.Salt).(Structure).b );
        end
        
        if Settings.Loss_Options.(Structure).c > tol
            Loss = Loss + Settings.Loss_Options.(Structure).c*...
                reg( (Min_dat.c - DFT.(Settings.Salt).(Structure).c)/DFT.(Settings.Salt).(Structure).c );
        end
        
        if isfield(Settings.Loss_Options.(Structure),'V') && Settings.Loss_Options.(Structure).V > tol
            Loss = Loss + Settings.Loss_Options.(Structure).V*...
                reg( (Min_dat.V - DFT.(Settings.Salt).(Structure).V)/DFT.(Settings.Salt).(Structure).V );
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
            Actual_gap = Settings.Minimization_Data{strcmpi(Structures,Gap_ref_structure)}.E - Min_dat.E;
        catch
            Actual_gap = nan;
        end
        
        % Is the target satified?        
        if compare_fun(Actual_gap,Target_Gap)
            Delta_Gap = 0;
        else
            Delta_Gap = Target_Gap - Actual_gap;
        end
        
        % Add to the total loss function
        Loss = Loss + Settings.Loss_Options.(Structure).Gap.Weight*reg(Delta_Gap);
    end
    
    % Finite Temperature properties
    %     Finite_T_Data.Fusion_dH = nan;
    %     Finite_T_Data.Fusion_dV = nan;
    %     Finite_T_Data.Liquid_V_MP = nan;
    %     Finite_T_Data.Solid_V_MP = nan;
    %     Finite_T_Data.MP = nan;
    
    %     Finite_T_Data.Exp_Fusion_dH % kJ/mol of Formula units
    %     Finite_T_Data.Exp_Fusion_dV % Volume change on melting at MP: Ang^3 / Forumla unit
    %     Finite_T_Data.Exp_Liquid_V_MP % Liquid Volume at MP in Ang^3 / Forumla unit
    %     Finite_T_Data.Exp_Solid_V_MP % Solid Volume at MP in Ang^3 / Forumla unit
    %     Finite_T_Data.Exp_MP % MP in K
    
    % Calculate loss for finite temperature properties
    if ~Settings.skip_finite_T
        
        % Fit the Fusion enthalpy
        if Settings.Loss_Options.Fusion_Enthalpy > tol
            rel_er = Settings.Loss_Options.Fusion_Enthalpy*...
                reg((Settings.Finite_T_Data.Fusion_dH - Settings.Finite_T_Data.Exp_Fusion_dH)/...
                Settings.Finite_T_Data.Exp_Fusion_dH);
            if isnan(rel_er)
                rel_er = Settings.Loss_Options.Fusion_Enthalpy*Settings.BadFcnLossPenalty;
            end
            Loss = Loss + rel_er;
        end
        
        % Fit the change in volume due to melting
        if Settings.Loss_Options.MP_Volume_Change > tol
            rel_er = Settings.Loss_Options.MP_Volume_Change*...
                reg((Settings.Finite_T_Data.Fusion_dV - Settings.Finite_T_Data.Exp_Fusion_dV)/...
                Settings.Finite_T_Data.Exp_Fusion_dV);
            if isnan(rel_er)
                rel_er = Settings.Loss_Options.MP_Volume_Change*Settings.BadFcnLossPenalty;
            end
            Loss = Loss + rel_er;
        end
        
        % Fit the liquid molar volume at the MP
        if Settings.Loss_Options.Liquid_MP_Volume > tol
            rel_er = Settings.Loss_Options.Liquid_MP_Volume*...
                reg((Settings.Finite_T_Data.Liquid_V_MP - Settings.Finite_T_Data.Exp_Liquid_V_MP)/...
                Settings.Finite_T_Data.Exp_Liquid_V_MP);
            if isnan(rel_er)
                rel_er = Settings.Loss_Options.Liquid_MP_Volume*Settings.BadFcnLossPenalty;
            end
            Loss = Loss + rel_er;
        end
        
        % Fit the solid molar volume at the MP
        if Settings.Loss_Options.Solid_MP_Volume > tol
            rel_er = Settings.Loss_Options.Solid_MP_Volume*...
                reg((Settings.Finite_T_Data.Solid_V_MP - Settings.Finite_T_Data.Exp_Solid_V_MP)/...
                Settings.Finite_T_Data.Exp_Solid_V_MP);
            if isnan(rel_er)
                rel_er = Settings.Loss_Options.Solid_MP_Volume*Settings.BadFcnLossPenalty;
            end
            Loss = Loss + rel_er;
        end
        
        
        % Fit the metal ion diffusivity at the MP
        if Settings.Loss_Options.Liquid_DM_MP > tol
            rel_er = Settings.Loss_Options.Liquid_DM_MP*...
                reg((Settings.Finite_T_Data.Liquid_DM_MP - Settings.Finite_T_Data.Exp_DM_MP)/...
                Settings.Finite_T_Data.Exp_DM_MP);
            if isnan(rel_er)
                rel_er = Settings.Loss_Options.Liquid_DM_MP*Settings.BadFcnLossPenalty;
            end
            Loss = Loss + rel_er;
        end
        
        
        % Fit the experimental melting point
        if Settings.Loss_Options.MP > tol
            rel_er = Settings.Loss_Options.MP*...
                reg((Settings.Finite_T_Data.MP - Settings.Finite_T_Data.Exp_MP)/...
                Settings.Finite_T_Data.Exp_MP);
            if isnan(rel_er)
                rel_er = Settings.Loss_Options.MP*Settings.BadFcnLossPenalty;
            end
            Loss = Loss + rel_er;
        end
    end
end