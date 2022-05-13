function Loss = LiX_Loss(Loss_Options,Minimization_Data,Salt)
    tol = sqrt(eps);

    % Insure the correct structure order
    N = length(Minimization_Data);
    Structures = cell(1,N);
    for idx = 1:N
        Structures{idx} = Minimization_Data{idx}.Structure;
    end
    
    % Define the regularization function
    if isa(Loss_Options.regularization,'function_handle')
        reg = Loss_Options.regularization;
    else
        switch Loss_Options.regularization
            case 'L1'
                reg = @(x) abs(x);
            case 'L2'
                reg = @(x) (x).^2;
            otherwise
                error(['Unknown regularization function: ' Loss_Options.regularization]);
        end
    end

    % Load DFT/experimental data
    DFT = Load_Best_DFT_Data;
    if Loss_Options.Experimental_LE || Loss_Options.Experimental_LP
        Experimental = Load_Experimental_Data;
        
        if Loss_Options.Experimental_LE
            E_Correction = Experimental.(Salt).Rocksalt.E - DFT.(Salt).Rocksalt.Energy;        
            for idx = 1:N
                try
                    DFT.(Salt).(Structures{idx}).Energy = DFT.(Salt).(Structures{idx}).Energy + E_Correction;
                catch
                    DFT.(Salt).(Structures{idx}).Energy = nan;
                    DFT.(Salt).(Structures{idx}).a = nan;
                    DFT.(Salt).(Structures{idx}).b = nan;
                    DFT.(Salt).(Structures{idx}).c = nan;
                end
            end
        end
        if Loss_Options.Experimental_LP
            RS_a_Correction = Experimental.(Salt).Rocksalt.a_zero - DFT.(Salt).Rocksalt.a;
            RS_V_correction = Experimental.(Salt).Rocksalt.V_zero - DFT.(Salt).Rocksalt.V;
            WZ_a_Correction = Experimental.(Salt).Wurtzite.a_zero - DFT.(Salt).Wurtzite.a;
            WZ_c_Correction = Experimental.(Salt).Wurtzite.c_zero - DFT.(Salt).Wurtzite.c;
            WZ_V_correction = Experimental.(Salt).Wurtzite.V_zero - DFT.(Salt).Wurtzite.V;
            
            DFT.(Salt).Rocksalt.a = DFT.(Salt).Rocksalt.a + RS_a_Correction;
            DFT.(Salt).Rocksalt.V = DFT.(Salt).Rocksalt.V + RS_V_correction;
            if ~isnan(WZ_a_Correction)
                DFT.(Salt).Wurtzite.a = DFT.(Salt).Wurtzite.a + WZ_a_Correction;
                DFT.(Salt).Wurtzite.c = DFT.(Salt).Wurtzite.c + WZ_c_Correction;
                DFT.(Salt).Wurtzite.V = DFT.(Salt).Wurtzite.V + WZ_V_correction;
            end
            
        end
    end
    
    %% Calculate Loss function
    try
        Min_RS_E = Minimization_Data{strcmpi(Structures,'Rocksalt')}.E; % Model-Calculated RS lattice energy
    catch
        Min_RS_E = nan;
    end
    DFT_RS_E = DFT.(Salt).Rocksalt.Energy;                          % DFT/Experimental-calculated   RS lattice energy
    
    Loss = 0; % Initialize loss
    for idx = 1:N
        Structure = Minimization_Data{idx}.Structure;
        Min_dat = Minimization_Data{idx};

        if ~isfield(Loss_Options,Structure)
            % Skip structures that are not contained in the loss options
            continue
        end
        
        %% Incorporate the absolute error in lattice energy for each structure.
        if Loss_Options.(Structure).LE < tol
            % Skip this if the weight on LE is below tolerance
        elseif ~isfield(Loss_Options.(Structure),'LE_freedom') || ...
                Loss_Options.(Structure).LE_freedom < abs( Min_dat.E - DFT.(Salt).(Structure).Energy )
            
            Loss = Loss + Loss_Options.(Structure).LE*...
                reg((Min_dat.E - DFT.(Salt).(Structure).Energy)/DFT.(Salt).(Structure).Energy);
        end

        %% Incorporate the relative error in lattice energy for each structure.
        % Grab the relative lattice energies
        DFT_diff_E = DFT.(Salt).(Structure).Energy - DFT_RS_E; 
        Min_diff_E = Min_dat.E                     - Min_RS_E;
        
        if Loss_Options.(Structure).RLE < tol
            % Skip this if the weight on RLE is below tolerance
        elseif ~isfield(Loss_Options.(Structure),'RLE_freedom') || ...
                Loss_Options.(Structure).RLE_freedom < abs(DFT_diff_E - Min_diff_E)
            
            Loss = Loss + Loss_Options.(Structure).RLE*...
                reg( (Min_diff_E - DFT_diff_E)/DFT_RS_E );
        end

        %% Incorporate the relative absolute error in lattice parameter / volume for each structure.
        if Loss_Options.(Structure).a > tol
            Loss = Loss + Loss_Options.(Structure).a*...
                reg( (Min_dat.a - DFT.(Salt).(Structure).a)/DFT.(Salt).(Structure).a );
        end
        
        if Loss_Options.(Structure).b > tol
            Loss = Loss + Loss_Options.(Structure).b*...
                reg( (Min_dat.b - DFT.(Salt).(Structure).b)/DFT.(Salt).(Structure).b );
        end
        
        if Loss_Options.(Structure).c > tol
            Loss = Loss + Loss_Options.(Structure).c*...
                reg( (Min_dat.c - DFT.(Salt).(Structure).c)/DFT.(Salt).(Structure).c );
        end
        
        if isfield(Loss_Options.(Structure),'V') && Loss_Options.(Structure).V > tol
            Loss = Loss + Loss_Options.(Structure).V*...
                reg( (Min_dat.V - DFT.(Salt).(Structure).V)/DFT.(Salt).(Structure).V );
        end

        %% Incorporate energy gaps for each structure
        % Comparison function: less than, greater than, equal to, etc.
        % Less than:    the Actual_Gap must be less than    the Target_Gap or the loss is non-zero.
        % Greater than: the Actual_Gap must be greater than the Target_Gap or the loss is non-zero.
        % Equal to:     the Actual_Gap must be equal to     the Target_Gap or the loss is non-zero.
        compare_fun = Loss_Options.(Structure).Gap.Type;
        
        % The reference structure for the gap. 
        % Current structure is looped through.
        Gap_ref_structure = Loss_Options.(Structure).Gap.Ref;
        
        % The target gap between reference and current structure
        % Gap.Value < 0 -> "reference structure" is favoured
        % Gap.Value > 0 -> "current structure" is favoured
        % Gap.Value = 0 -> structures are equal in energy
        Target_Gap = Loss_Options.(Structure).Gap.Value;
        
        % Difference in energy between the reference structure and the current structure. 
        % Negative values means the "reference structure" is favoured. 
        % Positive values means the "current structure" is favoured.
        try
            Actual_gap = Minimization_Data{strcmpi(Structures,Gap_ref_structure)}.E - Min_dat.E;
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
        Loss = Loss + Loss_Options.(Structure).Gap.Weight*reg(Delta_Gap);

    end
    
    Loss = real(log1p(Loss));
end