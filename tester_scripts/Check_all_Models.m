
diary('Check_models.log')
model_loc = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';
fullopts = dir(fullfile(model_loc,'*fullopt.mat'));

for idx = 1:length(fullopts)
    %bayeopt_file = strrep(fullopts(idx).name,'fullopt','bayesopt');
    %OptParamdat = load(fullfile(model_loc,fullopts(idx).name));
    outp = regexp(fullopts(idx).name,'([a-zA-Z]{2,4})_([a-zA-Z]{2,4})_Model_([a-zA-Z0-9]{1,5})','tokens','once');
    
    Settings.Salt = outp{1};
    Settings.Theory = outp{2};
    Settings.Model = outp{3};
    [Settings,ModelFound] = Load_Model_Params(Settings);
    
    if ~ModelFound
        continue
    end
    
    Settings.Table_Length = 10; % nm
    Settings.Table_StepSize = 0.002;
    Settings.MaxRepWellDepth = 0; % kJ/mol
    Settings.MinExpWallHeight = 100; % kJ/mol
    Settings.BadFcnLossPenalty = 1000;
    Settings.MDP.RVDW_Cutoff = 2;
    Settings.MDP.vdw_modifier = 'potential-shift';
    
    if ~isfield(Settings,'CR_Damp')
        Settings.CR_Damp = Init_CRDamping_Object;
    end
    if ~isfield(Settings,'C6_Damp')
        Settings.C6_Damp = Init_C6Damping_Object;
    end
    if ~isfield(Settings,'TF_Paramset')
        Settings.TF_Paramset = 0;
    end
    
    switch outp{2}
        % Check for JC and BH parameters that are too small
        case 'JC'
            Settings.WaterModel = 'SPC/E';
            [UP.MX,UP.MM,UP.XX] = JC_Potential_Parameters(Settings);
            [U_MX, U_MM, U_XX] = JC_Potential_Generator(Settings,...
                'Startpoint',0.01,'ReturnAsStructure',true);
            
            for ID = ["MX" "MM" "XX"]
                if UP.(ID).sigma < 0 || abs(UP.(ID).sigma) < 1e-5
                    disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - Parameter ' char(ID) '.sigma = ' num2str(UP.(ID).sigma,'%10.8e')])
                end
                if UP.(ID).epsilon < 0 || abs(UP.(ID).epsilon) < 1e-5
                    disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - Parameter ' char(ID) '.epsilon = ' num2str(UP.(ID).epsilon,'%10.8e')])
                end
            end
            
        case 'BH'
            [UP.MX,UP.MM,UP.XX] = BH_Potential_Parameters(Settings);
            [U_MX, U_MM, U_XX] = BH_Potential_Generator(Settings,...
                'Startpoint',0.01,'ReturnAsStructure',true);
            
            for ID = ["MX" "MM" "XX"]
                if UP.(ID).B < 0 || abs(UP.(ID).B) < 1e-5
                    disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - Parameter ' char(ID) '.B = ' num2str(UP.(ID).B,'%10.8e')])
                end
                if UP.(ID).C < 0 || abs(UP.(ID).C) < 1e-5
                    disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - Parameter ' char(ID) '.C = ' num2str(UP.(ID).C,'%10.8e')])
                end
                if UP.(ID).alpha < 0 || abs(UP.(ID).alpha) < 1e-5
                    disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - Parameter ' char(ID) '.alpha = ' num2str(UP.(ID).alpha,'%10.8e')])
                end
            end
            
        case 'TF'
            [UP.MX,UP.MM,UP.XX] = TF_Potential_Parameters(Settings);
            [U_MX, U_MM, U_XX] = TF_Potential_Generator(Settings,...
                'Startpoint',0.01,'ReturnAsStructure',true);
            
            for ID = ["MX" "MM" "XX"]
                if UP.(ID).B < 0
                    disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - Parameter ' char(ID) '.B = ' num2str(UP.(ID).B,'%10.8e')])
                end
                if UP.(ID).C < 0
                    disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - Parameter ' char(ID) '.C = ' num2str(UP.(ID).C,'%10.8e')])
                end
                if UP.(ID).D < 0
                    disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - Parameter ' char(ID) '.D = ' num2str(UP.(ID).D,'%10.8e')])
                end
                if UP.(ID).alpha < 0
                    disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - Parameter ' char(ID) '.alpha = ' num2str(UP.(ID).alpha,'%10.8e')])
                end
            end
        otherwise
            continue
    end
    
    % Check tables
    %% Grab the peaks and valleys of the MX attractive potential
    peaks_idx = islocalmax(U_MX.Total);
    valleys_idx = islocalmin(U_MX.Total);
    
    maxima_U = U_MX.Total(peaks_idx);
    minima_U = U_MX.Total(valleys_idx);
    
    if isempty(minima_U) % If no well minimum exists in MX interaction
        disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - No well minimum in MX interaction.'])
    elseif isempty(maxima_U) % If no peak exists in MX interaction
        % Do nothing, this is normal for JC potential and some BH/TF potentials
    else % Otherwise, a well minimum exists and at least one peak exists
        % ensure peak - well height is greater than specified threshold
        Threshold = Settings.MinExpWallHeight; % kJ/mol
        dU = maxima_U - minima_U;
        %Loss_add = log(1 + max(Threshold - dU,0)*Settings.BadFcnLossPenalty/Threshold);
        if max(Threshold - dU,0) > 0
            disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - MX well depth = ' num2str(dU) ' kJ/mol'])
        end
    end
%     plot(U_MX.r,U_MX.Total)
%     hold on
%     scatter(U_MX.r(peaks_idx),U_MX.Total(peaks_idx))
%     scatter(U_MX.r(valleys_idx),U_MX.Total(valleys_idx))
%     ylim([-1000 1000])
    
    %% Grab the peaks and valleys of the MM/XX potentials
    Us = [U_MM,U_XX];
    IDs = ["MM" "XX"];
    
    for jdx = 1:2
        U = Us(jdx);
        ID = IDs(jdx);
        
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
            %Loss_add = Loss_add + log(1 + max(dU - Threshold,0)*Settings.BadFcnLossPenalty);
            if max(dU - Threshold,0) > 0
                disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - ' char(ID) ' well depth = ' num2str(dU) ' kJ/mol'])
            end
            
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
            %Loss_add = Loss_add + log(1 + max(dU - Threshold,0)*Settings.BadFcnLossPenalty);
            if max(dU - Threshold,0) > 0
                disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - ' char(ID) ' well depth = ' num2str(dU) ' kJ/mol'])
            end
            
        elseif length(minima_U) == 1 % well minima is available but no peaks are visible, there must be a hidden peak to the right
            Threshold = Settings.MaxRepWellDepth; % kJ/mol
            dU = U.Total(end) - minima_U;
            %Loss_add = Loss_add + log(1 + max(dU - Threshold,0)*Settings.BadFcnLossPenalty);
            if max(dU - Threshold,0) > 0
                disp([outp{1} ' ' outp{2} ' Model ' outp{3} ' - ' char(ID) ' well depth = ' num2str(dU) ' kJ/mol'])
            end
        else
            % This should never be reached...
            warning('Possible issue with the potential!')
        end
    end
diary off

end