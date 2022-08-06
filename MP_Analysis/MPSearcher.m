function Tm = MPSearcher(fun,T0,lb,ub,Settings)

if T0 < lb
    disp(['Initial T = ' num2str(T0,'%.4f') ' K is lower than specified lower bound of '  num2str(lb,'%.4f') ' K.' ...
        ' Setting initial T to lower bound.'])
    T = lb;
elseif T0 > ub
    disp(['Initial T = ' num2str(T0,'%.4f') ' K is higher than specified upper bound of '  num2str(ub,'%.4f') ' K.' ...
        ' Setting initial T to upper bound.'])
    T = ub;
else
    T = T0;
end

% Make initial function evalulation at T
[~,~,T_dat] = fun(T);

% Reclassify any Freeze points that are above known melting points
% Or any melt points that are below known freeze points
[T_Freeze_Max,TFMidx] = max(T_dat.T_Trace(logical(T_dat.Freeze_Trace)));
[T_Melt_Min,TMMidx] = min(T_dat.T_Trace(logical(T_dat.Melt_Trace)));
Potential_Tms = T_dat.T_Trace(~T_dat.Freeze_Trace & ~T_dat.Melt_Trace);

if isempty(T_Melt_Min)
    T_Melt_Min = ub;
end
if isempty(T_Freeze_Max)
    T_Freeze_Max = lb;
end

if any(T_dat.T_Freeze_Trace > T_Melt_Min) || any(Potential_Tms > T_Melt_Min)...
    || any(T_dat.T_Melt_Trace < T_Freeze_Max) || any(Potential_Tms < T_Freeze_Max)

    disp('Warning: updating classification of at least one freezing/melting point.')
    disp('This occurs when a new melted T is found to be lower than a previous freezing or unknown T')
    disp('or when a new freezing T is found to be higher than a previous melting or unknown T.')
    T_dat.Melt_Trace = T_dat.T_Trace >= T_Melt_Min;
    T_dat.Freeze_Trace = T_dat.T_Trace <= T_Freeze_Max;

    T_dat.T_Melt_Trace = T_dat.T_Trace(T_dat.Melt_Trace);
    T_dat.T_Freeze_Trace = T_dat.T_Trace(T_dat.Freeze_Trace);

    if isempty(TFMidx)
        T_dat.dT(1) = lb;
        T_dat.df_bracket(1) = nan;
    else
        T_dat.dT(1) = T_Freeze_Max;
        T_dat.df_bracket(1) = T_dat.df_Trace(TFMidx);
    end

    if isempty(TMMidx)
        T_dat.dT(2) = ub;
        T_dat.df_bracket(2) = nan;
    else
        T_dat.dT(2) = T_Melt_Min;
        T_dat.df_bracket(2) = T_dat.df_Trace(TMMidx);
    end

    copyfile(Settings.CurrentTFile,Settings.PrevTFile)
    save(Settings.CurrentTFile,'T_dat')
end

% Update the MP bracket
MP_check = (~T_dat.Freeze_Trace & ~T_dat.Melt_Trace);
if any(MP_check)
    highest_MP = max(T_dat.T_Trace(MP_check));
    lowest_MP = min(T_dat.T_Trace(MP_check));
    dT_melt_side = abs(T_dat.dT(2) - highest_MP);
    dT_freeze_side = abs(T_dat.dT(1) - lowest_MP);

    TBracket = max(dT_melt_side,dT_freeze_side)*2;
else
    TBracket = diff(T_dat.dT);
end

% Check for abort calculation flag
if T_dat.Alt_Structure
    Tm = T;
    return
end

while TBracket > Settings.BracketThreshold

    disp(['Current melting point bracket: ' num2str(TBracket,'%.4f') ...
        ' K, between T[Freeze] = ' num2str(T_dat.dT(1),'%.4f') ' K and T[Melt] = ' ...
        num2str(T_dat.dT(2),'%.4f') ' K.'])
    % Info
    disp(['Checked temperatures that melted: ' num2str(sort(T_dat.T_Melt_Trace),'%.4f K    ')])
    disp(['Checked temperatures that froze:  ' num2str(sort(T_dat.T_Freeze_Trace,'descend'),'%.4f K    ')])
    disp(['Checked indeterminate points:             ' num2str(sort(T_dat.T_Trace(MP_check)),'%.4f K    ')])
    
    if any(MP_check)

        T_Freeze_Max = T_dat.dT(1);
        T_Melt_Min = T_dat.dT(2);

        MP_Max = max(T_dat.T_Trace(MP_check));
        MP_Min = min(T_dat.T_Trace(MP_check));

        % Close the upper bound first
        if T_Melt_Min - MP_Max > Settings.BracketThreshold/2
            T = MP_Max + Settings.BracketThreshold/2;
        % If upper bound is closed, close the lower bound
        elseif MP_Min - T_Freeze_Max > Settings.BracketThreshold/2
            T = MP_Min - Settings.BracketThreshold/2;
        else
            error('It should not be possible to reach this!')
        end

    % If there are no known indeterminate points, attempt to narrow the bracket
    else
        T_Freeze_Max = max(T_dat.T_Trace(logical(T_dat.Freeze_Trace)));
        T_Melt_Min = min(T_dat.T_Trace(logical(T_dat.Melt_Trace)));

        if isempty(T_Freeze_Max) % If only melting points have been found so far
            % Update the stepsize
            StepSize = (Settings.MeshSizeMultiplier^(length(T_dat.T_Melt_Trace)))*Settings.InitialMeshSize;
            T = max(T_Melt_Min - StepSize,lb);
        elseif isempty(T_Melt_Min) % If only freezing points have been found so far
            % Update the stepsize
            StepSize = (Settings.MeshSizeMultiplier^(length(T_dat.T_Freeze_Trace)))*Settings.InitialMeshSize; 
            T = min(T_Freeze_Max + StepSize,ub);
        else
            % Once a melting point bracket has been established, narrow it
            if Settings.UseDerivativeWeighting && ~any(isnan(T_dat.df_bracket))
                W = (1/abs(T_dat.df_bracket(1))) + (1/abs(T_dat.df_bracket(2)));
                T = ((1/abs(T_dat.df_bracket(1)))*T_dat.dT(1) + (1/abs(T_dat.df_bracket(2)))*T_dat.dT(2))/W; % derivative-weighted average
            else
                T = mean(T_dat.dT); % Simple average
            end
        end
    end

    % Evaluate at the next test point
    [~,df,T_dat] = fun(T);

    % Reclassify any Freeze points that are above known melting points
    % Or any melt points that are below known freeze points
    [T_Freeze_Max,TFMidx] = max(T_dat.T_Trace(logical(T_dat.Freeze_Trace)));
    [T_Melt_Min,TMMidx] = min(T_dat.T_Trace(logical(T_dat.Melt_Trace)));
    Potential_Tms = T_dat.T_Trace(~T_dat.Freeze_Trace & ~T_dat.Melt_Trace);

    if isempty(T_Melt_Min)
        T_Melt_Min = ub;
    end
    if isempty(T_Freeze_Max)
        T_Freeze_Max = lb;
    end

    if any(T_dat.T_Freeze_Trace > T_Melt_Min) || any(Potential_Tms > T_Melt_Min)...
        || any(T_dat.T_Melt_Trace < T_Freeze_Max) || any(Potential_Tms < T_Freeze_Max)
        
        disp('Warning: updating classification of at least one freezing/melting point.')
        disp('This occurs when a new melted T is found to be lower than a previous freezing or unknown T')
        disp('or when a new freezing T is found to be higher than a previous melting or unknown T.')
        T_dat.Melt_Trace = T_dat.T_Trace >= T_Melt_Min;
        T_dat.Freeze_Trace = T_dat.T_Trace <= T_Freeze_Max;
        
        T_dat.T_Melt_Trace = T_dat.T_Trace(T_dat.Melt_Trace);
        T_dat.T_Freeze_Trace = T_dat.T_Trace(T_dat.Freeze_Trace);

        if isempty(TFMidx)
            T_dat.dT(1) = lb;
            T_dat.df_bracket(1) = nan;
        else
            T_dat.dT(1) = T_Freeze_Max;
            T_dat.df_bracket(1) = T_dat.df_Trace(TFMidx);
        end

        if isempty(TMMidx)
            T_dat.dT(2) = ub;
            T_dat.df_bracket(2) = nan;
        else
            T_dat.dT(2) = T_Melt_Min;
            T_dat.df_bracket(2) = T_dat.df_Trace(TMMidx);
        end

        copyfile(Settings.CurrentTFile,Settings.PrevTFile)
        save(Settings.CurrentTFile,'T_dat')
    end

    if T == ub && (df < 0) % derivative less than 0 indicates the system froze
        disp('System froze at upper T bound. Unable to bracket melting point within bounds!')
        Tm = T;
        return
    elseif T == lb && df > 0  % derivative greater than 0 indicates the system melted
        disp('System melted at lower T bound. Unable to bracket melting point within bounds!')
        Tm = T;
        return
    end
    if T_dat.Alt_Structure % Check for abort calculation flag
        Tm = T;
        return
    end

    % Update the MP bracket
    MP_check = (~T_dat.Freeze_Trace & ~T_dat.Melt_Trace);
    if any(MP_check)
        highest_MP = max(T_dat.T_Trace(MP_check));
        lowest_MP = min(T_dat.T_Trace(MP_check));
        dT_melt_side = abs(T_dat.dT(2) - highest_MP);
        dT_freeze_side = abs(T_dat.dT(1) - lowest_MP);

        TBracket = max(dT_melt_side,dT_freeze_side)*2;
    else
        TBracket = diff(T_dat.dT);
    end
end

if any(MP_check)
    Tms = sort(T_dat.T_Trace(MP_check));
    idx = max(floor(length(Tms)/2),1);

    Tm = Tms(idx);
else
    Tm = T_dat.dT(1); % Report the freezing side of the melting point range by convention
end


disp(['Melting point bracketed to within: ' num2str(TBracket,'%.4f') ' K.'])
disp(['This is smaller than the specified M.P. temperature bracket threshold of ' ...
    num2str(Settings.BracketThreshold,'%.4f') ' K.'])
disp(repmat('*',1,40));
disp(['Estimated melting point is: ' num2str(mean(T_dat.dT),'%.4f') ' K.'])
disp(repmat('*',1,40));
disp('Summary of melting point trajectory:')
disp(repmat('*',1,40));
for idx = 1:length(T_dat.T_Trace)
    if T_dat.df_Trace(idx) < 0
        cls = 'freezing';
    elseif T_dat.df_Trace(idx) > 0
        cls = 'melting';
    else
        cls = 'indeterminate';
    end

    if T_dat.Freeze_Trace(idx)
        cls2 = 'freezing';
    elseif T_dat.Melt_Trace(idx)
        cls2 = 'melting';
    else
        cls2 = 'indeterminate';
    end

    disp(['At T = ' num2str(T_dat.T_Trace(idx),'%.4f') ' K, detected as ' cls ' and classified as ' cls2 '.'])
end
disp(repmat('*',1,40));


    
end