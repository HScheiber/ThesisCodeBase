function Tm = MPSearcher_old(fun,T,lb,ub,Settings)

% Initialize parameters
StepSize = Settings.InitialMeshSize;

if T < lb
    disp(['Initial T = ' num2str(T,'%.4f') ' K is lower than specified lower bound of '  num2str(lb,'%.4f') ' K.' ...
        ' Setting initial T to lower bound.'])
    T = lb;
elseif T > ub
    disp(['Initial T = ' num2str(T,'%.4f') ' K is higher than specified upper bound of '  num2str(ub,'%.4f') ' K.' ...
        ' Setting initial T to upper bound.'])
    T = ub;
end

%% Step 1: Make initial bracket of melting point
% Evaluate the input function at initial point
% Function returns a value, a derivative, and a user data structure
[f,df,T_dat] = fun(T);
df_prev = df;

% Check for calculation completion
[calc_finished,Tm] = MPCheckFun(T_dat);
if calc_finished
    return
end

while sign(df) == sign(df_prev)
    
    if T == ub && (df < 0) % derivative less than 0 indicates the system froze
        disp('System froze at upper T bound. Unable to bracket melting point within bounds!')
        Tm = T;
        return
    elseif T == lb && df > 0  % derivative greater than 0 indicates the system melted
        disp('System melted at lower T bound. Unable to bracket melting point within bounds!')
        Tm = T;
        return
    end
    
    % Update temperature based on the sign of the derivative
    if df < 0 % derivative less than 0 indicates the system froze
        T = min(T + StepSize,ub); % New temperature
    elseif df > 0  % derivative greater than 0 indicates the system melted
        T = max(T - StepSize,lb); % New temperature
    end
    df_prev = df;
    [f,df,T_dat] = fun(T);
    
    [calc_finished,Tm] = MPCheckFun(T_dat);
    if calc_finished
        return
    end
    if (1/f) < Settings.CheckTime
        StepSize = StepSize*Settings.MeshSizeMultiplier; % Update stepsize
    else
        StepSize = Settings.InitialMeshSize;
    end
        
end

%% Step 2: Initial bracket is generated. Check midpoint until bracket is small
TBracket = diff(T_dat.dT);
disp(['Current melting point bracket: ' num2str(TBracket,'%.4f') ...
    ' K, between T[Freeze] = ' num2str(T_dat.dT(1),'%.4f') ' K and T[Melt] = ' num2str(T_dat.dT(2),'%.4f') ' K.'])
disp('Beginning melting point bracket reduction algorithm.')

while TBracket > Settings.BracketThreshold
    % Estimate a new melting point guess between the current MP bracket
    T_prev = T;
    if Settings.UseDerivativeWeighting
        W = (1/abs(T_dat.df_bracket(1))) + (1/abs(T_dat.df_bracket(2)));
        T = ((1/abs(T_dat.df_bracket(1)))*T_dat.dT(1) + (1/abs(T_dat.df_bracket(2)))*T_dat.dT(2))/W; % derivative-weighted average
    else
        T = mean(T_dat.dT); % Simple average
    end
    
    % Ensure step size is not less than MinStepSize parameter
    if abs(T_prev - T) < Settings.MinStepSize
        if df < 0 % derivative less than 0 indicates the system froze on previous step
            T = min(T_prev + Settings.MinStepSize,ub);
        elseif df > 0  % derivative greater than 0 indicates the system melted on previous step
            T = max(T_prev - Settings.MinStepSize,lb);
        end
    end
    
    [f,df,T_dat] = fun(T); % Evaluate function at new guess temperature
    TBracket = diff(T_dat.dT); % Update MP bracket
    
    % Check if melting point is found
    [calc_finished,Tm] = MPCheckFun(T_dat);
    if calc_finished
        return
    end
end

% Check if melting point is found (this may be reached on a job restart)
[calc_finished,Tm] = MPCheckFun(T_dat);
if calc_finished
    return
end

%% Once here, the bracket on melting temperature will be smaller than the bracket threshold, but no melting point is yet found
if Settings.UseDerivativeWeighting
    W = (1/abs(T_dat.df_bracket(1))) + (1/abs(T_dat.df_bracket(2)));
    Tm = ((1/abs(T_dat.df_bracket(1)))*T_dat.dT(1) + (1/abs(T_dat.df_bracket(2)))*T_dat.dT(2))/W; % derivative-weighted average
else
    Tm = mean(T_dat.dT); % Simple average
end

disp(['Melting point bracketed to within: ' num2str(TBracket,'%.4f') ' K.'])
disp(['This is smaller than the specified M.P. temperature bracket threshold of ' ...
    num2str(Settings.BracketThreshold,'%.4f') ' K.'])
disp(['Estimated melting point is: ' num2str(Tm,'%.4f') ' K.'])
Tm = T_dat.dT(1); % Use lower bound as estimated melting point
% disp('Performing final check at estimated melting point')
% Settings.IgnoreBounds = true;
% fun = @(T)Melting_Point_Check(T,Settings);
% fun(Tm);
% disp('Final M.P. check complete at estimated melting point.')

function [calc_finished,Tm] = MPCheckFun(T_dat)
    MP_check = (~T_dat.Freeze_Trace & ~T_dat.Melt_Trace);
    if any(MP_check)
        MPs = T_dat.T_Trace(MP_check);
        Tm = MPs(end);
        for kdx = 1:length(MPs)
            disp(['Melting point found at Tm = ' num2str(MPs(kdx),'%.4f') ' K.'])
        end
        disp(['Current melting point bracket: ' num2str(diff(T_dat.dT),'%.4f') ...
            ' K, between T[Freeze] = ' num2str(T_dat.dT(1),'%.4f') ' K and T[Melt] = ' num2str(T_dat.dT(2),'%.4f') ' K.'])
        Tm_max = max(MPs);
        Tm_min = min(MPs);

        if T_dat.dT(2) - Tm_max > Settings.BracketThreshold/2 || Tm_min - T_dat.dT(1) > Settings.BracketThreshold/2
            disp('Melting point bracket is not yet sufficiently small.')
            disp('Reducing melting point bracket further...')
            [T_dat] = ReduceMPBracket(fun,Settings,T_dat);
            MP_check = (~T_dat.Freeze_Trace & ~T_dat.Melt_Trace);
            MPs = T_dat.T_Trace(MP_check);
            if isempty(MPs)
                Tm = nan;
                calc_finished = false;
                return
            else
                Tm = MPs(end);
            end
        end
        disp(['Melting point bracketed to within: ' num2str(diff(T_dat.dT),'%.4f') ' K.'])
        disp(['The user-specified M.P. temperature bracket threshold of ' ...
            num2str(Settings.BracketThreshold,'%.4f') ' K has been reached, based on all calculated melting points.'])
        for kdx = 1:length(MPs)
            disp(['Final melting point(s): Tm = ' num2str(MPs(kdx),'%.4f') ' K.'])
        end
        calc_finished = true;
    else
        Tm = nan;
        calc_finished = false;
    end
end

function [T_dat] = ReduceMPBracket(fun,Settings,T_dat)
    
    MP_check = (~T_dat.Freeze_Trace & ~T_dat.Melt_Trace );
    MPs = T_dat.T_Trace(MP_check);
    Tm_min = min(MPs);
    Tm_max = max(MPs);
    
    % Narrow in the lower bound
    while (Tm_min - T_dat.dT(1)) - Settings.BracketThreshold/2 > sqrt(eps) % System has not been sufficiently bracketed on freezing side of MP
        T_LB_Test = Tm_min - Settings.BracketThreshold/2;
        
        [f,df,T_dat] = fun(T_LB_Test);
        
        % Three possible cases:
        if f <= Settings.FunctionTolerance % System did not melt or freeze
            disp('Additional MP found.')
        elseif df < 0 % derivative less than 0 indicates the system froze on previous step
            break
        elseif df > 0  % derivative greater than 0 indicates the system melted on previous step
            % This shouldn't be reachable, and indicates the system is
            % not at the thermodynamic limit...
            disp('Warning: System melted after lowering temperature from MP!')
            disp(['Changing any freezing points above ' num2str(T_LB_Test,'%.4f') ' K to melting points.'])      
            T_dat.Melt_Trace((~T_dat.Freeze_Trace & ~T_dat.Melt_Trace) & (T_dat.T_Trace > T_LB_Test)) = true;
            T_dat.T_Melt_Trace = T_dat.T_Trace(logical(T_dat.Melt_Trace));
            T_dat.df_Trace(logical(T_dat.Melt_Trace) & T_dat.df_Trace <= 0) = 1;
            if isempty(max(T_dat.T_Freeze_Trace))
                dTl = lb;
            else
                dTl = max(T_dat.T_Freeze_Trace);
            end
            if isempty(min(T_dat.T_Melt_Trace))
                dTu = ub;
            else
                dTu = min(T_dat.T_Melt_Trace);
            end
            T_dat.dT = [dTl dTu];
            
            copyfile(Settings.CurrentTFile,Settings.PrevTFile)
            save(Settings.CurrentTFile,'T_dat')
        end
        % re-calculate Tm-min
        MP_check = (~T_dat.Freeze_Trace & ~T_dat.Melt_Trace );
        MPs = T_dat.T_Trace(MP_check);
        Tm_min = min(MPs);
    end
    disp('Lower bracket on melting point established.')
    
    while (T_dat.dT(2) - Tm_max) - Settings.BracketThreshold/2 > sqrt(eps) % System has not been sufficiently bracketed on melting side of MP
        T_UB_Test = Tm_max + Settings.BracketThreshold/2;
        
        [f,df,T_dat] = fun(T_UB_Test);
        
        % Three possible cases:
        if f <= Settings.FunctionTolerance % System did not melt or freeze
            disp('Additional MP found.')
        elseif df < 0 % derivative less than 0 indicates the system froze on previous step
            % This shouldn't be reachable, and indicates a previous step
            % was delcared a melting point in error
            disp('Warning: System froze after raising temperature from MP!')
            disp(['Changing any melting points below ' num2str(T_UB_Test,'%.4f') ' K to freezing points.'])         
            T_dat.Freeze_Trace((~T_dat.Freeze_Trace & ~T_dat.Melt_Trace) & (T_dat.T_Trace < T_UB_Test)) = true;
            T_dat.T_Freeze_Trace = T_dat.T_Trace(logical(T_dat.Freeze_Trace));
            T_dat.df_Trace(logical(T_dat.Freeze_Trace) & T_dat.df_Trace >= 0) = -1;
            if isempty(max(T_dat.T_Freeze_Trace))
                dTl = lb;
            else
                dTl = max(T_dat.T_Freeze_Trace);
            end
            if isempty(min(T_dat.T_Melt_Trace))
                dTu = ub;
            else
                dTu = min(T_dat.T_Melt_Trace);
            end
            T_dat.dT = [dTl dTu];
            
            copyfile(Settings.CurrentTFile,Settings.PrevTFile)
            save(Settings.CurrentTFile,'T_dat')
        elseif df > 0  % derivative greater than 0 indicates the system melted on previous step
            break
        end
        % re-calculate Tm-max
        MP_check = (~T_dat.Freeze_Trace & ~T_dat.Melt_Trace);
        MPs = T_dat.T_Trace(MP_check);
        Tm_max = max(MPs);
    end
    disp('Upper bracket on melting point established.')
end

end