% Run_D4_Minimization
% Maintain_Symmetry = false; % If true, maintains the space group symmetry of the unit cell. If false, set the space group to S1

% opt_type: 'FULLOPTG' = Optimize both unit cell vectors and atom positions (fixed symmetry), 'ATOMONLY' = Atoms only, 'CELLONLY' = cell only, 'ITATOCEL' = iteratively cell, atoms, cell, atoms ...
% opt_conv = Convergence critera for optimization (10^-X hartrees)
% MaxCycles =  Maximum number of steps allowed for optimization
% Gradient_Tol_RMS = Convergence criterion on the Root-Mean-Squared of the gradient between optimization steps

% function Run_D4_Minimization(Salt,Structure,Theory,opt_type,P1_Symmetry,opt_conv,...
%     MaxCycles,Gradient_Tol_RMS,OptDir,JobName,Cry,Dispersion,Crystal_Input_text,Crystal_Geom_Template',OldCoordinates)
function Run_D4_Minimization(varargin)

if nargin > 0 % optionl directory input
    cd(varargin{1});
end

% Find mat file (or newest mat file if mulitple exist)
MATfiles = dir('*.mat');
if length(MATfiles) > 1
    times = [MATfiles(:).datenum];
    [~,time_idx] = max(times);
    MATfiles = MATfiles(time_idx);
end

load(MATfiles.name,'Salt','Structure',...
    'Theory','opt_type','P1_Symmetry','opt_conv','MaxCycles','Gradient_Tol_RMS',...
    'OptDir','JobName','Cry','Dispersion','Crystal_Input_text','Crystal_Geom_Template',...
    'Crystal_Input_Template','Coordinates');

% Gradient_Tol_RMS is input in terms of a.u. (Ha/Bohr) need in terms of
% Ha/Angstrom
Bohr_Ang = 1/0.529177210904; % Bohr/Angstrom
Gradient_Tol_RMS = Gradient_Tol_RMS*Bohr_Ang; % Gradient RMS convergence tolerance in Ha/Angstrom

%% Settings
Energy_Tol = 10^(-1*opt_conv); % a.u. Convergence criteria 1: must have change in energy between cycles is less than this value
Gradient_Tol_Max = 1.5*Gradient_Tol_RMS; % a.u. Convergence criteria 3: must have maximum cell parameter gradient less than this value for convergence
Gamma_Init = 1; % Initial step size coefficient (multiplied by numerical derivative along each dimension) used for moving along gradient
Gamma_Multiplier = 2; %1.5; % Scale Gamma by this value after each cycle (should be greater than 1)
Max_Step_size = 0.5; % Maximum possible step size in angstroms
alpha = 1/2;
beta = 0.1; % Reduce Gamma by this multiple when step size is too large.
MaxLineTries = 10; % Maximum number of tries to compute lower energy point before decreasing numerical derivative step size
Max_cycle_restarts = 10; % Maximum number of cycle restarts
OptPos_SingleStage = true; % When true, optimize positions and lattice parameters simultaneously by numerical gradients
Maintain_Symmetry = ~P1_Symmetry;

%% Begin Code
logfile = fullfile(OptDir,[JobName '.out']); % crystal log file
backuplogfile = fullfile(OptDir,[JobName '.out.old']); % backup log file
cd(OptDir); % Change to optimization directory

% Get Metal and Halide info from Current Salt
[Metal,Halide] = Separate_Metal_Halide(Salt);

% Grab the conventional atomic numbers from OldCoordinates
NumCoord = str2double(strtrim(regexp(Coordinates,'^[0-9]+\n','match','once')));
Coordinates = regexprep(Coordinates,'^[0-9]+\n','');
AtNums = textscan(Coordinates,'%f %*f %*f %*f');
MetalCAN = AtNums{1}(1);
HalideCAN = AtNums{1}(2);

% Select symmetry settings
if Maintain_Symmetry
    switch Structure
        case 'BetaBeO'
            DOF = {'a' 'c'};
            DOF_Units = {[' ' char(0197)] [' ' char(0197)]};
        case 'CsCl'
            DOF = {'a'};
            DOF_Units = {[' ' char(0197)]};
        case 'FiveFive'
            DOF = {'a' 'b' 'c'};
            DOF_Units = {[' ' char(0197)] [' ' char(0197)] [' ' char(0197)]};
        case 'Sphalerite'
            DOF = {'a'};
            DOF_Units = {[' ' char(0197)]};
        case 'NiAs'
            DOF = {'a' 'c'};
            DOF_Units = {[' ' char(0197)] [' ' char(0197)]};
        case 'Rocksalt'
            DOF = {'a'};
            DOF_Units = {[' ' char(0197)]};
        case 'Wurtzite'
            DOF = {'a' 'c'};
            DOF_Units = {[' ' char(0197)] [' ' char(0197)]};
    end
    DOF_txt = DOF;
else
    DOF = {'a' 'b' 'c' 'alpha' 'beta' 'gamma'};
    DOF_txt = {'a' 'b' 'c' char(945) char(946) char(947)};
    DOF_Units = {[' ' char(0197)] [' ' char(0197)] [' ' char(0197)] ...
        [' ' char(0176)] [' ' char(0176)] [' ' char(0176)]};
end

% Extra DOF for Optimization of Fractional Coordinates in 1 step
FC_DOFs = {[Metal '_x'] [Metal '_y'] [Metal '_z'] [Halide '_x'] [Halide '_y'] [Halide '_z']};
if OptPos_SingleStage && strcmp(opt_type,'FULLOPTG')
    DOF(end+1:end+6) = FC_DOFs;
    DOF_txt(end+1:end+6) = FC_DOFs;
    DOF_Units(end+1:end+6) = {'' '' '' '' '' ''};
end

% How many Degrees of freedom in gradient
N_DOF = length(DOF);

% Generate crystal17 run command
if ispc
    runcry17i = ['wsl source ~/.bashrc; runcry17 ' JobName];
    runcry17 = ['wsl source ~/.bashrc; runcry17 ' JobName ' ' JobName];
else
    usr = getenv('USERNAME');
    if strcmpi(usr,'user') % Lab PC
        runcry17i = ['source ~/.bashrc; runcry17 ' JobName];
        runcry17 = ['source ~/.bashrc; runcry17 ' JobName ' ' JobName];
    else % On server
        runcry17i = ['runcry17 ' JobName];
        runcry17 = ['runcry17 ' JobName ' ' JobName];
    end
end

TotalTimer = tic;
disp(['Beginning ' Salt ' ' Structure ' ' Theory '-' Dispersion ' Optimization...'])
disp('********************** Convergence Requirements ***********************')
disp(['Max ' char(916) 'E between cycles: ' num2str(Energy_Tol,'%4.4E') ' a.u.']);
disp(['Max RMS ' char(8711) 'E: ' num2str(Gradient_Tol_RMS,'%4.4E') ' a.u. / ' char(0197)]);
disp(['Max component ' char(8711) 'E: ' num2str(Gradient_Tol_Max,'%4.4E') ' a.u. / ' char(0197)]);

%% Begin Optimization Loop
% Loop broken when convergence criteria is met or max cycles reached
skip_results = false;
auto_h = true;
Cycle_restarts = 0;
Gamma = Gamma_Init;
ReCalc_Initial_E = true;
ResetIntegralClass = false; % Keeps track of whether or not the integral classification has been reset during the previous cycle.
for Index = 1:MaxCycles
    Restart_cycle = false;
    % Calculate reasonable step size
    if auto_h
        h = zeros(1,N_DOF);

        if Index == 1 
            % 1st crystal job
            [errcode,~] = system(runcry17i);
            
        elseif ResetIntegralClass && ReCalc_Initial_E
            % Restart with new classification
            [errcode,~] = system(runcry17);
            
        elseif ReCalc_Initial_E % New cycle
            Reset_Crystal_Input(NumCoord,MetalCAN,HalideCAN,OptDir,JobName,...
                Cry,Crystal_Input_text,Crystal_Geom_Template,Crystal_Input_Template,false);
            % Run crystal job
            [errcode,~] = system(runcry17);
        else % Don't need to recalculate initial energy when lower energy point found
            Reset_Crystal_Input(NumCoord,MetalCAN,HalideCAN,OptDir,JobName,...
                Cry,Crystal_Input_text,Crystal_Geom_Template,Crystal_Input_Template,false);
            errcode = 0;
        end
        
        if errcode ~= 0
            error(['Problem with CRYSTAL17, see log file at: ' OptDir filesep JobName '.out'])
        end
        
        % Grab DFT energy
        E_DFT = GrabEnergyDFT(OptDir,JobName,Cry.NF); % in a.u.
        
        % Grab D4 energy
        E_D4 = GrabEnergyD4(Salt,OptDir,JobName,Cry,Theory,Dispersion); % in a.u.
        
        % Total Energy
        E = E_DFT + E_D4; % DFT + D4 energy

        disp(['********************Cycle ' num2str(Index) ' Initial Conditions********************']);
        for Didx = 1:N_DOF
            if ~ismember(DOF{Didx},FC_DOFs)
                h(Didx) = nuderst(Cry.(DOF{Didx}));
                disp(['Lattice Parameter ' DOF_txt{Didx} ' = ' ...
                    num2str(Cry.(DOF{Didx}),'%2.8f') DOF_Units{Didx} '.' ...
                    ' Num. Der. Step Size: ' char(948) '(' DOF_txt{Didx} ') = ' num2str(h(Didx),'%2.8E') DOF_Units{Didx} '.']);
            else
                h(Didx) = nuderst(1);
            end
        end
        disp(['Fractional Coordinates for ' Metal ': ']);
        disp(num2str(Cry.FC_Metal(:,:),'%2.8f '))

        disp(['Fractional Coordinates for ' Halide ': ']);
        disp(num2str(Cry.FC_Halide(:,:),'%2.8f '))
        disp(['Initial E = ' num2str(E,'%4.10f') ' a.u.']);
        disp(['Step size coefficient: ' num2str(Gamma,'%4.4E') ])
        Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
        disp(['Time Elapsed: ' Telap])
        disp('******************************************************************')
        
        % Save log file
        copyfile(logfile,backuplogfile);
    else
        auto_h = true;
    end

    % Loop through each DOF
    Gradient = nan(1,N_DOF);
    for Didx = 1:N_DOF

        CurDOF = DOF{Didx};
        CurDOFtxt = DOF_txt{Didx};
        CurDOFunit = DOF_Units{Didx};
        delta = h(Didx);

        % Get energy at minus step
        CryMinus = Cry;

        if ismember(DOF{Didx},FC_DOFs)
            % For FC degrees of freedom
            MX = DOF{Didx}(1:2);
            xyz = strcmp(DOF{Didx}(end),{'x' 'y' 'z'});
            if strcmp(MX,pad(Metal,2,'_'))
                CryMinus.FC_Metal(1,xyz) = mod(CryMinus.FC_Metal(1,xyz) - delta,1);
            else
                CryMinus.FC_Halide(1,xyz) = mod(CryMinus.FC_Halide(1,xyz) - delta,1);
            end

            % Grab asymmetric unit FC
            AsymFC_Metal = CryMinus.FC_Metal(1,:);
            AsymFC_Halide = CryMinus.FC_Halide(1,:);

            % Reseat fractional coordinates outside of asymmetric unit
            [CryMinus.FC_Metal,CryMinus.FC_Halide] = ...
                UnitCell_FractionalCoords(AsymFC_Metal,AsymFC_Halide,Structure);
        else
            % For lattice parameter degrees of freedom
            CryMinus.(CurDOF) = CryMinus.(CurDOF) - delta;

            % Update transformation matrix
            CryMinus.Transform = GenTransformMatrix(CryMinus);

            % Maintain symmetry if set
            if Maintain_Symmetry && strcmp(CurDOF,'a')
                CryMinus = SymmetryAdapt(CryMinus,Structure); %#ok<*UNRCH>
            end
        end
        
        % Reset DFT input
        Reset_Crystal_Input(NumCoord,MetalCAN,HalideCAN,OptDir,JobName,...
            CryMinus,Crystal_Input_text,Crystal_Geom_Template,Crystal_Input_Template,false);
        
        % Run crystal job at new point
        [errcode,~] = system(runcry17);
        
        if errcode ~= 0
            error(['Problem with CRYSTAL17, see log file at: ' OptDir filesep JobName '.out'])
        end
        
        % Grab DFT energy
        E_DFT_Minus = GrabEnergyDFT(OptDir,JobName,CryMinus.NF); % in a.u.
        
        % Grab D4 energy
        E_D4_Minus = GrabEnergyD4(Salt,OptDir,JobName,CryMinus,Theory,Dispersion); % in a.u.
        
        % Total Energy
        E_Minus = E_DFT_Minus + E_D4_Minus; % DFT + D4 energy

        Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
        disp(['Finite step ' CurDOFtxt ' -' num2str(delta,'%2.2E') ...
            CurDOFunit ' Calculated. ' char(916) 'E(' CurDOFtxt '-' char(948) ...
            ') = ' num2str(E_Minus-E,'%+4.4E') ' a.u. Time Elapsed: ' Telap]);

        % Get energy at plus step
        CryPlus = Cry;

        if ismember(DOF{Didx},FC_DOFs)
            % For FC degrees of freedom
            MX = DOF{Didx}(1:2);
            xyz = strcmp(DOF{Didx}(end),{'x' 'y' 'z'});
            if strcmp(MX,pad(Metal,2,'_'))
                CryPlus.FC_Metal(1,xyz) = mod(CryPlus.FC_Metal(1,xyz) + delta,1);
            else
                CryPlus.FC_Halide(1,xyz) = mod(CryPlus.FC_Halide(1,xyz) + delta,1);
            end

            % Grab asymmetric unit FC
            AsymFC_Metal = CryPlus.FC_Metal(1,:);
            AsymFC_Halide = CryPlus.FC_Halide(1,:);

            % Reseat fractional coordinates outside of asymmetric unit
            [CryPlus.FC_Metal,CryPlus.FC_Halide] = ...
                UnitCell_FractionalCoords(AsymFC_Metal,AsymFC_Halide,Structure);
        else
            % For lattice parameter degrees of freedom
            CryPlus.(CurDOF) = CryPlus.(CurDOF) + delta;

            % Update transformation matrix
            CryPlus.Transform = GenTransformMatrix(CryPlus);

            % Maintain symmetry if set
            if Maintain_Symmetry && strcmp(CurDOF,'a')
                CryPlus = SymmetryAdapt(CryPlus,Structure); %#ok<*UNRCH>
            end
        end

        % Generate Energy at new point       
        % Reset DFT input
        Reset_Crystal_Input(NumCoord,MetalCAN,HalideCAN,OptDir,JobName,...
            CryPlus,Crystal_Input_text,Crystal_Geom_Template,Crystal_Input_Template,false);
        
        % Run crystal job
        [errcode,~] = system(runcry17);
        
        if errcode ~= 0
            error(['Problem with CRYSTAL17, see log file at: ' OptDir filesep JobName '.out'])
        end
        
        % Grab DFT energy
        E_DFT_Plus = GrabEnergyDFT(OptDir,JobName,CryPlus.NF); % in a.u.
        
        % Grab D4 energy
        E_D4_Plus = GrabEnergyD4(Salt,OptDir,JobName,CryPlus,Theory,Dispersion); % in a.u.
        
        % Total Energy
        E_Plus = E_DFT_Plus + E_D4_Plus; % DFT + D4 energy
        
        Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
        disp(['Finite step ' CurDOFtxt ' +' num2str(delta,'%2.2E') ...
            CurDOFunit ' Calculated. ' char(916) 'E(' CurDOFtxt '+' char(948) ...
            ') = ' num2str(E_Plus-E,'%+4.4E') ' a.u. Time Elapsed: ' Telap]);
        
        % 3-Point Numerical derivative wrt current DOF
        Gradient(Didx) = (1/(2*delta))*(-E_Minus + E_Plus);
    end
    disp([char(8711) 'E = ' num2str(Gradient,'%+4.4E  ') ' a.u. / ' char(0197)])
    
    % Check for convergence on cycle 1
    if ( Index == 1 || ResetIntegralClass ) && (rms(Gradient) < Gradient_Tol_RMS) && ...
            (max(abs(Gradient)) < Gradient_Tol_Max) && ~strcmp(opt_type,'FULLOPTG')
        % If gradient convergence criteria are met, end loop
        disp(['Energy convergence reached on cycle ' num2str(Index) ' gradient test.'])
        
        % save crystal log file of energy minimum
        delete(logfile)
        copyfile(backuplogfile,logfile);
        break
    end

    % Move one step in direction of steepest descent
    for StepInd = 1:MaxLineTries

        CryNew = Cry;

        for Didx = 1:N_DOF

            CurDOF = DOF{Didx};
            CurGrad = Gradient(Didx);

            % Step size for this degree of freedom along gradient
            % component
            ShiftStep = -sign(CurGrad)*min(abs(Gamma*CurGrad),Max_Step_size);

            if ismember(CurDOF,FC_DOFs)
                % For FC degrees of freedom
                MX = CurDOF(1:2);
                xyz = strcmp(CurDOF(end),{'x' 'y' 'z'});
                if strcmp(MX,pad(Metal,2,'_'))
                    CryNew.FC_Metal(1,xyz) = mod(CryNew.FC_Metal(1,xyz) + ShiftStep,1);
                else
                    CryNew.FC_Halide(1,xyz) = mod(CryNew.FC_Halide(1,xyz) + ShiftStep,1);
                end
            else
                % For lattice parameter degrees of freedom
                CryNew.(CurDOF) = Cry.(CurDOF) + ShiftStep;
            end
        end

        if OptPos_SingleStage
            % Grab asymmetric unit FC
            AsymFC_Metal = CryNew.FC_Metal(1,:);
            AsymFC_Halide = CryNew.FC_Halide(1,:);

            % Reseat fractional coordinates outside of asymmetric unit
            [CryNew.FC_Metal,CryNew.FC_Halide] = ...
                UnitCell_FractionalCoords(AsymFC_Metal,AsymFC_Halide,Structure);
        end

        if Maintain_Symmetry
            CryNew = SymmetryAdapt(CryNew,Structure);
        end

        % Update transformation matrix
        CryNew.Transform = GenTransformMatrix(CryNew);
        
        % Generate Energy at new point
        % Reset DFT input
        Reset_Crystal_Input(NumCoord,MetalCAN,HalideCAN,OptDir,JobName,...
            CryNew,Crystal_Input_text,Crystal_Geom_Template,Crystal_Input_Template,false);
        
        % Run crystal job
        [errcode,~] = system(runcry17);
        
        if errcode ~= 0
            error(['Problem with CRYSTAL17, see log file at: ' OptDir filesep JobName '.out'])
        end
        
        % Grab DFT energy
        E_DFT_New = GrabEnergyDFT(OptDir,JobName,CryNew.NF); % in a.u.
        
        % Grab D4 energy
        E_D4_New = GrabEnergyD4(Salt,OptDir,JobName,CryNew,Theory,Dispersion); % in a.u.
        
        % Total Energy
        E_New = E_DFT_New + E_D4_New; % DFT + D4 energy

        if E_New > E - min(alpha*Gamma*norm(Gradient)^2,50)
            Gamma = beta*Gamma;

            Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
            disp(['Optimizing step size... Time Elapsed: ' Telap])
        else
            Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
            disp(['Finite step in direction of steepest descent calculated. Time Elapsed: ' Telap]);
            disp([char(916) 'E = ' num2str(E_New-E,'%+4.4E') ' a.u.']);
            disp(['RMS ' char(8711) 'E = ' num2str(rms(Gradient),'%4.4E') ' a.u. / ' char(0197)])
            disp(['|Max(' char(8711) 'E)| = ' num2str(max(abs(Gradient)),'%4.4E') ' a.u. / ' char(0197)])
            break
        end

        if StepInd >= MaxLineTries
            h = max(h./10,eps*1000);
            Gamma = Gamma_Init;
            Restart_cycle = true;
            break
        end
    end

    if Restart_cycle
        Cycle_restarts = Cycle_restarts+1;
        if Cycle_restarts >= Max_cycle_restarts
            disp(['Unable to find lower energy point after ' num2str(Max_cycle_restarts) ' attempts.'])
            disp('Stopping here.')
            if (rms(Gradient) < Gradient_Tol_RMS) && (max(abs(Gradient)) < Gradient_Tol_Max)
                disp('Gradient convergence criteria fulfilled.')
                disp(['Energy convergence reached after ' num2str(Index) ' cycles.'])
            else
                disp('Gradient convergence criteria NOT fulfilled.')
                disp(['Search haulted after ' num2str(Index) ' cycles.'])
                disp('Check initial conditions. Poor initial conditions may cause this.')
                disp('Removing output files.')
                skip_results = true;
            end

            break
        else
            auto_h = false;
            disp('Unable to find a point of lower energy in chosen direction, starting next cycle with modified numerical derivative step size.');
            continue
        end
    end

    Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
    disp('************************************')
    disp(['Cycle ' num2str(Index) ' complete.'])
    disp(['Time Elapsed: ' Telap])
    disp('************************************')

    % Check if energy converged AND integal classifications reset: Job
    % finished
    if (abs(E_New - E) < Energy_Tol) && ...
            (rms(Gradient) < Gradient_Tol_RMS) && ...
            (max(abs(Gradient)) < Gradient_Tol_Max) && ResetIntegralClass
        % If all convergence criteria are met, end loop
        Cry = CryNew;
        E = E_New;
        E_DFT = E_DFT_New;
        E_D4 = E_D4_New;        
        disp(['Energy convergence reached after ' num2str(Index) ' cycles.'])
        break
    % Check if energy converged NOT integral classifications reset: Restart
    % with new integral classification
    elseif (abs(E_New - E) < Energy_Tol) && ...
            (rms(Gradient) < Gradient_Tol_RMS) && ...
            (max(abs(Gradient)) < Gradient_Tol_Max) && ~ResetIntegralClass
        % If all convergence criteria are met, end loop
        Cry = CryNew;
        E = E_New;
        E_DFT = E_DFT_New;
        E_D4 = E_D4_New;
        ReCalc_Initial_E = true;
        
        disp(['Energy convergence conditions met.' newline 'Restarting with renewed integral classification...'])
        
        % Reset the d12 file with new central point for integral classification
        Crystal_Input_text = Reset_Crystal_Input(NumCoord,MetalCAN,HalideCAN,OptDir,JobName,...
            Cry,Crystal_Input_text,Crystal_Geom_Template,Crystal_Input_Template,true);
        
        Cycle_restarts = 0;
        ResetIntegralClass = true;
    % If energy converged but not gradients
    elseif (abs(E_New - E) < Energy_Tol)
        Gamma = Gamma_Init; % Re-initialize Gamma
        E = E_New;
        E_DFT = E_DFT_New;
        E_D4 = E_D4_New;
        Cry = CryNew;
        ResetIntegralClass = false;
        ReCalc_Initial_E = false;
        Cycle_restarts = 0;
    % Otherwise a normal decrease in energy
    elseif E_New < E
        Gamma = Gamma*Gamma_Multiplier;
        E = E_New;
        E_DFT = E_DFT_New;
        E_D4 = E_D4_New;
        Cry = CryNew;
        ResetIntegralClass = false;
        ReCalc_Initial_E = false;
        Cycle_restarts = 0;
    % If energy increases
    elseif E_New > E
        disp('Warning: Total energy increased after geometry optimization step.')
        disp('Reverting to previous conditions.');
        Gamma = Gamma*Gamma_Multiplier;
        ResetIntegralClass = false;
        ReCalc_Initial_E = true;
        Cycle_restarts = 0;
    end

    if Index < MaxCycles
        disp(['Beginning Cycle ' num2str(Index+1) '.'])
    else
        disp(['Convergence NOT reached after ' num2str(Index) ' cycles. Stopping.'])
        skip_results = true;
    end
end

Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
disp(['Completed: ' Salt ' ' Structure ' ' Theory '-' Dispersion ' Geometry Optimization. Time ' Telap])
if ~skip_results
    disp(['Final Optimized Energy is ' num2str(E,'%4.10f') ' a.u.'])
    for Didx = 1:N_DOF
        if ~ismember(DOF{Didx},FC_DOFs)
            disp(['Lattice Parameter ' DOF_txt{Didx} ' = ' ...
                num2str(Cry.(DOF{Didx}),'%2.8f') DOF_Units{Didx} '.']);
        end
    end
    disp(['Fractional Coordinates for ' Metal ': ']);
    disp(num2str(Cry.FC_Metal(:,:),'%2.8f '))
    disp(['Fractional Coordinates for ' Halide ': ']);
    disp(num2str(Cry.FC_Halide(:,:),'%2.8f '))
    disp('Components of Final Gradient:')
    disp(num2str(Gradient,'%+4.4E  '))
    
    % Save results to structure
    Data = Cry;
	Data.Toldee = opt_conv;
    
    % Want energy output per primitive cell
    if contained_in_cell(Structure,{'CsCl' 'Rocksalt' 'Sphalerite' 'Pair'})
        Data.E = E; % Energy per formula unit (a.u.)
        Data.E_DFT = E_DFT;
        Data.E_D4 = E_D4;
    elseif contained_in_cell(Structure,{'Wurtzite' 'NiAs'})
        Data.E = 2*E; % Energy per formula unit (a.u.)
        Data.E_DFT = 2*E_DFT;
        Data.E_D4 = 2*E_D4;
    elseif contained_in_cell(Structure,{'BetaBeO' 'FiveFive'})
        Data.E = 4*E; % Energy per formula unit (a.u.)
        Data.E_DFT = 4*E_DFT;
        Data.E_D4 = 4*E_D4;
    end
    % Save to .mat file
    save([OptDir filesep JobName '_OUT.mat'],'Data')
end

end

