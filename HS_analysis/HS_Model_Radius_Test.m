Salt = 'MX';
Structures = {'Wurtzite'};
Model = 'HS';
TFParamset = 0;
C6_Damp = Init_C6Damping_Object;
Scaling = Init_Scaling_Object;

CR_Damp.MX.r_d = 0.05; % This is the value of the sigmoid's midpoint in nm. Set to < 0 to disable close range damping
CR_Damp.MX.b = 1e6; % sigmoid "steepness" for damping. Set to < 0 to disable close range damping
CR_Damp.MM.r_d = 0.05; % 0.15 0.21; 0.26;
CR_Damp.MM.b  = 1e6;% 100 75
CR_Damp.XX.r_d = 0.05; 
CR_Damp.XX.b  = 1e6;


% Repulsive Wall Height (HS only)
% Units: kJ/mol
Scaling.Z.All = 1e10;
Scaling.Z.MX = 1;
Scaling.Z.MM = 1;
Scaling.Z.XX = 1;

% Repulsive Steepness Parameter b (HS only)
% Units: inverse nanometers
Scaling.b.All = Inf;
Scaling.b.MX = 1;
Scaling.b.MM = 1;
Scaling.b.XX = 1;

% Hard Sphere Distance Parameter r_d (HS only)
% Units are nanometers
Scaling.r_d.All = 1;
Scaling.r_d.MM = 0.09;
rM_over_rX = 0.1:0.005:1.0;
XX_Size = (Scaling.r_d.MM./rM_over_rX);

% Redo if some failed
Redo = true;
Redo_rM_over_rX = [0.1:0.005:0.13 0.795 0.82 0.85 0.915 0.93:0.005:0.945 0.98];
Continuation = false; % set to true to restart from the final geometry of the last run
tol = 1e-5;

Parallel_Struct_Min = false;
OptPos = false;
Maintain_Symmetry = true;
[GAdjust_MX,GAdjust_MM,GAdjust_XX] = Init_GAdjust_Object;

N_Scale = length(XX_Size);
N_Structure = length(Structures);

load('Rocksalt_vs_NiAs_HS_Model_M0.09.mat','-mat')


p = gcp('nocreate');
if isempty(p)
    ppool = parpool(8);
elseif ~(p.Connected && p.NumWorkers == 8)
    delete(p)
    ppool = parpool(8);
else
    ppool = p;
end

%% Generate Data
for idx = 1:N_Structure
    Structure = Structures{idx};
    disp('*************************************')
    disp(['CURRENT STRUCTURE: ' Structure])
    disp('*************************************')
    
    CoulSR    = zeros(1,N_Scale);
    CoulLR    = zeros(1,N_Scale);
    Pot       = zeros(1,N_Scale);
    CoulSR_MM = zeros(1,N_Scale);
    CoulSR_MX = zeros(1,N_Scale);
    CoulSR_XX = zeros(1,N_Scale);
    a = zeros(1,N_Scale);
    b = zeros(1,N_Scale);
    c = zeros(1,N_Scale);
    
    if ~Redo
        E.(Structure).CoulSR = zeros(1,N_Scale);
        E.(Structure).CoulLR = zeros(1,N_Scale);
        E.(Structure).Pot    = zeros(1,N_Scale);
        E.(Structure).CoulSR_MM = zeros(1,N_Scale);
        E.(Structure).CoulSR_MX = zeros(1,N_Scale);
        E.(Structure).CoulSR_XX = zeros(1,N_Scale);
        E.(Structure).a = zeros(1,N_Scale);
        E.(Structure).b = zeros(1,N_Scale);
        E.(Structure).c = zeros(1,N_Scale);
        
        N_Scale_redo = N_Scale;
        redo_indeces = 1:N_Scale;
    else
        N_Scale_redo = length(Redo_rM_over_rX);
        LocA = ~ismembertol(rM_over_rX,Redo_rM_over_rX,tol);
        redo_indeces = find(~LocA);
        
        if length(E.(Structure).CoulSR) < N_Scale
            E.(Structure).CoulSR(end+1:N_Scale) = 0;
            E.(Structure).CoulLR(end+1:N_Scale) = 0;
            E.(Structure).Pot(end+1:N_Scale) = 0;
            E.(Structure).CoulSR_MM(end+1:N_Scale) = 0;
            E.(Structure).CoulSR_MX(end+1:N_Scale) = 0;
            E.(Structure).CoulSR_XX(end+1:N_Scale) = 0;
            E.(Structure).a(end+1:N_Scale) = 0;
            E.(Structure).b(end+1:N_Scale) = 0;
            E.(Structure).c(end+1:N_Scale) = 0;
        end
        
        CoulSR(LocA)    = E.(Structure).CoulSR(LocA);
        CoulLR(LocA)    = E.(Structure).CoulLR(LocA);
        Pot(LocA)       = E.(Structure).Pot(LocA);
        CoulSR_MM(LocA) = E.(Structure).CoulSR_MM(LocA);
        CoulSR_MX(LocA) = E.(Structure).CoulSR_MX(LocA);
        CoulSR_XX(LocA) = E.(Structure).CoulSR_XX(LocA);
        a(LocA)         = E.(Structure).a(LocA);
        b(LocA)         = E.(Structure).b(LocA);
        c(LocA)         = E.(Structure).c(LocA);
    end
    
    
    f(1:N_Scale_redo) = parallel.FevalFuture;
    
    for jdx = 1:N_Scale_redo
        E_idx = redo_indeces(jdx);
        Scale = Scaling;
        Scale.r_d.XX = XX_Size(E_idx);
        Scale.r_d.MX = mean([Scale.r_d.MM Scale.r_d.XX]);
        Cur_rM_over_rX = Scale.r_d.MM / Scale.r_d.XX;
        
        if Redo && Continuation
            switch Structure
                case {'Rocksalt' 'Sphalerite'}
                    x0 = [E.(Structure).a(E_idx)/sqrt(2) ...
                          E.(Structure).b(E_idx)/sqrt(2) ...
                          E.(Structure).c(E_idx)/sqrt(2)];
                case 'FiveFive'
                    x0 = [E.(Structure).b(E_idx)/(3/sqrt(3)) ...
                          E.(Structure).c(E_idx) ...
                          E.(Structure).a(E_idx)];
                otherwise
                    x0 = [E.(Structure).a(E_idx) ...
                          E.(Structure).b(E_idx) ...
                          E.(Structure).c(E_idx)];
            end
        else
            x0 = [];
        end
        
        disp(['Current rM/rX: ' num2str(Cur_rM_over_rX)])
            
        f(jdx) = parfeval(ppool,@Structure_Minimization,1,Salt,...
            Structure,Model,TFParamset,C6_Damp,CR_Damp,Scale,...
            Parallel_Struct_Min,OptPos,Maintain_Symmetry,...
            GAdjust_MX,GAdjust_MM,GAdjust_XX,'Continuation',x0);
    end
    wait(f);
    
    % Update and save
    for jdx = 1:N_Scale_redo
        E_idx = redo_indeces(jdx);
        
        Output = f(jdx).fetchOutputs;
        CoulSR(E_idx)    = Output.CoulSR;
        CoulLR(E_idx)    = Output.CoulLR;
        Pot(E_idx)       = Output.Pot;
        CoulSR_MM(E_idx) = Output.CoulSR_MM;
        CoulSR_MX(E_idx) = Output.CoulSR_MX;
        CoulSR_XX(E_idx) = Output.CoulSR_XX;
        a(E_idx)         = Output.a;
        b(E_idx)         = Output.b;
        c(E_idx)         = Output.c;
    end
    
    
    % Update and save
    E.(Structure).CoulSR    = CoulSR;
    E.(Structure).CoulLR    = CoulLR;
    E.(Structure).Pot       = Pot;
    E.(Structure).CoulSR_MM = CoulSR_MM;
    E.(Structure).CoulSR_MX = CoulSR_MX;
    E.(Structure).CoulSR_XX = CoulSR_XX;
    E.(Structure).a = a;
    E.(Structure).b = b;
    E.(Structure).c = c;
    save('Rocksalt_vs_NiAs_HS_Model_M0.09.mat',"E")
end