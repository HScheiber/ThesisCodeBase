% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETERS AND kJ/mol


% GAdjust are N x 3 arrays of gaussian parameters
% (i , 1) is the Gaussian height of the ith adjustment (may be negative or
% positive)
% (i , 2) is the center point of the ith Gaussian (should be positive)
% (i , 3) is the standard deviation or width (negative and positive values
% are the same)
function U = HS_Potential_Generator(Settings,varargin)

% Optional inputs
p = inputParser;
p.FunctionName = 'HS_Potential_Generator';
addOptional(p,'PlotType','full')
addOptional(p,'Startpoint',0);
addOptional(p,'Plotswitch',false);
addOptional(p,'MDP_Minimize',false);
parse(p,varargin{:});
PlotType = p.Results.PlotType;
Startpoint = p.Results.Startpoint;
Plotswitch = p.Results.Plotswitch;
% Allowed plot types: 'full', 'lj', 'full-derivative', 'lj-derivative',
% 'dispersion', 'dispersion-derivative', 'repulsive',
% 'repulsive-derivative'

if p.Results.MDP_Minimize
    MDP = 'MinMDP';
else
    MDP = 'MDP';
end

tol_steep = 1e4; % If steepness is higher than this tolerance, then perfect hard spheres are used

%% Gaussian adjustments
G_a.MM = Settings.GAdjust_MM(:,1);
G_b.MM = Settings.GAdjust_MM(:,2);
G_c.MM = Settings.GAdjust_MM(:,3);

G_a.XX = Settings.GAdjust_XX(:,1);
G_b.XX = Settings.GAdjust_XX(:,2);
G_c.XX = Settings.GAdjust_XX(:,3);

G_a.MX = Settings.GAdjust_MX(:,1);
G_b.MX = Settings.GAdjust_MX(:,2);
G_c.MX = Settings.GAdjust_MX(:,3);

%% Conversion factors and fundamental constants
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1

%% Generate range (r) in nm
U.r = Startpoint:Settings.Table_StepSize:Settings.Table_Length;

%% Calculate Repulsive Wall Height Z: MX = +-   MM = ++     XX = --
Z.MX = Settings.S.Z.All*Settings.S.Z.MX; % kJ/mol
Z.MM = Settings.S.Z.All*Settings.S.Z.MM;
Z.XX = Settings.S.Z.All*Settings.S.Z.XX;

%% Calculate Repulsive Steepness Parameter B: MX = +-   MM = ++     XX = --
B.MX = Settings.S.b.All*Settings.S.b.MX; % nm^-1
B.MM = Settings.S.b.All*Settings.S.b.MM;
B.XX = Settings.S.b.All*Settings.S.b.XX;

%% Calculate Hard Sphere Distance Parameter R: MX = +-   MM = ++     XX = --
R.MX = Settings.S.r_d.All*Settings.S.r_d.MX; % nm
R.MM = Settings.S.r_d.All*Settings.S.r_d.MM;
R.XX = Settings.S.r_d.All*Settings.S.r_d.XX;

%% Coulomb Charges Q: MX = +-   MM = ++     XX = --
Q.MX = -Settings.S.Q;
Q.MM = Settings.S.Q;
Q.XX = Settings.S.Q;

%% Build Interactions
for interaction = {'MX' 'XX' 'MM'}
    int = interaction{1};
    
    [U.(int).f0,U.(int).df0] = Coulomb_Potential(Settings,U.r,int);
    
    Coulomb  =  k_0*(e_c^2)*Q.(int).*U.(int).f0;
    dCoulomb = -k_0*(e_c^2)*Q.(int).*U.(int).df0;
    
    %% If Damping at close range, affects only attractive coulombic interactions
    if strcmp(int,'MX') && ( Settings.CR_Damp.(int).r_d >= 0 && Settings.CR_Damp.(int).b >= 0 )
        r_d = Settings.CR_Damp.(int).r_d;
        sb  = Settings.CR_Damp.(int).b;
        
        f_r = 1./(1 + exp(-sb.*(U.r - r_d))); % sigmoid damping function
        df_r = (sb.*exp(-sb.*(U.r - r_d)))./((1 + exp(-sb.*(U.r - r_d))).^2); % sigmoid damping function derivative
        f_cutoff = 1/(1 + exp(-sb*(Settings.(MDP).RVDW_Cutoff - r_d))); % damping function value at vdw
        
        % Close Range Attractive Damping
        Coul_tot  =   f_r.*Coulomb;
        Coul_dtot =   f_r.*dCoulomb + df_r.*Coulomb;
        Coul_h    = - Coulomb  + f_r.*Coulomb;
        Coul_dh   = - dCoulomb + f_r.*dCoulomb + df_r.*Coulomb;
        h_cutoff  = - k_0*(e_c^2)*Q.(int)/(Settings.(MDP).RVDW_Cutoff) ...
           + f_cutoff*k_0*(e_c^2)*Q.(int)/(Settings.(MDP).RVDW_Cutoff);
        
    else
    %% No damping
        Coul_tot  = Coulomb;
        Coul_dtot = dCoulomb;
        Coul_h    = 0;
        Coul_dh   = 0;
        h_cutoff  = 0;
    end
    
    %% Optional: modify potential with Gaussian Adjustments
    G_r = zeros(1,length(U.r));
    dG_r = zeros(1,length(U.r));
    G_r_Cutoff = 0;
    for i = 1:length(G_a.(int))
        if G_a.(int)(i) == 0
            G_r = zeros(1,length(U.r));
            G_r_Cutoff = 0;
            dG_r = zeros(1,length(U.r));
        else
            G_r = G_r + G_a.(int)(i).*exp((-(U.r - G_b.(int)(i)).^2)./(2.*(G_c.(int)(i).^2)));
            G_r_Cutoff = G_r_Cutoff + G_a.(int)(i)*exp((-(Settings.(MDP).RVDW_Cutoff - G_b.(int)(i))^2)/(2*(G_c.(int)(i)^2)));
            dG_r = dG_r + (G_a.(int)(i).*(U.r - G_b.(int)(i))).*(exp((-(U.r - G_b.(int)(i)).^2)./(2.*(G_c.(int)(i).^2))))/(G_c.(int)(i).^2);
        end
    end
    
    %% Build repulsive walls
    if B.(int) < tol_steep
        g_r = Z.(int)./(1 + exp(B.(int).*(U.r - R.(int)))); % Hard wall
        dg_r = - Z.(int).*(B.(int).*exp(B.(int).*(U.r - R.(int))))./((1 + exp(B.(int).*(U.r - R.(int)))).^2); % Hard wall derivative
        g_cutoff = Z.(int)./(1 + exp(B.(int)*(Settings.(MDP).RVDW_Cutoff - R.(int)))); % Hard wall function value at vdw cutoff
    else
        g_r = zeros(size(U.r));
        dg_r = zeros(size(U.r));
        g_r(U.r <= R.(int)) = Z.(int);
        
        if Settings.(MDP).RVDW_Cutoff <= R.(int)
            g_cutoff  = Z.(int);
        else
            g_cutoff = 0;
        end
    end
    
    %% Build PES
    
    % Total potential
    U.(int).Total = Coul_tot + g_r + G_r;
    
    % Potential components

    U.(int).g = G_r;                  % Dispersion g(r) (includes any Gaussian augments)
    U.(int).h = g_r + Coul_h;         % Short range repulsion h(r) (includes coulomb damping)    

    % Total derivative
    U.(int).dTotal = Coul_dtot + dg_r - dG_r;

    % (Negative of) Derivative components
    U.(int).dg =   dG_r;                     % Dispersion -dg(r)/dr
    U.(int).dh = - dg_r - Coul_dh;           % Short range repulsion -dh(r)/dr
    
    if contains(Settings.(MDP).vdw_modifier,'potential-shift','IgnoreCase',true)
        EVDW_Cutoff = g_cutoff + h_cutoff + G_r_Cutoff;

        % Shift by the dispersion energy at vdw cutoff radius. only affects one
        % energy component, not derivatives (i.e. forces)
        U.(int).Total = U.(int).Total - EVDW_Cutoff;
        U.(int).g     = U.(int).g     - EVDW_Cutoff;
    end

    % remove infinities
    U.(int) = Remove_Infinities(U.(int));
    
    % Print
    U.(int).f  = k_0*(e_c^2)*q.(int(1))*q.(int(2)).*U.(int).f0;
    U.(int).df = k_0*(e_c^2)*q.(int(1))*q.(int(2)).*U.(int).df0;
    U.(int).Total = U.(int).f + U.(int).g + U.(int).h;
    U.(int).dTotal = -(U.(int).df + U.(int).dg + U.(int).dh);
end

%% PLOT if plotswitch chosen
if Plotswitch
    figure;
    % Options
    lw=3;
    fs=35;
    
    h = cell(1,8);
    hold on
    switch lower(PlotType)
        case 'full'
            h{1} = plot(U.r.*10,k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.f + U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.f + U.MM.g + U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.f + U.XX.g + U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Full Potential';
        case 'full-derivative'
            h{1} = plot(U.r.*10,k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.df + U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.df + U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.df + U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Derivative of Full Potential';
        case 'lj'
            h{1} = plot(U.r.*10,U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.g + U.MM.h,'Color','g','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.g + U.XX.h,'Color','b','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Lennard-Jones Potential';
        case 'lj-derivative'
            h{1} = plot(U.r.*10,U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Lennard-Jones Potential';
        case 'dispersion'
            h{1} = plot(U.r.*10,U.MX.g,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.g,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.g,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Dispersion Potential';
        case 'dispersion-derivative'
            h{1} = plot(U.r.*10,U.MX.dg,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.dg,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.dg,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Dispersion Potential';
        case 'repulsive'
            h{1} = plot(U.r.*10,U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Repulsive Potential';
        case 'repulsive-derivative'
            h{1} = plot(U.r.*10,U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Repulsive Potential';
    end
    
    title(['Plot of ' ttxt ' for Charged HS Model'],...
       'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('Separation [\AA]','fontsize',fs,'Interpreter','latex');
    ylabel('Potential Energy [kJ mol$^{-1}$]','fontsize',fs,'Interpreter','latex');

    ylim(yl);
    xlim([Startpoint Settings.Table_Length]);

    % Blank line
    hline = refline([0 0]);
    hline.Color = 'k';
    hline.LineWidth = lw-1;
    hline.LineStyle = '--';
    leg1 = legend([h{:}],{'HS - MX' 'HS - MM' 'HS - XX'});
    leg1.Interpreter = 'latex';
end

end