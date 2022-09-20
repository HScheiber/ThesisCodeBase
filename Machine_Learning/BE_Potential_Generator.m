function [U_MX_out, U_MM_out, U_XX_out] = BE_Potential_Generator(Settings,varargin)

% Optional inputs
p = inputParser;
p.FunctionName = 'BD_Potential_Generator';
addOptional(p,'PlotType','full')
addOptional(p,'ReturnAsStructure',false);
addOptional(p,'Startpoint',0);
addOptional(p,'Plotswitch',false);
addOptional(p,'MDP_Minimize',false);
parse(p,varargin{:});
PlotType = p.Results.PlotType;
ReturnAsStructure = p.Results.ReturnAsStructure;
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

%% Conversion factors and fundamental constants
kj_per_erg = 1e-10; % kJ per erg
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
nm_per_Ang = 0.1; % nm per Angstrom
kJ_per_kcal = 4.184; % kj per kcal

[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);

%% JC Ion sigma/epsilon Parameters in SPC/E water
Param.Li.sigma = (0.791*nm_per_Ang)*(2^(5/6)); % nm
Param.Li.epsilon = (0.3367344)*kJ_per_kcal; % kJ/mol

Param.Na.sigma = (1.212*nm_per_Ang)*(2^(5/6)); % nm
Param.Na.epsilon = (0.3526418)*kJ_per_kcal; % kJ/mol

Param.K.sigma = (1.593*nm_per_Ang)*(2^(5/6)); % nm
Param.K.epsilon = (0.4297054)*kJ_per_kcal; % kJ/mol

Param.Rb.sigma = (1.737*nm_per_Ang)*(2^(5/6)); % nm
Param.Rb.epsilon = (0.4451036)*kJ_per_kcal; % kJ/mol

Param.Cs.sigma = (2.021*nm_per_Ang)*(2^(5/6)); % nm
Param.Cs.epsilon = (0.0898565)*kJ_per_kcal; % kJ/mol

Param.F.sigma = (2.257*nm_per_Ang)*(2^(5/6)); % nm
Param.F.epsilon = (0.0074005)*kJ_per_kcal; % kJ/mol

Param.Cl.sigma = (2.711*nm_per_Ang)*(2^(5/6)); % nm
Param.Cl.epsilon = (0.0127850)*kJ_per_kcal; % kJ/mol

Param.Br.sigma = (2.751*nm_per_Ang)*(2^(5/6)); % nm
Param.Br.epsilon = (0.0269586)*kJ_per_kcal; % kJ/mol

Param.I.sigma = (2.919*nm_per_Ang)*(2^(5/6)); % nm
Param.I.epsilon = (0.0427845)*kJ_per_kcal; % kJ/mol

%% Parameter: q (charge)
q.Li =  Settings.S.Q; % atomic
q.Na =  Settings.S.Q; % atomic
q.K  =  Settings.S.Q; % atomic
q.Rb =  Settings.S.Q; % atomic
q.Cs =  Settings.S.Q; % atomic

q.F  = -Settings.S.Q; % atomic
q.Cl = -Settings.S.Q; % atomic
q.Br = -Settings.S.Q; % atomic
q.I  = -Settings.S.Q; % atomic

%% TF Repulsive Size Parameter sigma (AKA r+/-): P = +   M = -
% Metals
sigma.Li = 0.816*nm_per_Ang; % nm
sigma.Na = 1.170*nm_per_Ang; % nm
sigma.K  = 1.463*nm_per_Ang; % nm
sigma.Rb = 1.587*nm_per_Ang; % nm
sigma.Cs = 1.720*nm_per_Ang; % nm

% Halides
sigma.F  = 1.179*nm_per_Ang; % nm
sigma.Cl = 1.585*nm_per_Ang; % nm
sigma.Br = 1.716*nm_per_Ang; % nm
sigma.I  = 1.907*nm_per_Ang; % nm

%% TF Parameter: Number of Valence electrons (for Pauling Coefficient Calculation)
% Metals
valence.Li = 2;
valence.Na = 8;
valence.K = 8;
valence.Rb = 8;
valence.Cs = 8;

% Halides
valence.F = 8;
valence.Cl = 8;
valence.Br = 8;
valence.I = 8;

%% TF Hardness Parameter Rho
rho.LiF = 0.299*nm_per_Ang; % nm
rho.LiCl = 0.342*nm_per_Ang; % nm
rho.LiBr = 0.353*nm_per_Ang; % nm
rho.LiI = 0.430*nm_per_Ang; % nm

rho.NaF = 0.330*nm_per_Ang; % nm
rho.NaCl = 0.317*nm_per_Ang; % nm
rho.NaBr = 0.340*nm_per_Ang; % nm
rho.NaI = 0.386*nm_per_Ang; % nm

rho.KF = 0.338*nm_per_Ang; % nm
rho.KCl = 0.337*nm_per_Ang; % nm
rho.KBr = 0.335*nm_per_Ang; % nm
rho.KI = 0.355*nm_per_Ang; % nm

rho.RbF = 0.328*nm_per_Ang; % nm
rho.RbCl = 0.318*nm_per_Ang; % nm
rho.RbBr = 0.335*nm_per_Ang; % nm
rho.RbI = 0.337*nm_per_Ang; % nm

rho.CsF = 0.282*nm_per_Ang; % nm
rho.CsCl = 0.272*nm_per_Ang; % nm THIS IS ESTIMATED
rho.CsBr = 0.289*nm_per_Ang; % nm THIS IS ESTIMATED
rho.CsI = 0.291*nm_per_Ang; % nm THIS IS ESTIMATED

%% Huggins-Mayer potential parameter b (same for all salts)
b = (0.338e-12)*kj_per_erg*NA; % kJ/mol

%% Calculate Pauling Coefficients beta: MX = +-   MM = ++     XX = --
beta.MM = 1 + 2/valence.(Metal);   % Unitless
beta.XX = 1 - 2/valence.(Halide); % Unitless
beta.MX = 1 + 1/valence.(Metal) - 1/valence.(Halide); % Unitless

%% Calculate TF Repulsive Exponential Parameter alpha: MX = +-   MM = ++     XX = --
alpha.MM = Settings.S.A.All*Settings.S.A.MM/rho.(Settings.Salt); % nm^-1
alpha.MX = Settings.S.A.All*Settings.S.A.MX/rho.(Settings.Salt); % nm^-1
alpha.XX = Settings.S.A.All*Settings.S.A.XX/rho.(Settings.Salt); % nm^-1

%% Calculate TF Repulsive Scaling Parameter B: MX = +-   MM = ++     XX = -- (Including scaling)
B.MM = Settings.S.R.All*Settings.S.R.MM*beta.MM*b*exp(2*sigma.(Metal)/rho.(Settings.Salt));
B.XX = Settings.S.R.All*Settings.S.R.XX*beta.XX*b*exp(2*sigma.(Halide)/rho.(Settings.Salt));
B.MX = Settings.S.R.All*Settings.S.R.MX*beta.MX*b*exp((sigma.(Metal) + sigma.(Halide))/rho.(Settings.Salt));

%% Calculate parameters of interest for LJ potential
sigma_MM = Param.(Metal).sigma;
sigma_XX = Param.(Halide).sigma;

epsilon_MM = Param.(Metal).epsilon;
epsilon_XX = Param.(Halide).epsilon;

% Change parameteters into C6/r6 format and apply mixing rules
CMM_pre = 4*epsilon_MM*(sigma_MM^6);
CXX_pre = 4*epsilon_XX*(sigma_XX^6);
C.MM = Settings.S.D.All*Settings.S.D.MM*CMM_pre;
C.XX = Settings.S.D.All*Settings.S.D.XX*CXX_pre;
C.MX = Settings.S.D.All*Settings.S.D.MX*sqrt(CMM_pre*CXX_pre);

%% Generate range (r) in nm
r = Startpoint:Settings.Table_StepSize:Settings.Table_Length;

%% If Damping at close range, affects all attractive interactions
for interaction = {'MX' 'XX' 'MM'}
    int = interaction{1};
    if strcmp(int,'MM')
        Y1 = Metal;
        Y2 = Metal;
    elseif strcmp(int,'XX')
        Y1 = Halide;
        Y2 = Halide;
    elseif strcmp(int,'MX')
        Y1 = Metal;
        Y2 = Halide;
    end
    
    %% Build PES
    
    % Components of potential
    U.(int).f = 1./r; % Electrostatics function f(r)
    U.(int).g = - C.(int)./(r.^6); % Dispersion g(r)
    U.(int).h = B.(int)*exp(-alpha.(int).*r); % Short range repulsion (with possible close-range coulomb damping)
    
    % Negative components of derivative
    U.(int).df = 1./(r.^2); % Electrostatics function (not including Coulomb constant or charges)
    U.(int).dg = -C.(int).*6./(r.^7); % Dispersion -dg(r)/dr
    U.(int).dh = alpha.(int)*B.(int)*exp(-alpha.(int).*r); % Short range repulsion
    
    % Shift the potential to zero at the cutoff
    if contains(Settings.(MDP).vdw_modifier,'potential-shift','IgnoreCase',true)
        EVDW_Cutoff = B.(int)*exp(-alpha.(int)*Settings.(MDP).RVDW_Cutoff) ...
                      -C.(int)./(Settings.(MDP).RVDW_Cutoff.^6);

        % Shift by the dispersion energy at vdw cutoff radius. only affects one
        % energy component, not derivatives (i.e. forces)
        U.(int).g = U.(int).g - EVDW_Cutoff;
    end
    
    %% Grab the inflection point of the potential
    U_Total  =  k_0*(e_c^2).*q.(Y1)*q.(Y2).*U.(int).f  + U.(int).g  + U.(int).h; % Full potential
    U_Coul   =  k_0*(e_c^2).*q.(Y1)*q.(Y2).*U.(int).f; % Coulomb potential
    dU_Total = -k_0*(e_c^2).*q.(Y1)*q.(Y2).*U.(int).df - U.(int).dg - U.(int).dh; % Full potential derivative
    dU_Coulomb = -k_0*(e_c^2).*q.(Y1)*q.(Y2).*U.(int).df;
    
    peaks_idx = [false islocalmax(U_Total(2:end),'MinProminence',1e-8)];
    peak_r = r(peaks_idx);
    if numel(peak_r) > 1
        peak_r = peak_r(1);
    end
    
    inflex_idx = find([false islocalmax(dU_Total(2:end),'MinProminence',1e-8) | islocalmin(dU_Total(2:end),'MinProminence',1e-8)]);
    if ~isempty(peak_r) && ~isempty(inflex_idx) % Ensure peak exists
        inflex_idx(r(inflex_idx) <= peak_r) = [];
        inflex_idx = inflex_idx(1);
        inflex_r = r(inflex_idx);
        dU_infl = dU_Total(inflex_idx);
        
        % Generate a steep repulsion beyond the inflection
        below_infl_idx = (r < inflex_r); % values of r below the inflection point
        r_g = r(below_infl_idx); % nm
        D = -dU_infl*(1/alpha.(int))*exp(alpha.(int)*inflex_r); % Wall prefactor chosen to match the total derivative at the inflection point
        
        fwall = D.*exp(-alpha.(int).*r_g) - D.*exp(-alpha.(int).*inflex_r) ...
            + U_Total(inflex_idx) - U_Coul(below_infl_idx); % Repulsive wall shifted to its value at the inflection point
        dfwall = alpha.(int)*D.*exp(-alpha.(int).*r_g) + dU_Coulomb(below_infl_idx); % Wall -derivative

        % Kill the attractive interaction beyond the peak
        %g_at_infl = U_LJ(inflex_idx);
        N_belowpeak = inflex_idx-1;
        U.(int).g(below_infl_idx) = zeros(1,N_belowpeak);
        U.(int).dg(below_infl_idx) = zeros(1,N_belowpeak);

        % Remove infinity at 0
        fwall(1) = fwall(2)*2;
        dfwall(1) = dfwall(2);

        % Add this repulsion to the repulsive part of the function
        U.(int).h(below_infl_idx) = fwall;
        U.(int).dh(below_infl_idx) = dfwall;
    end
        
%         %% Visualization
%         nm_per_m = 1e+9; % nm per m
%         NA = 6.0221409e23; % Molecules per mole
%         e_c = 1.60217662e-19; % Elementary charge in Coulombs
%         epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
%         k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
%         
%         if idx == 1
%             U.Total = -k_0*(e_c^2).*(Settings.S.Q^2).*U.f0 + U.h + U.g ;
%         else
%             U.Total =  k_0*(e_c^2).*(Settings.S.Q^2).*U.f0 + U.h + U.g ;
%         end
%         hold on
%         plot(U.r(2:end).*10,U.Total(2:end),'-k','Linewidth',3)
%         scatter(inflex_r.*10,U.Total(U.r == inflex_r),100,'r','Linewidth',3,'MarkerEdgeColor','r')
%         %%
        
%         %% Testing visualization
%         if idx == 1
%             U.Total = -k_0*(e_c^2).*(Settings.S.Q^2).*U.f0 + U.h + U.g ;
%         else
%             U.Total =  k_0*(e_c^2).*(Settings.S.Q^2).*U.f0 + U.h + U.g ;
%         end
%         plot(U.r.*10,U.Total,':r','Linewidth',3)
%         ylim([-1000 4000])
%         xlim([0 5])
%         set(gca,'Fontsize',32,'TickLabelInterpreter','latex','XTick',0:1:5)
%         xlabel(gca,'$r_{ij}$ [\AA]','fontsize',32,'interpreter','latex')
%         ylabel(gca,'$u_{ij}$ [kJ mol$^{-1}$]','fontsize',32,'interpreter','latex')
%         exportgraphics(gca,'Augmented_Potential.eps')
%         %%

    % remove infinities
    U.(int) = Remove_Infinities(U.(int));
    
    % Print
    U_out = [r ; U.(int).f ; U.(int).df ; U.(int).g ; U.(int).dg ; U.(int).h ; U.(int).dh];
    U.(int).out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],U_out(:)) );
end

if ReturnAsStructure
	U.MX.f0 = U.MX.f;
    U.MM.f0 = U.MM.f;
    U.XX.f0 = U.XX.f;
    
	U.MX.df0 = U.MX.df;
    U.MM.df0 = U.MM.df;
    U.XX.df0 = U.XX.df;
    
	U.MX.f = k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.f;
    U.MM.f = k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.f;
    U.XX.f = k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.f;
    
	U.MX.df = k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.df;
    U.MM.df = k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.df;
    U.XX.df = k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.df;
    
    U.MX.Total = U.MX.f + U.MX.g + U.MX.h;
    U.MM.Total = U.MM.f + U.MM.g + U.MM.h;
    U.XX.Total = U.XX.f + U.XX.g + U.XX.h;
    
    U.MX.dTotal = -(U.MX.df + U.MX.dg + U.MX.dh);
    U.MM.dTotal = -(U.MM.df + U.MM.dg + U.MM.dh);
    U.XX.dTotal = -(U.XX.df + U.XX.dg + U.XX.dh);
    
    U_MX_out = U.MX;
    U_MM_out = U.MM;
    U_XX_out = U.XX;

    U_MX_out.r = r;
    U_MM_out.r = r;
    U_XX_out.r = r;
else
    U_MX_out = U.MX.out;
    U_MM_out = U.MM.out;
    U_XX_out = U.XX.out;
end

%% PLOT if plotswitch chosen
if Plotswitch
    figure;
    % Options
    lw=2;
    fs=25;

    h = cell(1,8);
    hold on
    switch lower(PlotType)
        case 'full'
            h{1} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.f + U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.f + U.MM.g + U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.f + U.XX.g + U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Full Potential';
        case 'full-derivative'
            h{1} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.df + U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.df + U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.df + U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Derivative of Full Potential';
        case 'lj'
            h{1} = plot(r.*10,U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.g + U.MM.h,'Color','g','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.g + U.XX.h,'Color','b','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Lennard-Jones Potential';
        case 'lj-derivative'
            h{1} = plot(r.*10,U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Lennard-Jones Potential';
        case 'dispersion'
            h{1} = plot(r.*10,U.MX.g,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.g,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.g,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Dispersion Potential';
        case 'dispersion-derivative'
            h{1} = plot(r.*10,U.MX.dg,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dg,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dg,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Dispersion Potential';
        case 'repulsive'
            h{1} = plot(r.*10,U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Repulsive Potential';
        case 'repulsive-derivative'
            h{1} = plot(r.*10,U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Repulsive Potential';
    end
    
    title(['Plot of ' ttxt ' for ' Settings.Salt ' BH Model'],...
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
    leg1 = legend([h{:}],{[Metal '$^{+}$' ' - ' Halide '$^{-}$'] [Metal '$^{+}$' ' - ' Metal '$^{+}$'] [Halide '$^{-}$' ' - ' Halide '$^{-}$']});
    leg1.Interpreter = 'latex';
end

end