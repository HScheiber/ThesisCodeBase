% Jung-Cheatham model parameters adapted for three water models with
% Lorenz-Berthelot combining rules
% Use Watermodel = 'SPC/E', 'TIP3P', or 'TIP4PEW'
% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETER AND
% kJ/mol

% C6Damping:
% 0 = no (default) damping. This is default of JC model.
% 1 = BJ/rational damping (same as in D3(BJ))
% 2 = Tang Damping (Essentially removes dispersion in JC)
% 3 = MMDRE Damping function (Very weak damping, damps mainly at mid range)
% 4 = PAMoC Damping function (fairly weak damping, damps mainly at mid range)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)

% GAdjust are N x 3 arrays of gaussian parameters
% (i , 1) is the Gaussian height of the ith adjustment (may be negative or
% positive)
% (i , 2) is the center point of the ith Gaussian (should be positive)
% (i , 3) is the standard deviation or width (negative and positive values
% are the same)

% CRDamping = Close Range Damping
function U_AB_out = ClayFF_Potential_Generator(Startpoint,Endpoint,Spacing,...
    Atom_A,Atom_B,plotswitch,S,vdw_modifier,RVDW_Cutoff,C6Damp,CRDamp,GAdjust_AB)

%% Gaussian adjustments
if exist('GAdjust_AB','var')
    G_a_AB = GAdjust_AB(:,1);
    G_b_AB = GAdjust_AB(:,2);
    G_c_AB = GAdjust_AB(:,3);
else
    G_a_AB = 0;
    G_b_AB = 0;
    G_c_AB = 1;
end

%% Conversion factors and fundamental constants
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
nm_per_Ang = 0.1; % nm per Angstrom
kJ_per_kcal = 4.184; % kj per kcal
R0_to_Sigma = 2^(-1/6); % Conversion factor from R0 to sigma

%% JC Ion Parameters in SPC/E water
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

% ClayFF Parameters
Param.Al.sigma = 4.7943*nm_per_Ang*R0_to_Sigma; % nm
Param.Al.epsilon = 1.3298e-6*kJ_per_kcal; % kJ/mol

Param.O.sigma = 3.5532*nm_per_Ang*R0_to_Sigma; % nm
Param.O.epsilon = 0.1554*kJ_per_kcal; % kJ/mol

%% Parameter: q (charge)
q.Li =  S.Q; % atomic
q.Na =  S.Q; % atomic
q.K  =  S.Q; % atomic
q.Rb =  S.Q; % atomic
q.Cs =  S.Q; % atomic

q.F  = -S.Q; % atomic
q.Cl = -S.Q; % atomic
q.Br = -S.Q; % atomic
q.I  = -S.Q; % atomic

q.Al =  S.Q*1.5750;
q.O  = -S.Q*1.0500;

%% Calculate parameters of interest for LJ potential
sigma_AB = S.S.All*S.S.AB*( Param.(Atom_A).sigma + Param.(Atom_B).sigma )/2;
epsilon_AB = S.E.All*S.E.AB*sqrt(Param.(Atom_A).epsilon*Param.(Atom_B).epsilon);

% Change parameteters into A/r12 - B/r6 format
A.AB = S.R.All*S.R.AB*4*epsilon_AB*(sigma_AB^12);
B.AB = S.D.All*S.D.AB*4*epsilon_AB*(sigma_AB^6);

%% Generate range (r) in nm
r = Startpoint:Spacing:Endpoint;

%% If Damping at close range, affects all attractive interactions
if CRDamp.AB.r_d >= 0 && CRDamp.AB.b >= 0
    r_d = CRDamp.AB.r_d;
    b  = CRDamp.AB.b; 

    f_r.AB = 1./(1 + exp(-b.*(r - r_d))); % sigmoid damping function
    df_r.AB = (b.*exp(-b.*(r - r_d)))./((1 + exp(-b.*(r - r_d))).^2); % sigmoid damping function derivative
    f_cutoff.AB = 1./(1 + exp(-b.*(RVDW_Cutoff - r_d)));  % value of f at vdw cutoff
else
    f_r.AB = 1; % No damping
    df_r.AB = 0; % Damping derivative is zero
    f_cutoff.AB = 1; % value of f at vdw cutoff
end

% Apply C6 Damping
%% No Damping
if C6Damp.AB == 0
    % No-damping function
    f6.AB = f_r.AB;

    % Derivative parts of no-damping function   
    df6.AB = df_r.AB;

    % Values at vdw cutoff
    f6_cutoff.AB = f_cutoff.AB;

%% BJ-Type Rational Damping
elseif C6Damp.AB == 1

    % Generate conversion factors
    Bohr_nm = 0.0529177; % a_0 - > Angstrom
    c6conv = 1e-3/2625.4999/((0.052917726)^6); % J/mol nm^6 - > au (from sourcecode)
    J_kJ = 1e-3; % J - > kJ
    Ha_kJmol = 2625.4999; % Ha - > kJ/mol
    c6units = (1/c6conv)*J_kJ; % au - > kJ/mol nm^6
    c8units = (Ha_kJmol)*(Bohr_nm^8); % au - > kJ/mol nm^8

    % Factor used to calculate C8
    sqrt_Q.Li = 5.019869340000000;
    sqrt_Q.Na = 6.585855360000000;
    sqrt_Q.K  = 7.977627530000000;
    sqrt_Q.Rb = 9.554616980000000;
    sqrt_Q.Cs = 11.02204549000000;
    sqrt_Q.F  = 2.388252500000000;
    sqrt_Q.Cl = 3.729323560000000;
    sqrt_Q.Br = 4.590896470000000;
    sqrt_Q.I  = 5.533218150000000;

    % Calculate C8 (needed for cutoff radius)
    C8.AB = 3.0*(B.AB/c6units)*sqrt_Q.(Atom_A)*sqrt_Q.(Atom_B)*c8units; % in kJ/mol nm^8

    % Damping distances (no C8 term so define wrt crystal radii)
    R0.AB = sqrt(C8.AB/B.AB); % in Angstroms

    % Damping functions (unitless)
    f6.AB = f_r.AB./( 1 + ( R0.AB ./ r ).^6 );

    % Values of damping function at vdw cutoff
    f6_cutoff.AB = f_cutoff.AB/( 1 + ( R0.AB / RVDW_Cutoff )^6 );

    % Derivative of damping functions
    df6.AB = f_r.AB.*(6.*(R0.AB.^6).*(r.^5)./(((r.^6) + (R0.AB.^6)).^2)) + df_r.AB./( 1 + ( R0.AB ./ r ).^6 );

%% Tang and Toennies Damping function. Cite:
% "An improved simple model for the van der Waals potential based on universal damping functions for the dispersion coefficients."
% K. T. Tang, J. P. Toennies
% J. Chem. Phys. 1984, 80, 3726-3741.
elseif C6Damp.AB == 2

    % use the TF hardness parameters
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
    alpha = rho.(Salt);

    % C6 damping functions
    f6sum = 0;
    for k = 0:6
        f6sum = f6sum + ((alpha.*r).^k)./factorial(k); 
    end
    f6.AB = f_r.AB.*(1 - f6sum.*exp(-alpha.*r));

    % Calculate C6 damping derivatives
    df6sum = 0;
    for k = 1:6
        df6sum = df6sum + k.*(alpha.^k).*(r.^(k-1))./factorial(k);
    end
    df6.AB = f_r.AB.*((alpha.*exp(-alpha.*r)).*f6sum ...
        - (exp(-alpha.*r)).*df6sum) + df_r.AB.*(1 - f6sum.*exp(-alpha.*r));

    % Values of damping function at vdw cutoff
    f6_cutoff.AB = f_cutoff.AB*(1 - sum(((alpha.*RVDW_Cutoff).^(0:6))./factorial(0:6)).*exp(-alpha.*RVDW_Cutoff));

%% Family of other damping functions related by a single formula, requiring van der waals radii
elseif C6Damp.AB >= 3 && C6Damp.AB <= 6

    % Define crystal radii source:
    % https://en.wikipedia.org/wiki/Ionic_radius
    R0.Li = 0.09; % nm
    R0.Na = 0.116; % nm
    R0.K  = 0.152; % nm
    R0.Rb = 0.166; % nm
    R0.Cs = 0.181; % nm
    R0.F  = 0.119; % nm
    R0.Cl = 0.167; % nm
    R0.Br = 0.182; % nm
    R0.I  = 0.206; % nm

    %% MMDRE Damping function. Citation:
    % "Transferable ab initio intermolecular potentials. 1. Derivation from methanol dimer and trimer calculations"
    % W.T.M. Mooij, F.B. van Duijneveldt, J.G.C.M. van Duijneveldt-van de Rijdt, B.P. van Eijck
    % J. Phys. Chem. A 1999, 103, 9872-9882.
    if C6Damp.AB == 3
        a = -1;
        b = 7.19;
        m = 3;
        n = 2;

    %% PAMoC Damping function. Citation:
    % "Empirical correction to density functional theory for van der Waals interactions"
    % Q. Wu, W. Yang
    % J. Chem. Phys. 2002, 116, 515-524.
    elseif C6Damp.AB == 4
        a = -1;
        b = 3.54;
        m = 3;
        n = 2;

    %% EHFSK Damping function. Citation:
    % "Hydrogen bonding and stacking interactions of nucleic acid base pairs: A density-functional-theory based treatment"
    % M. Elstner, P. Hobza, T. Frauenheim, S. Suhai, E. Kaxiras
    % J. Chem. Phys. 2001, 114, 5149-5155.
    elseif C6Damp.AB == 5
        a = -1;
        b = 3;
        m = 7;
        n = 4;

    %% WY damping function. Citation:
    % "Empirical correction to density functional theory for van der Waals interactions"
    % Q. Wu, W. Yang
    % J. Chem. Phys. 2002, 116, 515-524.
    elseif C6Damp.AB == 6
        a = exp(23);
        b = 23;
        m = 1;
        n = -1;
    end

    R0.AB = R0.(Atom_A) + R0.(Atom_B);

    % Damping functions
    f6.AB = f_r.AB.*((1 + a.*exp(-b.*(r./R0.AB).^m) ).^n);

    % Derivative parts of damping functions
    df6.AB = f_r.AB.*(-(a.*b.*m.*n.*(r./R0.AB).^(m - 1).*exp(-b.*(r./R0.AB).^m).*(a.*exp(-b.*(r/R0.AB).^m) + 1).^(n - 1))./R0.AB) ...
        + df_r.AB.*((1 + a.*exp(-b.*(r./R0.AB).^m) ).^n);

    % Values at vdw cutoff
    f6_cutoff.AB = f_cutoff.AB*((1 + a*exp(-b*(RVDW_Cutoff/R0.AB)^m) )^n);
end

%% Modify potential with Gaussian Adjustments
G_r_AB = zeros(1,length(r));
dG_r_AB = zeros(1,length(r));
G_r_AB_Cutoff = 0;
for i = 1:length(G_a_AB)
    G_r_AB = G_r_AB + G_a_AB(i).*exp((-(r - G_b_AB(i)).^2)./(2.*(G_c_AB(i).^2)));
    G_r_AB_Cutoff = G_r_AB_Cutoff + G_a_AB(i)*exp((-(RVDW_Cutoff - G_b_AB(i))^2)/(2*(G_c_AB(i)^2)));
    dG_r_AB = dG_r_AB - (G_a_AB(i).*(r - G_b_AB(i))).*(exp((-(r - G_b_AB(i)).^2)./(2.*(G_c_AB(i).^2))))/(G_c_AB(i).^2);
end

%% Build PES: Opposite charges
if q.(Atom_A)*q.(Atom_B) < 0
    % Plus - Minus total potential
    U_AB.Total = f_r.AB.*k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)./(r) ...
        + A.AB./(r.^12) ...
        - f6.AB.*B.AB./(r.^6) ...
        + G_r_AB;

    % components
    U_AB.f = 1./r; % Electrostatics function f(r)
    U_AB.g = -f6.AB.*B.AB./(r.^6) + G_r_AB; % Dispersion g(r)
    U_AB.h = A.AB./(r.^12) - k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)./(r) + ...
        f_r.AB.*k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)./(r); % Short range repulsion

    % Plus - Minus total derivative
    U_AB.dTotal = -f_r.AB.*k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)./(r.^2) ...
        + df_r.AB.*k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)./(r)...
        - A.AB.*12./(r.^13) ...
        + f6.AB.*B.AB.*6./(r.^7) ...
        - df6.AB.*B.AB./(r.^6) ...
        + dG_r_AB;

    % components
    U_AB.df = 1./(r.^2);% Electrostatics function
    U_AB.dg = - f6.AB.*B.AB.*6./(r.^7) + df6.AB.*B.AB./(r.^6) - dG_r_AB; % Dispersion -dg(r)/dr
    U_AB.dh = + A.AB.*12./(r.^13) - k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)./(r.^2) ...
        + f_r.AB.*k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)./(r.^2) ...
        - df_r.AB.*k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)./(r);% Short range repulsion

    if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
        EVDW_Cutoff = A.AB./(RVDW_Cutoff.^12) ...
            - k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)/(RVDW_Cutoff) ...
            + f_cutoff.AB*k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)/(RVDW_Cutoff) ...
            - f6_cutoff.AB.*B.AB./(RVDW_Cutoff.^6) ...
            + G_r_AB_Cutoff;

        % Shift by the dispersion energy at vdw cutoff radius. only affects one
        % energy component, not derivatives (i.e. forces)
        U_AB.Total = U_AB.Total - EVDW_Cutoff;
        U_AB.g = U_AB.g - EVDW_Cutoff;
    end

    % remove infinities
    U_AB = Remove_Infinities(U_AB);

%% Build PES: Same charge
else
    % Plus - Plus total potential
    U_AB.Total = k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)./(r) ...
        + A.AB./(r.^12) ...
        - f6.AB.*B.AB./(r.^6) ...
        + G_r_AB;

    % components
    U_AB.f = 1./r; % Electrostatics function f(r)
    U_AB.g = -f6.AB.*B.AB./(r.^6) + G_r_AB; % Dispersion g(r)
    U_AB.h = A.AB./(r.^12);% Short range repulsion

    % Plus - Plus total derivative
    U_AB.dTotal = -k_0*(e_c^2)*q.(Atom_A)*q.(Atom_B)./(r.^2) ...
        - A.AB.*12./(r.^13) ...
        + f6.AB.*B.AB.*6./(r.^7) ...
        - df6.AB.*B.AB./(r.^6) ...
        + dG_r_AB;

    % components
    U_AB.df = 1./(r.^2);% Electrostatics function
    U_AB.dg = - f6.AB.*B.AB.*6./(r.^7) + df6.AB.*B.AB./(r.^6) - dG_r_AB; % Dispersion
    U_AB.dh = + A.AB.*12./(r.^13);% Short range repulsion

    if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
        EVDW_Cutoff = A.AB./(RVDW_Cutoff.^12) ...
        - f6_cutoff.AB.*B.AB./(RVDW_Cutoff.^6) ...
        + G_r_AB_Cutoff;

        % Shift by the dispersion energy at vdw cutoff radius. only affects one
        % energy component, not derivatives (i.e. forces)
        U_AB.Total = U_AB.Total - EVDW_Cutoff;
        U_AB.g = U_AB.g - EVDW_Cutoff;
    end

    % remove infinities
    U_AB = Remove_Infinities(U_AB);
end

%% Print
AB = [r ; U_AB.f ; U_AB.df ; U_AB.g ; U_AB.dg ; U_AB.h ; U_AB.dh];
U_AB_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],AB(:)) );


%% PLOT if plotswitch chosen
if plotswitch
    figure;
    % Options
    lw=3;
    fs=35;

    h = cell(1,8);
    hold on


%     h{1} = plot(r.*10,-k_0*(e_c^2).*U_AB.f + U_AB.h,'Color','r','LineWidth',lw,'LineStyle','-');
%     h{2} = plot(r.*10,k_0*(e_c^2).*U_AA.f + U_AA.h,'Color','b','LineWidth',lw,'Linestyle','-');
%     h{3} = plot(r.*10,k_0*(e_c^2).*U_BB.f + U_BB.h,'Color','g','LineWidth',lw,'Linestyle','-');
%     h{1} = plot(r.*10,U_AB.g,'Color','r','LineWidth',lw,'LineStyle','-');
%     h{2} = plot(r.*10,U_AA.g,'Color','b','LineWidth',lw,'Linestyle','-');
%     h{3} = plot(r.*10,U_BB.g,'Color','g','LineWidth',lw,'Linestyle','-');
%     h{1} = plot(r.*10,-k_0*(e_c^2).*U_AB.f + U_AB.g + U_AB.h,'Color','r','LineWidth',lw,'LineStyle','-');
%     h{2} = plot(r.*10,k_0*(e_c^2).*U_AA.f + U_AA.g + U_AA.h,'Color','b','LineWidth',lw,'Linestyle','-');
%     h{3} = plot(r.*10,k_0*(e_c^2).*U_BB.f + U_BB.g + U_BB.h,'Color','g','LineWidth',lw,'Linestyle','-');
    h{1} = plot(r.*10,U_AB.Total,'Color','r','LineWidth',lw,'LineStyle','-');
    h{2} = plot(r.*10,k_0*(e_c^2).*U_AB.f,'Color','r','LineWidth',lw,'LineStyle','-');

%     GOI = gradient(U_AB.Total,Spacing);

    %title(['Plot of JC Potentials for ' Salt],...
    %    'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('Separation [\AA]','fontsize',fs,'Interpreter','latex');
    ylabel('Potential Energy [kJ mol$^{-1}$]','fontsize',fs,'Interpreter','latex');

    ylim([-600 1000]);
%    ylim([-10 10]);
    xlim([0.5 6.5]);

    % Blank line
    hline = refline([0 0]);
    hline.Color = 'k';
    hline.LineWidth = lw-1;
    hline.LineStyle = '--';
    leg1 = legend([h{:}],{[Atom_A '$^{' num2str(q.(Atom_A),'%+.2f') '}$' ' - ' Atom_B '$^{' num2str(q.(Atom_B),'%+.2f') '}$'] ...
        [Atom_A '$^{' num2str(q.(Atom_A),'%+.2f') '}$' ' - ' Atom_B '$^{' num2str(q.(Atom_B),'%+.2f') '}$ Added']});
    leg1.Interpreter = 'latex';
end

end