function U = BD_Potential_Generator_vec(Settings)

%% Load BH Parameters
dat = load(fullfile(Settings.home,'data','BH_Default_Param.mat'));
Param = dat.Param;
b = dat.b;
rho = dat.rho;
sigma = dat.sigma;
valence = dat.valence;
QQ_prefactor = dat.QQ_prefactor;

%% Parameter: q (charge)
q.(Settings.Metal) =  Settings.S.Q; % atomic
q.(Settings.Halide)= -Settings.S.Q; % atomic

%% Calculate Pauling Coefficients beta: MX = +-   MM = ++     XX = --
beta.MM = 1 + 2/valence.(Settings.Metal); % Unitless
beta.XX = 1 - 2/valence.(Settings.Halide); % Unitless
beta.MX = 1 + 1/valence.(Settings.Metal) - 1/valence.(Settings.Halide); % Unitless

%% Calculate Repulsive Exponential Parameter alpha: MX = +-   MM = ++     XX = --
alpha.MM = Settings.S.A.All.*Settings.S.A.MM./rho.(Settings.Salt); % nm^-1
alpha.MX = Settings.S.A.All.*Settings.S.A.MX./rho.(Settings.Salt); % nm^-1
alpha.XX = Settings.S.A.All.*Settings.S.A.XX./rho.(Settings.Salt); % nm^-1

%% Calculate Repulsive Scaling Parameter B: MX = +-   MM = ++     XX = -- (Including scaling)
B.MM = Settings.S.R.All.*Settings.S.R.MM.*beta.MM.*b.*exp(2.*sigma.(Settings.Metal)./rho.(Settings.Salt));
B.XX = Settings.S.R.All.*Settings.S.R.XX.*beta.XX.*b.*exp(2.*sigma.(Settings.Halide)./rho.(Settings.Salt));
B.MX = Settings.S.R.All.*Settings.S.R.MX.*beta.MX.*b.*exp((sigma.(Settings.Metal) + sigma.(Settings.Halide))./rho.(Settings.Salt));

%% Calculate parameters of interest for LJ potential: change parameteters into C6/r6 format and apply mixing rules
CMM_pre = 4*Param.(Settings.Metal).epsilon*(Param.(Settings.Metal).sigma^6);
CXX_pre = 4*Param.(Settings.Halide).epsilon*(Param.(Settings.Halide).sigma^6);
C.MM = Settings.S.D.All.*Settings.S.D.MM.*CMM_pre;
C.XX = Settings.S.D.All.*Settings.S.D.XX.*CXX_pre;
C.MX = Settings.S.D.All.*Settings.S.D.MX.*sqrt(CMM_pre.*CXX_pre);

%% Generate range (r) in nm
U.r = Settings.Table_StepSize:Settings.Table_StepSize:Settings.Table_Length;

%% If Damping at close range, affects all attractive interactions
for interaction = {'MX' 'XX' 'MM'}
    int = interaction{1};
    switch int
        case 'MX'
            Y1 = Settings.Metal;
            Y2 = Settings.Halide;
        case 'MM'
            Y1 = Settings.Metal;
            Y2 = Settings.Metal;
        case 'XX'
            Y1 = Settings.Halide;
            Y2 = Settings.Halide;
    end
    
    % First build the potential
    U_LJ_all = B.(int).*exp(-alpha.(int).*U.r) - C.(int)./(U.r.^6);
    dU_LJ_all = -alpha.(int).*B.(int).*exp(-alpha.(int).*U.r) + 6.*C.(int)./(U.r.^7);
    
    % Locate the peaks in the LJ part of the potential
    peaks_idx = islocalmax(U_LJ_all,2,'MinProminence',1e-8);
    Num_peaks = sum(peaks_idx,2);
    np_idx = Num_peaks == 0;  % Potentials that contain no peak. do not need to modify these
    op_idx = ~np_idx;         % Potentials that contain at least 1 peak
    
    % For the subset of potentials that have peak(s), find the peak position
    U_LJ = U_LJ_all(op_idx,:);
    peaks_idx  = peaks_idx(op_idx,:);
    r = repmat(U.r,size(U_LJ,1),1);
    fp_idx = sum(cumprod(peaks_idx == 0, 2), 2) + 1;
    fpl_idx = sub2ind(size(U_LJ),(1:numel(fp_idx)).',fp_idx);
    peak_r = nan(size(U_LJ_all,1),1);
    peak_r(op_idx) = r(fpl_idx);
    
    % Find any inflection points that occur after the peak
    r = repmat(U.r,size(dU_LJ_all,1),1);
    inflex_idx = ( islocalmax(dU_LJ_all,2,'MinProminence',1e-8) | islocalmin(dU_LJ_all,2,'MinProminence',1e-8) ) & (r > peak_r); % Pick out inflection points after the peak
    
    % Exclude any with no inflection point or no peak
    Num_inflex = sum(inflex_idx,2);
    np_idx = Num_inflex == 0;  % Potentials that contain no inflection point or peak. exclude these...
    op_idx = ~np_idx;          % Potentials that contain a peak and inflection point
    inflex_idx = inflex_idx(op_idx,:);
    U_LJ = U_LJ_all(op_idx,:);
    dU_LJ = dU_LJ_all(op_idx,:);
    r = repmat(U.r,size(U_LJ,1),1);
    ip_idx = sum(cumprod(inflex_idx == 0, 2), 2) + 1;
    ipl_idx = sub2ind(size(U_LJ),(1:numel(ip_idx)).',ip_idx);
    inflex_r = r(ipl_idx); % inflection positions
    
    % Calculate a coefficient to match the derivative at the inflection point
    dU_infl = dU_LJ(ipl_idx); % value of derivative at inflection point
    D = -dU_infl.*(inflex_r.^13)/12; % coefficients

    % Generate a repulsion beyond the inflection point
    below_infl_idx = (r < inflex_r); % pick out values of r below the inflection point
    fwall = D./(r.^12) - D./(inflex_r.^12) + U_LJ(ipl_idx);
    
    % Add this repulsion to the repulsive part of the function
    U_LJ(below_infl_idx) = fwall(below_infl_idx);
    
    % Update
    U_LJ_all(op_idx,:) = U_LJ;
    
    % Build PES
    U.(int) = QQ_prefactor.*q.(Y1).*q.(Y2)./U.r + U_LJ_all;
end

end