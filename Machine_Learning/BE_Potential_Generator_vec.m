function U = BE_Potential_Generator_vec(Settings)

%% Load BH Parameters
dat = load(fullfile(Settings.home,'data','BH_Default_Param.mat'));
QQ_prefactor = dat.QQ_prefactor;

%% Parameter: q (charge)
q.M =  Settings.S.Q; % atomic
q.X = -Settings.S.Q; % atomic

%% Parameters of interest for vdW potential
sigma.MX = Settings.S.S.MX; % nm
sigma.MM = Settings.S.S.MM;
sigma.XX = Settings.S.S.XX;

epsilon.MX = Settings.S.E.MX; % kJ/mol
epsilon.MM = Settings.S.E.MM;
epsilon.XX = Settings.S.E.XX;

gamma.MX = Settings.S.G.MX; % unitless
gamma.MM = Settings.S.G.MM;
gamma.XX = Settings.S.G.XX;

%% Convert to B*exp - C6/r^6 form
B.MX = Settings.S.R.All.*Settings.S.R.MX.*(6*epsilon.MX./(gamma.MX - 6)).*exp(gamma.MX);
B.MM = Settings.S.R.All.*Settings.S.R.MM.*(6*epsilon.MM./(gamma.MM - 6)).*exp(gamma.MM);
B.XX = Settings.S.R.All.*Settings.S.R.XX.*(6*epsilon.XX./(gamma.XX - 6)).*exp(gamma.XX);

C.MX = Settings.S.D.All.*Settings.S.D.MX.*(epsilon.MX.*gamma.MX.*(sigma.MX.^6))./(gamma.MX - 6);
C.MM = Settings.S.D.All.*Settings.S.D.MM.*(epsilon.MM.*gamma.MM.*(sigma.MM.^6))./(gamma.MM - 6);
C.XX = Settings.S.D.All.*Settings.S.D.XX.*(epsilon.XX.*gamma.XX.*(sigma.XX.^6))./(gamma.XX - 6);

alpha.MX = Settings.S.A.All.*Settings.S.A.MX.*gamma.MX./sigma.MX; % nm^(-1)
alpha.MM = Settings.S.A.All.*Settings.S.A.MM.*gamma.MM./sigma.MM; % nm^(-1)
alpha.XX = Settings.S.A.All.*Settings.S.A.XX.*gamma.XX./sigma.XX; % nm^(-1)

%% Generate range (r) in nm
U.r = Settings.Table_StepSize:Settings.Table_StepSize:Settings.Table_Length;

%% If Damping at close range, affects all attractive interactions
for interaction = {'MX' 'XX' 'MM'}
    int = interaction{1};
    
    % First build the potential
    r_all = repmat(U.r,numel(alpha.(int)),1); % matrix of distances
    U_Total_all  = QQ_prefactor.*q.(int(1)).*q.(int(2))./r_all + B.(int).*exp(-alpha.(int).*U.r) - C.(int)./(U.r.^6);
    dU_Total_all = -QQ_prefactor.*q.(int(1)).*q.(int(2))./(r_all.^2) - alpha.(int).*B.(int).*exp(-alpha.(int).*U.r) + 6.*C.(int)./(U.r.^7);
    
    % Locate the peaks in the total potential
    peaks_idx = islocalmax(U_Total_all,2,'MinProminence',1e-8);
    Num_peaks = sum(peaks_idx,2);
    np_idx = Num_peaks == 0;  % Potentials that contain no peak. do not need to modify these
    op_idx = ~np_idx;         % Potentials that contain at least 1 peak
    
    % For the subset of potentials that have peak(s), find the peak position
    % Take out the subset of potentials that have peaks
    U_Total = U_Total_all(op_idx,:);
    peaks_idx  = peaks_idx(op_idx,:);
    r_fp = repmat(U.r,size(U_Total,1),1); % matrix of positions
    
    fp_idx  = sum(cumprod(peaks_idx == 0, 2), 2) + 1; % index of first peak
    fpl_idx = sub2ind(size(U_Total),(1:numel(fp_idx)).',fp_idx); % linear index of first peak
    peak_r  = nan(size(U_Total_all,1),1);
    peak_r(op_idx) = r_fp(fpl_idx); % Peak positions for those with a peak, NaN for those without a peak
    
    % Find any inflection points that occur after the peak
    inflex_idx = ( islocalmax(dU_Total_all,2,'MinProminence',1e-8) | islocalmin(dU_Total_all,2,'MinProminence',1e-8) ) & (r_all > peak_r); % Pick out inflection points after the peak
    
    % Exclude any with no inflection point or no peak
    Num_inflex = sum(inflex_idx,2);
    ni_idx = Num_inflex == 0;  % Potentials that contain no inflection point or peak. exclude these...
    oi_idx = ~ni_idx;          % Potentials that contain a peak and inflection point
    if sum(oi_idx) > 0
        inflex_idx = inflex_idx(oi_idx,:); % Inflection point index for potentials with a peak and inflection point
        U_Total = U_Total_all(oi_idx,:);   % Total Potential for those with a peak and inflection point
        dU_Total = dU_Total_all(oi_idx,:); % Total derivative for those with a peak and inflection point
        r_oi = repmat(U.r,size(U_Total,1),1); % matrix of positions for potentials with a peak and inflection point
        ip_idx = sum(cumprod(inflex_idx == 0, 2), 2) + 1; % indeces for inflection points
        ipl_idx = sub2ind(size(U_Total),(1:numel(ip_idx)).',ip_idx); % linear indeces for inflection points
        inflex_r = r_oi(ipl_idx); % inflection positions

        % Calculate a coefficient to match the total derivative at the inflection point
        dU_infl = dU_Total(ipl_idx); % value of total derivative at inflection point
        D = -dU_infl.*(1./alpha.(int)(oi_idx)).*exp(alpha.(int)(oi_idx).*inflex_r); % Calculate the wall prefactor

        % Generate a repulsion beyond the inflection point
        below_infl_idx = (r_oi < inflex_r); % pick out values of r below the inflection point
    
        fwall = D.*exp(-alpha.(int)(oi_idx).*r_oi) - D.*exp(-alpha.(int)(oi_idx).*inflex_r) ...
            + U_Total(ipl_idx); % Repulsive wall shifted to its value at the inflection point, with U_Coulomb subtracted out
        
        % Add this repulsion to the repulsive part of the function
        U_Total(below_infl_idx) = fwall(below_infl_idx);

        % Update
        U_Total_all(oi_idx,:) = U_Total;
    end
    
    % Build PES
    U.(int) = U_Total_all;
end

end