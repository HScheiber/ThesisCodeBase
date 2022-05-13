function [Ql,Wl] = calc_Wl(wb,set_l,pbc_on,Ref_traj,TM,sel_coords_t,SelRadius,...
    First_Shell_Switch,Second_Shell_Switch,N_Neigh,Max_R,Min_R)

    Ql = zeros(1,length(set_l));
    Wl = zeros(1,length(set_l));
    
    % Find PBC image closest to reference atom
    if pbc_on
        % Get a set of reference coordinates in fractional coordinates
        ref_coords = mod(Ref_traj/TM,1);

        % Recenter unit cell around the reference atom, applying PBC
        sel_coords = mod(sel_coords_t - ref_coords + (1/2),1);

        % Convert back to cartesian coordinates
        sel_coords = sel_coords*TM;
        ref_coords = [1/2 1/2 1/2]*TM;
    else
        ref_coords = Ref_traj;
        sel_coords = sel_coords_t;
    end

    % Calculate distance between ref and all others
    dR = vecnorm((ref_coords - sel_coords),2,2);

    % Remove reference atom from sel if present
    sel_coords(dR < 1e-5,:) = [];
    dR(dR < 1e-5) = [];

    % Take only those atoms within given radius
    if SelRadius || First_Shell_Switch || Second_Shell_Switch
        Clus_idx = ( (dR >= Min_R) & (dR <= Max_R) );

        R_Clus = dR(Clus_idx);
        N_Clus = length(R_Clus);
        if N_Clus < 1
            % Ignore escaping ion pairs
            Ql = nan(1,length(set_l));
            Wl = nan(1,length(set_l));
            return
        end
        sel_Clus = sel_coords(Clus_idx,:);
        R = dR(Clus_idx,:);

    % Take only those atoms up to given nearest neighbour number
    else
        [R,Sort_idx] = sort(dR);
        sel_coords_sort = sel_coords(Sort_idx,:);
        sel_Clus = sel_coords_sort(1:N_Neigh,:);
        R = R(1:N_Neigh);
    end

    % Cartesian Bond vectors of cluster w.r.t. reference atom
    sel_Bond = sel_Clus - ref_coords;

    %sel_Bond = rodrigues_rot(sel_Bond,[0 0 1],pi/2); %
    %rotate the frame of reference (testing)
    X = sel_Bond(:,1);
    Y = sel_Bond(:,2);
    Z = sel_Bond(:,3);

    % Azimuth = Phi
    % Inclination (polar angle) = theta
    % NOTE: The matlab function cart2sph outputs the angle THETA incorrectly (by pi/2) for use in the function harmonicY 
    % The Theta must be corrected by theta = pi/2-th
    %R = sqrt(X.^2 + Y.^2 + Z.^2);
    Theta = acos(Z./R);
    Phi = atan2(Y,X); % do NOT use atan for this. Must be defined in all four quadrants

    % Calculate Ql and Wl
    for lidx = 1:length(set_l) % Loop through values of l
        l = set_l(lidx);

        % Calculate Ql: Loop through values of m
        qml = 0;
        for m = -l:1:l
            qml_i = harmonicY(l,m,Theta,Phi,'phase',true);
            qml = qml + mean(qml_i)*mean(conj(qml_i));
        end
        Ql(lidx) = (qml*(4*pi/(2*l + 1)))^(1/2);
        
        % Calculate Wl: find all combinations of m that sum to zero
        % Get all triple permutations from the set of m's
        ms = -l:1:l;
        m_triplets = npermutek(ms,3);

        % Select permutations that sum to zero
        subset_zero = (sum(m_triplets,2) == 0);
        m_triplets_zero = m_triplets(subset_zero,:);

        Wl_sub = 0; % initialize
        for m_idx = 1:size(m_triplets_zero,1)
            m_set = m_triplets_zero(m_idx,:);

            qlm1 = harmonicY(l,m_set(1),Theta,Phi,'phase',true);
            qlm2 = harmonicY(l,m_set(2),Theta,Phi,'phase',true);
            qlm3 = harmonicY(l,m_set(3),Theta,Phi,'phase',true);
            w3j = Wigner3j([l l l],m_set);

            Qlm_product = mean(qlm1)*mean(qlm2)*mean(qlm3)*w3j;
            Wl_sub = Wl_sub + Qlm_product;

        end
        %qml_normalize = qml^(3/2);
        Wl(lidx) = real(Wl_sub).*1e2;%/qml_normalize;
        
    end
    increment(wb)
end