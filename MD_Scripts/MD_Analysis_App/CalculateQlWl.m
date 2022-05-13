% input a cluster of bond vectors in cartesian coordinates
% l_series is a series of l values for Ql and Wl to be calculated wrt to
function [Ql,Wl_norm] = CalculateQlWl(Cluster,l_series)

Ql = zeros(size(l_series));
Wl_norm = zeros(size(l_series));
if isempty(Cluster)
    return
end

% % Get angles w.r.t. cartesian coordinate system
% fileID = fopen('C:\Users\Hayden\Documents\Patey_Lab\FCC_Cluster13.xyz');
% Cluster = textscan(fileID,'%*s %f %f %f');
% Cluster = [Cluster{1} Cluster{2} Cluster{3}];
% fclose(fileID);

Cluster = rodrigues_rot(Cluster,[0 1 5],pi/6);
X = Cluster(:,1);
Y = Cluster(:,2);
Z = Cluster(:,3);

% Azimuth = Phi
% Inclination (polar angle) = theta
% NOTE: The matlab function cart2sph outputs the angle THETA incorrectly (by pi/2) for use in the function harmonicY 
% The Theta must be corrected by theta = pi/2-th
R = sqrt(X.^2 + Y.^2 + Z.^2);
Theta = acos(Z./R);
Phi = atan2(Y,X); % do NOT use atan for this. Must be defined in all four quadrants

% Calculate Ql and Wl
for lidx = 1:length(l_series) % Loop through values of l
    l = l_series(lidx);

    % Calculate Ql: Loop through values of m
    qml = 0;
    for m = -l:1:l
        qml_i = harmonicY(l,m,Theta,Phi,'phase',true);
        qml = qml + mean(qml_i)*mean(conj(qml_i));
    end
    Ql(lidx) = (qml*(4*pi/(2*l + 1)))^(1/2);
    
    %% Common Legendre Polynomial Version
%     % All combinations of bonds in cluster
%     combs_idx = nchoosek(1:length(Cluster),2);
%     
%     VecA  = Cluster(combs_idx(:,1),:);
%     VecB  = Cluster(combs_idx(:,2),:);
%     
%     cos_omegaij = cosd(VecAngle(VecA,VecB,2));
%     Pl = legendreP(l,cos_omegaij);
%     
%     Ql2 = 1/Nb + 2*sum(Pl)/(Nb^2);
%     Ql(lidx) = sqrt(abs(Ql2));
    
    %% Associated Legendre Polynomial Version
%     perms_idx = npermutek(1:length(Cluster),2);
% 
%     Phi_Sum = Phi(perms_idx(:,1)) + Phi(perms_idx(:,2));
%     cos_Theta_A = cos(Theta(perms_idx(:,1)));
%     cos_Theta_B = cos(Theta(perms_idx(:,2)));
%     
%     Polynoms = zeros(2*l+1,length(cos_Theta_A));
%     Polynoms(l+1:end,:) = legendre(l,cos_Theta_A).*legendre(l,cos_Theta_B);
%     Polynoms(1:l,:) = flipud(Polynoms(l+2:end,:));
%     
%     ms = -l:1:l;
%     for m_idx = 1:(2*l+1)
%         m = ms(m_idx);
%         cos_phi_m = cos(m.*Phi_Sum);
%         qml(m_idx) = sum(Polynoms(m_idx,:).*cos_phi_m');
%     end
%     
%     Ql2 = ( 4*pi/(2*l+1) )*sum(qml);
%     
%     Ql(lidx) = sqrt(abs(Ql2));

    %% Calculate Wl: find all combinations of m that sum to zero
    % Get all triple permutations from the set of m's
    ms = -l:1:l;
    m_triplets = npermutek(ms,3);

    % Select permutations that sum to zero
    subset_zero = (sum(m_triplets,2) == 0);
    m_triplets_zero = m_triplets(subset_zero,:);

    Wl = 0; % initialize
    for m_idx = 1:size(m_triplets_zero,1)
        m_set = m_triplets_zero(m_idx,:);
        
        m1 = harmonicY(l,m_set(1),Theta,Phi,'phase',true);
        m2 = harmonicY(l,m_set(2),Theta,Phi,'phase',true);
        m3 = harmonicY(l,m_set(3),Theta,Phi,'phase',true);
        w3j = Wigner3j([l l l],m_set);
        
        Qlm_product = mean(m1)*mean(m2)*mean(m3)*w3j;

        Wl = Wl + Qlm_product;
    end
    %qml_normalize = qml^(3/2);
    Wl_norm(lidx) = real(Wl).*1e2;%/qml_normalize;
end

end