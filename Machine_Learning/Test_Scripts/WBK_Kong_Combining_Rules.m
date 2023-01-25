%% Initialize
epsilon_ii = 10^(-6+rand*10);
gamma_ii = 9+rand*11;
sigma_ii = rand*0.5;

epsilon_jj = 10^(-6+rand*10);
gamma_jj = 9+rand*11;
sigma_jj = rand*0.5;

%% Define approximate tight form constants (valid as r->0)
A_ii = 6*epsilon_ii*exp(gamma_ii)/gamma_ii;
B_ii = gamma_ii/sigma_ii;
C_ii = 2*epsilon_ii*(gamma_ii + 3)/gamma_ii;

A_jj = 6*epsilon_jj*exp(gamma_jj)/gamma_jj;
B_jj = gamma_jj/sigma_jj;
C_jj = 2*epsilon_jj*(gamma_jj + 3)/gamma_jj;

%% Apply Kong mixing rules
A_ij = (1/2).*( A_ii.*(A_ii.*B_ii./(A_jj.*B_jj)).^(-B_ii./(B_ii + B_jj)) + ...
                A_jj.*(A_jj.*B_jj./(A_ii.*B_ii)).^(-B_jj./(B_ii + B_jj)) );
B_ij = 2.*B_ii.*B_jj./(B_ii + B_jj);
C_ij = sqrt(C_ii.*C_jj);


%% Convert back. May only be valid under some sets of parameters...
gamma_ij = -lambertw(-1,-3*C_ij/(A_ij*exp(3))) - 3;
sigma_ij = gamma_ij/B_ij;
epsilon_ij = gamma_ij*C_ij/(2*(gamma_ij + 3));
epsilon_ij = A_ij*gamma_ij*exp(-gamma_ij)/6;

disp(['Gamma_ii = ' num2str(gamma_ii) '   Gamma_jj = ' num2str(gamma_jj) '   Gamma_ij = ' num2str(real(gamma_ij))])
disp(['Sigma_ii = ' num2str(sigma_ii) '   Sigma_jj = ' num2str(sigma_jj) '   Sigma_ij = ' num2str(real(sigma_ij))])
disp(['Epsilon_ii = ' num2str(epsilon_ii) '   Epsilon_jj = ' num2str(epsilon_jj) '   Epsilon_ij = ' num2str(real(epsilon_ij))])
