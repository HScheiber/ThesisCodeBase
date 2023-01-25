%% Loose form of exp6
syms epsilonij gammaij sigmaij

epsilonii = rand*10;
gammaii = 7 + rand*10;
sigmaii = rand*10;

epsilonjj = rand*10;
gammajj = 7 + rand*10;
sigmajj = rand*10;


cr1 =  ((epsilonij.*gammaij.*exp(gammaij))./((gammaij - 6).*sigmaij )).^(2.*sigmaij ./ gammaij ) == ((( epsilonii.*gammaii.*exp(gammaii))./(( gammaii - 6).*sigmaii )).^( sigmaii ./ gammaii )).*((( epsilonjj.*gammajj.*exp( gammajj ))./(( gammajj -6).*sigmajj )).^( sigmajj ./ gammajj ));
cr2 = sigmaij/gammaij == 1/2*((sigmaii/gammaii) + (sigmajj/gammajj));
cr3 = (epsilonij*gammaij*sigmaij.^6)./(gammaij - 6) == ((epsilonii.*gammaii.*(sigmaii.^6))./(gammaii - 6).*(epsilonjj.*gammajj.*(sigmajj.^6))./(gammajj - 6)).^(1./2);

[epsilonij,gammaij,sigmaij] = vpasolve([cr1, cr2, cr3], [epsilonij,gammaij,sigmaij]);

epsilonij = double(epsilonij);
gammaij = double(gammaij);
sigmaij = double(sigmaij);

%% Convert back to tight form
Aii = 6*epsilonii*exp(gammaii)/(gammaii-6);
Ajj = 6*epsilonjj*exp(gammajj)/(gammajj-6);
Aij = 6*epsilonij*exp(gammaij)/(gammaij-6);

Bii = gammaii/sigmaii;
Bjj = gammajj/sigmajj;
Bij = gammaij/sigmaij;

Cii = gammaii*epsilonii*(sigmaii^6)/(gammaii-6);
Cjj = gammajj*epsilonjj*(sigmajj^6)/(gammajj-6);
Cij = gammaij*epsilonij*(sigmaij^6)/(gammaij-6);


%% Test combining rules

% Bij
Bij_b = (2.*Bii.*Bjj)./(Bii+Bjj);
test_B = Bij == Bij_b

% Aij in two forms
Aij_1 = (1./2).*(Aii.*((Aii.*Bii)./(Ajj.*Bjj)).^(-(Bii./(Bii + Bjj))) + Ajj.*((Ajj.*Bjj)./(Aii.*Bii)).^(-(Bjj./(Bii+Bjj))));
Aij_2 = ( (Aii.*Bii).^(Bij./(2.*Bii)).*(Ajj.*Bjj).^(Bij./(2.*Bjj)) )./Bij;

% Aij
Aij_1 - Aij
Aij_2 - Aij

% Cij
Cij_2 = sqrt(Cii*Cjj);
Cij_2 - Cij
