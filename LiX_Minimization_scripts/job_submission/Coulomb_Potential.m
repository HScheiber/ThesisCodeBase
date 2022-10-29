% Generates Coulomb potential for tables
function [f,df]= Coulomb_Potential(Settings,r,int)
    if Settings.GaussianCharge
        beta_ij = Settings.S.beta.(int);
        f = erf(beta_ij.*r)./r;
        f(r == 0) = beta_ij*2/sqrt(pi);
        
        df = erf(beta_ij.*r)./(r.^2) - (2*beta_ij/sqrt(pi)).*exp(-(beta_ij.*r).^2)./r;
        df(r == 0) = 0;
    else
        f = 1./r; % Electrostatics function f(r)
        df = 1./(r.^2); % Electrostatics function (not including Coulomb constant or charges)
    end
end