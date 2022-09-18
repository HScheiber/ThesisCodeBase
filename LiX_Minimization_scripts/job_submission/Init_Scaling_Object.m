% Damp_Types:
% 0 = no (default) damping. This is default of JC model.
% 1 = BJ/rational damping (same as in D3(BJ))
% 2 = Tang Damping (Essentially removes dispersion in JC)
% 3 = MMDRE Damping function (Very weak damping, damps mainly at mid range)
% 4 = PAMoC Damping function (fairly weak damping, damps mainly at mid range)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)

% TF Parameter sets for C6/C8 coefficients (for when using as a script)
% 0 = default TF Parameters
% 1 = D3 values
% 2 = Best literature values available
% 3 = D4 with C6 and C8 generated on the fly

function S = Init_Scaling_Object
    % Dispersion
    S.D.All = 1;
    S.D.MM = 1;
    S.D.XX = 1;
    S.D.MX = 1;
    S.D.AB = 1;
    
    % 1/R6 Dispersion (TF only)
    S.D6D.All = 1;
    S.D6D.MM = 1;
    S.D6D.XX = 1;
    S.D6D.MX = 1;
    S.D6D.AB = 1;
    
    % 1/R8 Dispersion (TF only)
    S.D8D.All = 1;
    S.D8D.MM = 1;
    S.D8D.XX = 1;
    S.D8D.MX = 1;
    S.D8D.AB = 1;
    
    % Repulsion
    S.R.All = 1;
    S.R.MM = 1;
    S.R.XX = 1;
    S.R.MX = 1;
    S.R.AB = 1;
    
    % Epsilon (LJ only)
    S.E.All = 1;
    S.E.MM = 1;
    S.E.XX = 1;
    S.E.MX = 1;
    S.E.AB = 1;
    
    % Sigma (LJ only)
    S.S.All = 1;
    S.S.MM = 1;
    S.S.XX = 1;
    S.S.MX = 1;
    S.S.AB = 1;
    
    % Alpha (TF only)
    S.A.All = 1;
    S.A.MM = 1;
    S.A.XX = 1;
    S.A.MX = 1;
    S.A.AB = 1;
    
    % Born parameter n (Mie only)
    S.n.MM = 12;
    S.n.XX = 12;
    S.n.MX = 12;
    
    % Zeta, the final wall potential (HS only)
    S.Z.All = 1;
    S.Z.MM = 1;
    S.Z.XX = 1;
    S.Z.MX = 1;
    S.Z.AB = 1;
    
    % b, the hard sphere wall steepness (HS only)
    S.b.All = 1;
    S.b.MM = 1;
    S.b.XX = 1;
    S.b.MX = 1;
    S.b.AB = 1;
    
    % r_d, the hard sphere wall distance (HS only)
    S.r_d.All = 1;
    S.r_d.MM = 1;
    S.r_d.XX = 1;
    S.r_d.MX = 1;
    S.r_d.AB = 1;
    
    % Scaling charges
    S.Q = 1;
    
    % TF Parameter set
    S.TFParamset = 0;
    
    % Switch to add extra 1/r^N terms to the potential    
    S.N.MM.Value = -1; % Sets the exponent. When value is less than 0, the interaction is turned off
    S.N.MM.Scale = 1;  % Scale of the interaction
    S.N.XX.Value = -1; % Sets the exponent. When value is less than 0, the interaction is turned off
    S.N.XX.Scale = 1;  % Scale of the interaction
    S.N.MX.Value = -1; % Sets the exponent. When value is less than 0, the interaction is turned off
    S.N.MX.Scale = 1;  % Scale of the interaction

end