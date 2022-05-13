classdef GromObject
% Damp_Types:
% 0 = no (default) damping. This is default of JC model.
% 1 = BJ/rational damping (same as in D3(BJ))
% 2 = Tang Damping (Essentially removes dispersion in JC)
% 3 = MMDRE Damping function (Very weak damping, damps mainly at mid range)
% 4 = PAMoC Damping function (fairly weak damping, damps mainly at mid range)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)

% TF Parameter sets for C6/C8 coefficients
% 0 = default TF Parameters
% 1 = D3 values
% 2 = Best literature values available
% 3 = D4 with C6 and C8 generated on the fly

% GAdjust are N x 3 arrays of gaussian parameters
% (i , 1) is the Gaussian height of the ith adjustment (may be negative or
% positive)
% (i , 2) is the center point of the ith Gaussian (should be positive)
% (i , 3) is the standard deviation or width (negative and positive values
% are the same)
    properties
        Structure = 'Rocksalt';
        Salt = 'LiF'; % LiF LiCl LiBr LiI and NaCl
        Model = 'JC'; % model(s) to use: JC, JC3P, JC4P, TF
        Damping = 0; % damping functions to use: 0 to 6 are defined
        TF_Paramset = 0; % choose from 0-3
        Scale_Dispersion = 1; % Works for both JC and TF
        Scale_Repulsion = 1; % Works for both JC and TF
        Scale_MM_Dispersion = 1; % Works for both JC and TF
        Scale_XX_Dispersion = 1; % Works for both JC and TF
        Scale_MX_Dispersion = 1; % Works for both JC and TF
        Scale_Epsilon = 1; % Scale all Epsilon (affects JC only)
        Scale_Sigma = 1; % Scale all Sigma (affects JC only)
        Scale_Alpha = 1; % Scale the repulsive exponential parameter alpha (affects TF only)
        Data_Type = 1; % 1 = cell optimization, 2 = full optimization
        GAdjust_MX = [0 0 1]; % Gaussian Adjustment to the potential
        GAdjust_MM = [0 0 1]; % Gaussian Adjustment to the potential
        GAdjust_XX = [0 0 1]; % Gaussian Adjustment to the potential
    end
    methods
        % Initialize some default parameters for all variables
%         function obj = Initialize(obj)
% 
%         end
    end
end