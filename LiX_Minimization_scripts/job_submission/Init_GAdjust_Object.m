function [GAdjust_MX,GAdjust_MM,GAdjust_XX] = Init_GAdjust_Object

    % GAdjust are N x 3 arrays of gaussian parameters
    % (i , 1) is the Gaussian height of the ith adjustment (may be negative or
    % positive)
    % (i , 2) is the center point of the ith Gaussian (should be positive)
    % (i , 3) is the standard deviation or width (negative and positive values
    % are the same)

    % Perturb the potential with Gaussians
    GAdjust_MX = [0 0 1];  %[0 0 1];
    GAdjust_MM = [0 0 1];
    GAdjust_XX = [0 0 1];
    
end