function Transform_Matrix = GenTransformMatrix(Input)
% Get angles from input structure
alpha = Input.alpha;
beta = Input.beta;
gamma = Input.gamma;

% Generate the unit transformation matrix, the resulting matrix is the
% crystallographic unit cell vectors in cartesian coordinates.
% For simplicity, it is chosen so that edge vector a in the positive x-axis direction, 
% edge vector b in the x-y plane with positive y-axis component, and
% edge vector c with positive z-axis component in the Cartesian-system.
% See https://en.wikipedia.org/wiki/Fractional_coordinates#In_crystallography
Omega = sqrt(1 - (cosd(alpha).^2) - (cosd(beta).^2) - (cosd(gamma).^2) + 2*cosd(alpha)*cosd(beta)*cosd(gamma));
cx= cosd(beta);
cy= (cosd(alpha) - cosd(gamma)*cosd(beta))./sind(gamma);
cz= (1./sind(gamma)).*Omega;

Transform_Matrix = [1           0           0; ...
                    cosd(gamma) sind(gamma) 0; ...
                    cx          cy          cz];

end
