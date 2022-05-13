% Function that outputs the angle in degrees between two vectors
function Angle = VecAngle(a_vec,b_vec,varargin)
    if nargin == 3
        dim = varargin{1};
        Angle = atan2d(vecnorm(cross(a_vec,b_vec,dim),2,dim),dot(a_vec,b_vec,dim));
    else
        Angle = atan2d(norm(cross(a_vec,b_vec)),dot(a_vec,b_vec));
    end
end
