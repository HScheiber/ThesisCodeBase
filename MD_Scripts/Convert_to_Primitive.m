% This function converts a unit cell from conventional to primitive
function Geometry = Convert_to_Primitive(Geometry,varargin)

    if ~Geometry.Conv
        return
    else
        Geometry.Conv = false;
    end
    
    p = inputParser;
    p.FunctionName = 'Convert_to_Primitive';
    addOptional(p,'Center_Coordinates',true,@(x)validateattributes(x,{'logical'},{'nonempty'}))
    parse(p,varargin{:});
    Center_Coordinates = p.Results.Center_Coordinates;
    
    switch lower(Geometry.Structure)
        case 'rocksalt'
            Geometry.a = Geometry.a/sqrt(2);
            Geometry.b = Geometry.b/sqrt(2);
            Geometry.c = Geometry.c/sqrt(2);
            Geometry.FC_Metal  = [0.0 0.0 0.0];
            Geometry.FC_Halide = [1/2 1/2 1/2];
            Geometry.N = 2;
            Geometry.alpha = 60;
            Geometry.beta = 60;
            Geometry.gamma = 60;
        case 'fivefive'
            a = Geometry.a;
            Geometry.a = (a/2)*1/(sqrt((1/4) + (1/9)*(sind(60)^2)));
            Geometry.b = (a/2)*1/(sqrt((1/4) + (1/9)*(sind(60)^2)));
            Geometry.c = a;
            Geometry.FC_Metal  = [1/3 2/3 1/2;...
                                2/3 1/3 0.0];
            Geometry.FC_Halide = [1/3 2/3 0.0;...
                                2/3 1/3 1/2];
            Geometry.N = 4;
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 120;
        case 'sphalerite'
            Geometry.a = Geometry.a/sqrt(2);
            Geometry.b = Geometry.b/sqrt(2);
            Geometry.c = Geometry.c/sqrt(2);
            Geometry.FC_Metal  = [0.0 0.0 0.0];
            Geometry.FC_Halide = [1/4 1/4 1/4];
            Geometry.N = 2;
            Geometry.alpha = 60;
            Geometry.beta = 60;
            Geometry.gamma = 60;
        otherwise
            % Conventional = Primitive (No change)
            return
    end
    
    Geometry.Transform = GenTransformMatrix(Geometry);
    Geometry.NF = Geometry.N/2;
    Geometry.SkewFactor = det(Geometry.Transform);
    Geometry.FC = [Geometry.FC_Metal; Geometry.FC_Halide];
    
    % Calculate skew factors of the transform matrix
    a_vec = Geometry.Transform(1,:);
    b_vec = Geometry.Transform(2,:);
    c_vec = Geometry.Transform(3,:);
    skew_ab = det([ a_vec([1 2])./norm(a_vec([1 2])); b_vec([1 2])./norm(b_vec([1 2])) ]);
    skew_bc = det([ b_vec([2 3])./norm(b_vec([2 3])); c_vec([2 3])./norm(c_vec([2 3])) ]);
    skew_ac = det([ a_vec([1 3])./norm(a_vec([1 3])); c_vec([1 3])./norm(c_vec([1 3])) ]);
    Geometry.Skew_a = min([skew_ab skew_ac]);
    Geometry.Skew_b = min([skew_ab skew_bc]);
    Geometry.Skew_c = min([skew_bc skew_ac]);
    
    if isnan(Geometry.NF)
        Geometry.AtomNames = {};
    else
        Geometry.AtomNames = [repmat({Geometry.Metal},Geometry.NF,1); repmat({Geometry.Halide},Geometry.NF,1)];
    end
    
    if Center_Coordinates
        Geometry = CenterCoordinates(Geometry,Geometry.Structure);
    end

end