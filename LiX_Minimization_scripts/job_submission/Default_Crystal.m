function Geometry = Default_Crystal(Settings,varargin)

p = inputParser;
p.FunctionName = 'Default_Crystal';
addOptional(p,'Scale',Init_Scaling_Object,@(x)validateattributes(x,{'struct'},{'nonempty'}))
addOptional(p,'Center_Coordinates',true,@(x)validateattributes(x,{'logical'},{'nonempty'}))
parse(p,varargin{:});

Scale = p.Results.Scale;
Center_Coordinates = p.Results.Center_Coordinates;


% Default scale factors
c_over_a = sqrt(8/3);
switch Settings.Salt
    case 'LiF'
        SF = 4.00;
    case 'LiCl'
        SF = 5.09;
    case 'LiBr'
        SF = 5.47;
    case 'LiI'
        SF = 5.96;
    case 'LiAt'
        SF = 6.10;
    case 'NaCl'
        SF = 5.56;
    case 'CsCl'
        SF = 7;
    case 'MX'
        
        if Scale.r_d.MM/Scale.r_d.XX >= 0.415
            SF = (9.0755*Scale.r_d.XX + 0.8162+0.001)*10*Scale.r_d.MM/0.77;
        else
            SF = (11.1111*Scale.r_d.XX + 0.001)*10*Scale.r_d.MM/0.77;
        end
        
        %SF = (10*2*Scale.r_d.XX)/0.77;  % 0.2 -> 4.2 while 1 -> 14
%         if Scale.r_d.MM/Scale.r_d.XX > 0.414
%             c_over_a = 1.75334;
%         else
%             c_over_a = sqrt(8/3);
%         end
    otherwise
        SF = 5;
end

%DFT = Load_Best_DFT_Data;

if Settings.Use_Conv_cell
    Geometry.Conv = true;
    switch lower(Settings.Structure)
        case 'betabeo'
            Geometry.a = 1.32*SF; 
            Geometry.b = Geometry.a;
            Geometry.c = (1/sqrt(3)).*Geometry.a;
            x = 0.336; %0.336 Experimental value
            y = 0.310; %0.310 Experimental value
            Geometry.FC_Metal  =  [x     1-x   0.000;...
                                 1-x   x     0.000;...
                                 x+0.5 x+0.5 0.500;...
                                 0.5-x 0.5-x 0.500];
            Geometry.FC_Halide =  [y     y     0.000;...
                                 1-y   1-y   0.000;...
                                 y+0.5 0.5-y 0.500;...
                                 0.5-y y+0.5 0.500];
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 90;
            Geometry.N = 8;
        case 'cscl'
            Geometry.a = 0.65*SF;
            Geometry.b = Geometry.a;
            Geometry.c = Geometry.a;
            Geometry.FC_Metal  = [0.0 0.0 0.0];
            Geometry.FC_Halide = [1/2 1/2 1/2];
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 90;
            Geometry.N = 2;
        case 'fivefive'
            Geometry.a = 1.04*SF;
            Geometry.b = (3/2)*Geometry.a; %; % Theoretical b/a = 3/2
            Geometry.c = (sqrt(3)/2).*Geometry.a; % Theoretical c/a = cosd(30)
            Geometry.FC_Metal  = [1/4 1/6 1/2];
            Geometry.FC_Halide = [1/4 1/3 0.0];
            [Geometry.FC_Metal,Geometry.FC_Halide] = ...
                UnitCell_FractionalCoords(Geometry.FC_Metal,...
                Geometry.FC_Halide,'FiveFive');
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 90;
            Geometry.N = 8;
        case 'nias'
            Geometry.a = 0.71*SF;
            Geometry.b = Geometry.a;
            Geometry.c = c_over_a*Geometry.a;
            Geometry.FC_Metal  = [0.0 0.0 0.0;...
                                0.0 0.0 1/2];
            Geometry.FC_Halide = [1/3 2/3 1/4;...
                                2/3 1/3 3/4];
            Geometry.b = Geometry.a;
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 120;
            Geometry.N = 4;
        case 'antinias'
            Geometry.a = 0.69*SF;
            Geometry.b = Geometry.a;
            Geometry.c = c_over_a*Geometry.a;
            Geometry.FC_Metal  = [1/3 2/3 1/4;...
                                2/3 1/3 3/4];
            Geometry.FC_Halide = [0.0 0.0 0.0;...
                                0.0 0.0 1/2];
            Geometry.b = Geometry.a;
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 120;
            Geometry.N = 4;
        case 'rocksalt'
            Geometry.a = SF;
            Geometry.b = Geometry.a;
            Geometry.c = Geometry.a;
            Geometry.FC_Metal  = [0.0 0.0 0.0];
            Geometry.FC_Halide = [1/2 1/2 1/2];
            
            [Geometry.FC_Metal,Geometry.FC_Halide] = ...
                UnitCell_FractionalCoords(Geometry.FC_Metal,...
                Geometry.FC_Halide,'Rocksalt');
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 90;
            Geometry.N = 8;
        case 'sphalerite'
            Geometry.a = 1.08*SF;
            Geometry.b = Geometry.a;
            Geometry.c = Geometry.a;
            Geometry.FC_Metal  = [0.0 0.0 0.0];
            Geometry.FC_Halide = [1/4 1/4 1/4];
            [Geometry.FC_Metal,Geometry.FC_Halide] = ...
                UnitCell_FractionalCoords(Geometry.FC_Metal,...
                Geometry.FC_Halide,'Sphalerite');
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 90;
            Geometry.N = 8;
        case 'wurtzite'
            Geometry.a = 0.77*SF;
            Geometry.b = Geometry.a;
            Geometry.c = c_over_a*Geometry.a; % Perfect Wurtzite c/a = sqrt(8/3);
            Geometry.FC_Metal  = [1/3 2/3 0.0];
            Geometry.FC_Halide = [1/3 2/3 3/8];
            [Geometry.FC_Metal,Geometry.FC_Halide] = ...
                UnitCell_FractionalCoords(Geometry.FC_Metal,...
                Geometry.FC_Halide,'Wurtzite');
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 120;
            Geometry.N = 4;
        case {'liquid' 'previous'}
            Geometry.a = nan;
            Geometry.b = nan;
            Geometry.c = nan;
            Geometry.FC_Metal  = [];
            Geometry.FC_Halide = [];
            Geometry.alpha = 90;
            Geometry.beta  = 90;
            Geometry.gamma = 90;
            Geometry.Transform = eye(3);
            Geometry.N = nan;
    end
    
else % primitive cell
    Geometry.Conv = false;
    switch lower(Settings.Structure)
        case 'betabeo'
            Geometry.a = 1.32*SF; 
            Geometry.b = Geometry.a;
            Geometry.c = (1/sqrt(3)).*Geometry.a;
            x = 0.336; %0.336 Experimental value
            y = 0.310; %0.310 Experimental value
            Geometry.FC_Metal  =  [x     1-x   0.000;...
                                 1-x   x     0.000;...
                                 x+0.5 x+0.5 0.500;...
                                 0.5-x 0.5-x 0.500];
            Geometry.FC_Halide =  [y     y     0.000;...
                                 1-y   1-y   0.000;...
                                 y+0.5 0.5-y 0.500;...
                                 0.5-y y+0.5 0.500];
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 90;
            Geometry.N = 8;
        case 'cscl'
            Geometry.a = 0.65*SF;
            Geometry.b = Geometry.a;
            Geometry.c = Geometry.a;
            Geometry.FC_Metal  = [0.0 0.0 0.0];
            Geometry.FC_Halide = [1/2 1/2 1/2];
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 90;
            Geometry.N = 2;
        case 'fivefive'
            Geometry.a = 0.83*SF;
            Geometry.b = Geometry.a; %; % Theoretical b/a = 3/2
            Geometry.c = (2/sqrt(3))*Geometry.a; % Theoretical c/a = cosd(30)
            Geometry.FC_Metal  = [1/3 2/3 1/2;...
                                2/3 1/3 0.0];
            Geometry.FC_Halide = [1/3 2/3 0.0;...
                                2/3 1/3 1/2];
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 120;
            Geometry.N = 4;
        case 'nias'
            Geometry.a = 0.71*SF;
            Geometry.b = Geometry.a;
            Geometry.c = c_over_a*Geometry.a;
            Geometry.FC_Metal  = [0.0 0.0 0.0;...
                                0.0 0.0 1/2];
            Geometry.FC_Halide = [1/3 2/3 1/4;...
                                2/3 1/3 3/4];
            Geometry.b = Geometry.a;
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 120;
            Geometry.N = 4;
        case 'antinias'
            Geometry.a = 0.69*SF;
            Geometry.b = Geometry.a;
            Geometry.c = c_over_a*Geometry.a;
            Geometry.FC_Metal  = [1/3 2/3 1/4;...
                                2/3 1/3 3/4];
            Geometry.FC_Halide = [0.0 0.0 0.0;...
                                0.0 0.0 1/2];
            Geometry.b = Geometry.a;
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 120;
            Geometry.N = 4;
        case 'rocksalt'
            Geometry.a = 0.71*SF;
            Geometry.b = Geometry.a;
            Geometry.c = Geometry.a;
            Geometry.FC_Metal  = [0.0 0.0 0.0];
            Geometry.FC_Halide = [1/2 1/2 1/2];
            Geometry.alpha = 60;
            Geometry.beta = 60;
            Geometry.gamma = 60;
            Geometry.N = 2;
        case 'sphalerite'
            Geometry.a = 0.76*SF;
            Geometry.b = Geometry.a;
            Geometry.c = Geometry.a;
            Geometry.FC_Metal  = [0.0 0.0 0.0];
            Geometry.FC_Halide = [1/4 1/4 1/4];
            Geometry.alpha = 60;
            Geometry.beta = 60;
            Geometry.gamma = 60;
            Geometry.N = 2;
        case 'wurtzite'
            Geometry.a = 0.77*SF;
            Geometry.b = Geometry.a;
            Geometry.c = c_over_a*Geometry.a; % Perfect Wurtzite c/a = sqrt(8/3);
            Geometry.FC_Metal  = [1/3 2/3 0.0];
            Geometry.FC_Halide = [1/3 2/3 3/8];
            [Geometry.FC_Metal,Geometry.FC_Halide] = ...
                UnitCell_FractionalCoords(Geometry.FC_Metal,...
                Geometry.FC_Halide,'Wurtzite');
            Geometry.alpha = 90;
            Geometry.beta = 90;
            Geometry.gamma = 120;
            Geometry.N = 4;
        case {'liquid' 'previous'}
            Geometry.a = nan;
            Geometry.b = nan;
            Geometry.c = nan;
            Geometry.FC_Metal  = [];
            Geometry.FC_Halide = [];
            Geometry.alpha = 90;
            Geometry.beta  = 90;
            Geometry.gamma = 90;
            Geometry.Transform = eye(3);
            Geometry.N = nan;
    end
end
Geometry.Transform = GenTransformMatrix(Geometry);
Geometry.NF = Geometry.N/2;
Geometry.Salt = Settings.Salt;
Geometry.Structure = Settings.Structure;
Geometry.Label = StructureLabel(Settings.Structure);
[Geometry.Metal,Geometry.Halide] = Separate_Metal_Halide(Settings.Salt);
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
    Geometry = CenterCoordinates(Geometry,Settings.Structure);
end

end