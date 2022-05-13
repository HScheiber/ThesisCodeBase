% Initializes primitive cells
function Structure_out = initialize_crystal(Salt,Structure)

% Default scale factors
switch Salt
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
end

switch lower(Structure)
    case 'betabeo'
        Structure_out.a = 1.32*SF; 
        Structure_out.b = Structure_out.a;
        Structure_out.c = (1/sqrt(3)).*Structure_out.a;
        x = 0.336; %0.336 Experimental value
        y = 0.310; %0.310 Experimental value
        Structure_out.FC_Metal  = [x     1-x   0.000;...
                                 1-x   x     0.000;...
                                 x+0.5 x+0.5 0.500;...
                                 0.5-x 0.5-x 0.500];
        Structure_out.FC_Halide = [y     y     0.000;...
                                 1-y   1-y   0.000;...
                                 y+0.5 0.5-y 0.500;...
                                 0.5-y y+0.5 0.500];
        Structure_out.alpha = 90;
        Structure_out.beta = 90;
        Structure_out.gamma = 90;
        Structure_out.N = 8;
        Structure_out.Label = 'B';
    case 'cscl'
        Structure_out.a = 0.65*SF;
        Structure_out.b = Structure_out.a;
        Structure_out.c = Structure_out.a;
        Structure_out.FC_Metal  = [0.0 0.0 0.0];
        Structure_out.FC_Halide = [1/2 1/2 1/2];
        Structure_out.alpha = 90;
        Structure_out.beta = 90;
        Structure_out.gamma = 90;
        Structure_out.N = 2;
        Structure_out.Label = 'C';
    case 'fivefive'
        Structure_out.a = 0.83*SF;
        Structure_out.b = Structure_out.a; %; % Theoretical b/a = 3/2
        Structure_out.c = (2/sqrt(3))*Structure_out.a; % Theoretical c/a = cosd(30)
        Structure_out.FC_Metal  = [2/3 2/3 1/2;...
                                  1/3 1/3 0.0];
        Structure_out.FC_Halide = [2/3 2/3 0.0;...
                                  1/3 1/3 1/2];
        Structure_out.alpha = 90;
        Structure_out.beta = 90;
        Structure_out.gamma = 60;
        Structure_out.N = 4;
        Structure_out.Label = 'F';
    case 'nias'
        Structure_out.a = 0.71*SF;
        Structure_out.b = Structure_out.a;
        Structure_out.c = sqrt(8/3)*Structure_out.a;
        Structure_out.FC_Metal  = [0.0 0.0 0.0;...
                              0.0 0.0 1/2];
        Structure_out.FC_Halide = [1/3 1/3 1/4;...
                              2/3 2/3 3/4];
        Structure_out.b = Structure_out.a;
        Structure_out.alpha = 90;
        Structure_out.beta = 90;
        Structure_out.gamma = 60;
        Structure_out.N = 4;
        Structure_out.Label = 'N';
        
    case 'antinias'
        Structure_out.a = 0.69*SF;
        Structure_out.b = Structure_out.a;
        Structure_out.c = sqrt(8/3)*Structure_out.a;
        Structure_out.FC_Metal = [1/3 1/3 1/4;...
                                  2/3 2/3 3/4];
        Structure_out.FC_Halide  = [0.0 0.0 0.0;...
                                    0.0 0.0 1/2];
        Structure_out.b = Structure_out.a;
        Structure_out.alpha = 90;
        Structure_out.beta = 90;
        Structure_out.gamma = 60;
        Structure_out.N = 4;
        Structure_out.Label = 'N';
    case 'rocksalt'
        Structure_out.a = 0.71*SF;
        Structure_out.b = Structure_out.a;
        Structure_out.c = Structure_out.a;
        Structure_out.FC_Metal  = [0.0 0.0 0.0];
        Structure_out.FC_Halide = [1/2 1/2 1/2];
        Structure_out.alpha = 60;
        Structure_out.beta = 60;
        Structure_out.gamma = 60;
        Structure_out.N = 2;
        Structure_out.Label = 'R';
    case 'sphalerite'
        Structure_out.a = 0.76*SF;
        Structure_out.b = Structure_out.a;
        Structure_out.c = Structure_out.a;
        Structure_out.FC_Metal  = [0.0 0.0 0.0];
        Structure_out.FC_Halide = [1/4 1/4 1/4];
        Structure_out.alpha = 60;
        Structure_out.beta = 60;
        Structure_out.gamma = 60;
        Structure_out.N = 2;
        Structure_out.Label = 'S';
    case 'wurtzite'
        Structure_out.a = 0.77*SF;
        Structure_out.b = Structure_out.a;
        Structure_out.c = sqrt(8/3)*Structure_out.a; % Perfect Wurtzite c/a = sqrt(8/3);
        Structure_out.FC_Metal  = [2/3 2/3 3/8;...
                                  1/3 1/3 7/8];
        Structure_out.FC_Halide = [2/3 2/3 0.0;...
                                  1/3 1/3 1/2];
        Structure_out.alpha = 90;
        Structure_out.beta = 90;
        Structure_out.gamma = 60;
        Structure_out.N = 4;
        Structure_out.Label = 'W';
end

Structure_out.SF = Structure_out.a;
end