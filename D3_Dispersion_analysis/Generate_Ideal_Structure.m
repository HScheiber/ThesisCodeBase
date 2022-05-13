function [CellParams] = Generate_Ideal_Structure(CellParams,Structure,BL)

% [Metal,Halide] = Separate_Metal_Halide(Salt);

switch lower(Structure)
    case "rocksalt"
%         CellParams.alpha = 90;
%         CellParams.beta = 90;
%         CellParams.gamma = 90;
        CellParams.BL = BL;
        CellParams.a = 2*BL;
        CellParams.b = CellParams.a;
        CellParams.c = CellParams.a;
%         CellParams.N = 8;
%         CellParams.FC_Metal  = [0.0 0.0 0.0];
%         CellParams.FC_Halide = [1/2 1/2 1/2];   
    case "wurtzite"
%         CellParams.alpha = 90;
%         CellParams.beta = 90;
%         CellParams.gamma = 120;
        CellParams.BL = BL;
        CellParams.a = sqrt(8/3)*BL;
        CellParams.b = CellParams.a;
        CellParams.c = sqrt(8/3)*CellParams.a;
%         CellParams.N = 4;
%         CellParams.FC_Metal  = [1/3 2/3 0.0];
%         CellParams.FC_Halide = [1/3 2/3 3/8];
    case "betabeo"
%         CellParams.alpha = 90;
%         CellParams.beta = 90;
%         CellParams.gamma = 90;
        x = CellParams.x; %0.336 = Experimental value
        y = CellParams.y; %0.310 = Experimental value
        CellParams.BL = BL;
        CellParams.a = (1/sqrt(2*(x^2) + 2*(y^2) - 2*x -2*y + 1))*BL;
        CellParams.b = CellParams.a;
        CellParams.c = (sqrt(3)/2)*CellParams.a;
%         CellParams.N = 8;
%         CellParams.x = 0.336; %0.336; % Experimental value
%         CellParams.y = 0.310; %0.310; % Experimental value
%         CellParams.FC_Metal  = [CellParams.x 1-CellParams.x 0.000];
%         CellParams.FC_Halide = [CellParams.y CellParams.y   0.000];
    case "cscl"
%         CellParams.alpha = 90;
%         CellParams.beta = 90;
%         CellParams.gamma = 90;
        CellParams.BL = BL;
        CellParams.a = (2/sqrt(3))*BL;
        CellParams.b = CellParams.a;
        CellParams.c = CellParams.a;
%         CellParams.FC_Metal  = [0.0 0.0 0.0];
%         CellParams.FC_Halide = [1/2 1/2 1/2];
%         CellParams.N = 2;
    case "fivefive"
%         CellParams.alpha = 90;
%         CellParams.beta = 90;
%         CellParams.gamma = 90;
        CellParams.BL = BL;
        CellParams.a = 2*BL;
        CellParams.b = CellParams.a*(3/2);
        CellParams.c = CellParams.a*sqrt(3)/2;
%         CellParams.N = 8;
%         CellParams.FC_Metal  = [1/4 1/6 1/2];
%         CellParams.FC_Halide = [1/4 1/3 0.0];
    case "nias"
%         CellParams.alpha = 90;
%         CellParams.beta = 90;
%         CellParams.gamma = 120;
        x = CellParams.x;
        CellParams.BL = BL;
        CellParams.a = (1/sqrt((1/3) + ((x^2)/16)))*BL;
        CellParams.b = CellParams.a;
        CellParams.c = CellParams.a*x;
%         CellParams.FC_Metal  = [0.0 0.0 0.0];
%         CellParams.FC_Halide = [1/3 2/3 1/4];
%         CellParams.N = 4;
    case "sphalerite"
%         CellParams.alpha = 90;
%         CellParams.beta = 90;
%         CellParams.gamma = 90;
        CellParams.BL = BL;
        CellParams.a = (4/sqrt(3))*BL;
        CellParams.b = CellParams.a;
        CellParams.c = CellParams.a;
%         CellParams.N = 8;
%         CellParams.FC_Metal  = [0.0 0.0 0.0];
%         CellParams.FC_Halide = [1/4 1/4 1/4];     
    otherwise
        error(['Unknown structure: ' Structure])
end

% [CellParams.FC_Metal,CellParams.FC_Halide] = ...
%     UnitCell_FractionalCoords(CellParams.FC_Metal,...
%     CellParams.FC_Halide,Structure);

% Frac_Coordinates = [CellParams.FC_Metal; CellParams.FC_Halide];
% 
% Identities = [repmat(string(Metal),CellParams.N/2,1);
%     repmat(string(Halide),CellParams.N/2,1)];
end