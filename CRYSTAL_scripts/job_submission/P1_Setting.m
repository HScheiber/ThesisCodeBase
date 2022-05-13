function Cry = P1_Setting(Structure,FC_Metal,FC_Halide)

%% Fixed Structure Settings (Do not Change)
% BetaBeO
if strcmp(Structure,'BetaBeO')
    Cry.alpha = 90;
    Cry.beta = 90;
    Cry.gamma = 90;
    Cry.N = 8;
    Cry.Transform = eye(3);
    [Cry.FC_Metal,Cry.FC_Halide] = ...
        UnitCell_FractionalCoords(FC_Metal,...
        FC_Halide,'BetaBeO');

% CsCl
elseif strcmp(Structure,'CsCl')
    Cry.alpha = 90;
    Cry.beta = 90;
    Cry.gamma = 90;
    Cry.N = 2;
    Cry.Transform = eye(3);
    [Cry.FC_Metal,Cry.FC_Halide] = ...
        UnitCell_FractionalCoords(FC_Metal,...
        FC_Halide,'CsCl');

% FiveFive
elseif strcmp(Structure,'FiveFive')
    Cry.alpha = 90;
    Cry.beta = 90;
    Cry.gamma = 90;
    Cry.N = 8;
    Cry.Transform = eye(3);
    [Cry.FC_Metal,Cry.FC_Halide] = ...
        UnitCell_FractionalCoords(FC_Metal,...
        FC_Halide,'FiveFive');

% NiAs
elseif strcmp(Structure,'NiAs')
    Cry.alpha = 90;
    Cry.beta = 90;
    Cry.gamma = 120;
    Cry.N = 4;
    Cry.Transform =  [1        0        0; ...
                      -sind(30) cosd(30) 0; ...
                       0        0        1];
    [Cry.FC_Metal,Cry.FC_Halide] = ...
        UnitCell_FractionalCoords(FC_Metal,...
        FC_Halide,'NiAs');

% Rocksalt
elseif strcmp(Structure,'Rocksalt')
    Cry.alpha = 90;
    Cry.beta = 90;
    Cry.gamma = 90;
    Cry.N = 8;
    Cry.Transform = eye(3);
    [Cry.FC_Metal,Cry.FC_Halide] = ...
        UnitCell_FractionalCoords(FC_Metal,...
        FC_Halide,'Rocksalt');

% Sphalerite
elseif strcmp(Structure,'Sphalerite')
    Cry.alpha = 90;
    Cry.beta = 90;
    Cry.gamma = 90;
    Cry.N = 8;
    Cry.Transform = eye(3);
    [Cry.FC_Metal,Cry.FC_Halide] = ...
        UnitCell_FractionalCoords(FC_Metal,...
        FC_Halide,'Sphalerite');

% Wurtzite
elseif strcmp(Structure,'Wurtzite')
    Cry.alpha = 90;
    Cry.beta = 90;
    Cry.gamma = 120;
    Cry.N = 4;
    Cry.Transform =  [1        0        0; ...
                      -sind(30) cosd(30) 0; ...
                       0        0        1];
    [Cry.FC_Metal,Cry.FC_Halide] = ...
        UnitCell_FractionalCoords(FC_Metal,...
        FC_Halide,'Wurtzite');
end
end