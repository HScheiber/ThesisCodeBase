function [Na,Nb,Nc] = SuperCellSize(Settings,Geometry)

Na = ((Geometry.c/Geometry.a)*Settings.N_atoms/(Geometry.N*Settings.c_over_a))^(1/3);
Nc = (Settings.N_atoms/Geometry.N)/(Na^2);
Na = round(Na);
Nb = Na;
Nc = round(Nc);

end