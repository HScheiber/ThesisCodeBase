function [CellParams] = Scale_Structure(CellParams,SF)

b_a = CellParams.b./CellParams.a;
c_a = CellParams.c./CellParams.a;

CellParams.a = CellParams.a*SF;
CellParams.b = CellParams.a*b_a;
CellParams.c = CellParams.a*c_a;
CellParams.BL = CellParams.BL*SF;
CellParams.BL2 = CellParams.BL2*SF;
end