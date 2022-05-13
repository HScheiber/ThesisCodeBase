% Input information
SF = 4.032;
A = SF.*[1.0000000000	0.0000000000	0.0000000000];
B = SF.*[-0.5000000000	0.8660254038	0.0000000000];
C = SF.*[0.0000000000	0.0000000000	1.6329931616];
M_input    = [A;  B;  C];
Li1_Basis1 = [1/3 2/3 0.375];
Li2_Basis1 = [2/3 1/3 0.875];
Cl1_Basis1 = [1/3 2/3 0.000];
Cl2_Basis1 = [2/3 1/3 0.500];
disp('Basis 1')
disp([num2str(mod(Li1_Basis1,1),'%1.12f ') newline num2str(mod(Li2_Basis1,1),'%1.12f ') ...
    newline num2str(mod(Cl1_Basis1,1),'%1.12f ') newline num2str(mod(Cl2_Basis1,1),'%1.12f ')]);

% Convert to cartesian
Li1_Cartesian = Li1_Basis1*M_input;
Li2_Cartesian = Li2_Basis1*M_input;
Cl1_Cartesian = Cl1_Basis1*M_input;
Cl2_Cartesian = Cl2_Basis1*M_input;
disp('Cartesian')
disp([num2str(Li1_Cartesian,'%1.12f ') newline num2str(Li2_Cartesian,'%1.12f ') ...
    newline num2str(Cl1_Cartesian,'%1.12f ') newline num2str(Cl2_Cartesian,'%1.12f ')]);


% Convert back to basis 1
Li1_Basis1b = Li1_Cartesian/M_input;
Li2_Basis1b = Li2_Cartesian/M_input;
Cl1_Basis1b = Cl1_Cartesian/M_input;
Cl2_Basis1b = Cl2_Cartesian/M_input;
disp('Back to basis 1')
disp([num2str(mod(Li1_Basis1b,1),'%1.12f ') newline num2str(mod(Li2_Basis1b,1),'%1.12f ') ...
    newline num2str(mod(Cl1_Basis1b,1),'%1.12f ') newline num2str(mod(Cl2_Basis1b,1),'%1.12f ')]);

disp('Basis 2')



% Output information
A_out = SF.*[1       0.000           0.000];
B_out = SF.*[0.500   0.8660254038    0.000];
C_out = SF.*[0.000	 0.000           1.6329931616];
M_output = [A_out; B_out; C_out];
disp('New structure')
disp('1.000')
disp(num2str(M_output,'%1.12f '))

% Convert into basis 2
Li1_Basis2 = Li1_Cartesian/M_output;
Li2_Basis2 = Li2_Cartesian/M_output;
Cl1_Basis2 = Cl1_Cartesian/M_output;
Cl2_Basis2 = Cl2_Cartesian/M_output;
disp(newline)
disp([num2str(mod(Li1_Basis2,1),'%1.12f ') newline num2str(mod(Li2_Basis2,1),'%1.12f ') ...
    newline num2str(mod(Cl1_Basis2,1),'%1.12f ') newline num2str(mod(Cl2_Basis2,1),'%1.12f ')]);