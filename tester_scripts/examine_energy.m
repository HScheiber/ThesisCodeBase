% @ s0 legend "LJ (SR)"
% @ s1 legend "Coulomb (SR)"
% @ s2 legend "Potential"
% @ s3 legend "Temperature"
% @ s4 legend "Volume"
data = import_xvg('energy.xvg');
NF = 15680/2;
plot(data(:,1),data(:,6)./NF)
ylim([-1200 100])




data = load_gro_file('10K-CA2_R_BH_Model_GG2_NPT_1321.3000.gro');

N_sol = 7200;
data.res_number = data.res_number(N_sol+1:end) - N_sol;
data.res_name = data.res_name(N_sol+1:end);
data.atom_name = data.atom_name(N_sol+1:end);
data.atom_number = data.atom_number(N_sol+1:end) - N_sol;
data.xyz = data.xyz(N_sol+1:end,:);
data.vel = data.vel(N_sol+1:end,:);
data.N_atoms = data.N_atoms - N_sol;
data.Salt = 'LiF';
data.N = data.N_atoms;
SaveGroFile('test.gro',data,true);

min(pdist(data.xyz))



View_Gmx_Table('10K-CA2_R_BH_Model_GG2_NPT_Table.xvg',-1)