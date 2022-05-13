Sim_atoms = {'Li' 'I' 'Al' 'O'};
Filename_base = 'LiI_L_Sapphire_Slab_001';

Table_Length = 4.01; % How far should tabulated potentials extend in nm
Table_StepSize = 0.0005; % nm
vdW_modifier = 'Potential-shift'; % Potential-shift-Verlet, Potential-shift, None, Force-switch, Potential-switch
RvdW_Cutoff = 1.2; % 1.2 nm. note that rlist ? rCoulomb = RVDW when using Verlet and VerletBT = -1

Scale = Init_Scaling_Object;
% Dispersion
Scale.D.All = 1;
Scale.D.MM = 1;%910;
Scale.D.XX = 1;%0.23;
Scale.D.MX = 1;
% Repulsion
Scale.R.All = 1;
Scale.R.MM = 1;
Scale.R.XX = 1;
Scale.R.MX = 1;
% Epsilon (JC only)
Scale.E.All = 1;
Scale.E.MM = 1;
Scale.E.XX = 1;
Scale.E.MX = 1;
% Sigma (JC only)
Scale.S.All = 1;
Scale.S.MM = 1;
Scale.S.XX = 1;
Scale.S.MX = 1;
% Alpha (TF only)
Scale.A.All = 1;
Scale.A.MM = 1;
Scale.A.XX = 1;
Scale.A.MX = 1;
% Close-Range Sigmoid Damping of attractive potentials
CRDamping.MX.r_d = -1; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
CRDamping.MX.b = -1; % sigmoid "steepness" for damping
CRDamping.MM.r_d = -1;%0.26; % liCl = 0.21, LiI = 0.26
CRDamping.MM.b  = -1;%75; % 75
CRDamping.XX.r_d = -1; 
CRDamping.XX.b  = -1;
CRDamping.AB.r_d = -1; 
CRDamping.AB.b  = -1;
% C6/C8 Dispersion damping
C6Damping.MX = 0; % Set this to 0 for no close range dispersion damping!
C6Damping.MM = 0;
C6Damping.XX = 0;
C6Damping.AB = 0;

% Scale Static Charges (default = 1 for all LiX)
Scale.Q = 1;


Atom_Pairs = nmultichoosek(Sim_atoms,2);

for idx = 1:size(Atom_Pairs,1)
    Atom_A = Atom_Pairs{idx,1};
    Atom_B = Atom_Pairs{idx,2};
    if idx == 1
        TableFile = [Filename_base '_Table.xvg'];
    elseif idx == 2
        energygrp_table = [Atom_A ' ' Atom_B];
        TableFile = [Filename_base '_Table_' Atom_A '_' Atom_B '.xvg'];
    else
        energygrp_table = [energygrp_table '  ' Atom_A ' ' Atom_B];
        TableFile = [Filename_base '_Table_' Atom_A '_' Atom_B '.xvg'];
    end
    
    S = Scale;
    CRDamp = CRDamping;
    C6Damp = C6Damping;
    if strcmp(Atom_A,Atom_B) && strcmpi(Atom_A,'Li')
        S.D.AB = S.D.MM;
        CRDamp.AB = CRDamping.MM;
    elseif strcmp(Atom_A,Atom_B) && strcmpi(Atom_A,'I')
        S.D.AB = S.D.XX;
        CRDamp.AB = CRDamping.XX;
    end
    
    disp(['Potential for: ' Atom_A '-' Atom_B])
    U_AB_out = ClayFF_Potential_Generator(0,Table_Length,Table_StepSize,...
        Atom_A,Atom_B,true,S,vdW_modifier,RvdW_Cutoff,C6Damp,CRDamp);
    
    % Save tables into current directory
    fidMX = fopen(TableFile,'wt');
    fwrite(fidMX,regexprep(U_AB_out,'\r',''));
    fclose(fidMX);
    
end