Sim_atoms = {'Li' 'I' 'Al' 'O'};
Salt = 'LiI';
fix_slab = false;
Scale = Init_Scaling_Object;
% Dispersion
Scale.D.All = 1;
Scale.D.MM = 1; % 910
Scale.D.XX = 1; % 0.23
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
CRDamping.MM.r_d = -1; % liCl = 0.21, LiI = 0.26
CRDamping.MM.b  = -1;
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
Atom_Pairs = [Atom_Pairs repmat({[]},length(Atom_Pairs),2)]; % Column for parameters

for idx = 1:size(Atom_Pairs,1)
    Atom_A = Atom_Pairs{idx,1};
    Atom_B = Atom_Pairs{idx,2};
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

    % 3rd column sigma, 4th column epsilon
    
    [Atom_Pairs{idx,3},Atom_Pairs{idx,4}] = grab_clayff(Atom_A,Atom_B,Salt);
    if fix_slab
        switch Atom_A
            case {'Al' 'O'}
                switch Atom_B
                    case {'Al' 'O'}
                        Atom_Pairs{idx,3} = 0.0;
                        Atom_Pairs{idx,4} = 0.0;
                end
        end
    end
    
    disp(['  ' pad(Atom_A,2) '  ' pad(Atom_B,2) '   1     ' num2str(Atom_Pairs{idx,3},'%.15e') '    ' num2str(Atom_Pairs{idx,4},'%.15e')])
end


