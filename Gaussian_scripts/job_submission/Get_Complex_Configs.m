function [Mol,Multiplicities] = Get_Complex_Configs(Mol)

% Properties:
% Metal_Charge
% Multiplicities
% Basis_Route_Rough
% Basis_Set_Rough
% Basis_Rough_Aux
% Basis_Route_Tight
% Basis_Set_Tight
% Basis_Route_TD
% Basis_Set_TD

if ispc
    Basis_Set_Loc = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\basis_sets';
else
    Basis_Set_Loc = '/home/scheiber/ThesisCodeBase/basis_sets';
end

switch Mol.Metal_Ion
    case 'Ga'
        Mol.Metal_Charge = 3;
        Mol.Basis_Route_Rough = 'Def2SVP';
        Mol.Basis_Set_Rough = '';
        Mol.Basis_Rough_Aux = 'W06';

        Mol.Basis_Route_Tight = 'Def2TZVP';
        Mol.Basis_Set_Tight = '';
        
        Mol.Basis_Route_TD = 'Gen';
        Mol.Basis_Set_TD = ['@' fullfile(Basis_Set_Loc,'def2-TZVPD.gbs') newline];
        
        Multiplicities = 1;
    case 'Fe'
        Mol.Metal_Charge = 3;
        Mol.Basis_Route_Rough = 'Def2SVP';
        Mol.Basis_Set_Rough = '';
        Mol.Basis_Rough_Aux = 'W06';
        
        Mol.Basis_Route_Tight = 'Def2TZVP';
        Mol.Basis_Set_Tight = '';
        
        Mol.Basis_Route_TD = 'Gen';
        Mol.Basis_Set_TD = ['@' fullfile(Basis_Set_Loc,'def2-TZVPD.gbs') newline];

        Multiplicities = [2 6];
    case 'In'
        Mol.Metal_Charge = 3;
        Mol.Basis_Route_Rough = 'Def2SVP';
        Mol.Basis_Set_Rough = '';
        Mol.Basis_Rough_Aux = 'W06';
        
        Mol.Basis_Route_Tight = 'Def2TZVP'; 
        Mol.Basis_Set_Tight = '';
        
        Mol.Basis_Route_TD = 'GenECP';
        Basis_Set_file = ['@' fullfile(Basis_Set_Loc,'def2-TZVPD.gbs')];
        ECP_file = ['@' fullfile(Basis_Set_Loc,'In_def2-ECP.ecp')];
        Mol.Basis_Set_TD = ...
        [Basis_Set_file newline ...
         newline ...
         ECP_file newline];
     
        Multiplicities = 1;
    case 'Lu'
        Mol.Metal_Charge = 3;
        Mol.Basis_Route_Rough = 'GenECP';
        Basis_Set_file = ['@' fullfile(Basis_Set_Loc,'def2-SVP.gbs')];
        ECP_file = ['@' fullfile(Basis_Set_Loc,'Lu_def2-ECP.ecp')];
        Mol.Basis_Set_Rough = ...
        [Basis_Set_file newline ...
         newline ...
         ECP_file newline newline];
        Mol.Basis_Rough_Aux = 'auto';

        Mol.Basis_Route_Tight = 'GenECP';        
        Basis_Set_file = ['@' fullfile(Basis_Set_Loc,'def2-TZVP.gbs')];
        ECP_file = ['@' fullfile(Basis_Set_Loc,'Lu_def2-ECP.ecp')];
        Mol.Basis_Set_Tight = ...
        [Basis_Set_file newline ...
         newline ...
         ECP_file newline newline];
     
        Mol.Basis_Route_TD = 'GenECP';
        Basis_Set_file = ['@' fullfile(Basis_Set_Loc,'def2-TZVPD.gbs')];
        ECP_file = ['@' fullfile(Basis_Set_Loc,'Lu_def2-ECP.ecp')];
        Mol.Basis_Set_TD = ...
        [Basis_Set_file newline ...
         newline ...
         ECP_file newline];

        Multiplicities = 1;
    case 'Sc'
        Mol.Metal_Charge = 3;
        Mol.Basis_Route_Rough = 'Def2SVP';
        Mol.Basis_Set_Rough = '';
        Mol.Basis_Rough_Aux = 'W06';

        Mol.Basis_Route_Tight = 'Def2TZVP';
        Mol.Basis_Set_Tight = '';
        
        Mol.Basis_Route_TD = 'Gen';
        Mol.Basis_Set_TD = ['@' fullfile(Basis_Set_Loc,'def2-TZVPD.gbs') newline];

        Multiplicities = 1;
    otherwise
        error(['Unknown Metal Ion: ' Mol.Metal_Ion])
end

 
end