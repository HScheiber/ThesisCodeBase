function Basis_out = cp2k_basis_set(Basis_Name,Atom,Increase_Radial_Grid)

% Possible basis sets:
% pob-TZVP
% def2-TZVPD
% X2C-TZVPall
% X2C-TZVPPall
% ANO-RCC-VTZP
% ANO-RCC-QTZP
% Jorge-DZP-DKH
% Jorge-TZP-DKH
% Jorge-QZP-DKH
% Sapporo-TZP
% Sapporo-TZP-diffuse
% Sapporo-TZP-2012
% Sapporo-TZP-2012-diffuse
% Sapporo-QZP
% Sapporo-QZP-diffuse
% Sapporo-QZP-2012
% Sapporo-QZP-2012-diffuse

if nargin <= 2
    Increase_Radial_Grid = 0;
end

switch Basis_Name
    case 'pob-TZVP'
        Basis_set = ['		BASIS_SET pob-TZVP' newline];
    otherwise
        if ispc
            Basis_Set_Loc = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\basis_sets';
        else
            [~,Server] = system('hostname');
            if ~isempty(regexp(Server(1:3),'se[0-9]','ONCE')) || strcmpi(Server(1:3),'log')
                Basis_Set_Loc = '/home/haydensc/ThesisCodeBase/basis_sets';
            elseif strncmpi(Server,'pat',3)
                Basis_Set_Loc = '/home/user/ThesisCodeBase/basis_sets';
            else
                Basis_Set_Loc = '/home/scheiber/ThesisCodeBase/basis_sets';
            end
        end
        Basis_Set_file = fullfile(Basis_Set_Loc,[Basis_Name '.cp2k']);
        bs_txt = fileread(Basis_Set_file);
        bs_atom = regexpi(bs_txt,['\n' Atom ' ' strrep(Basis_Name,'-DKH3','(-DKH3)*') '.*?(\n|\r){2}'],'match','once');
        Basis_set = regexprep(bs_atom(2:end-1),'.+?\n','','once');
        Basis_set = [...
            '		&BASIS' newline ...
            Basis_set ...
            '		&END BASIS' newline];
end

switch Atom
    case {'Br' 'I' 'At'}
        Radial_Grid = num2str(150 + Increase_Radial_Grid);
    otherwise
        Radial_Grid = num2str(80 + Increase_Radial_Grid);
end

Basis_out = [...
'    &KIND ' Atom newline ...
Basis_set ...
'		POTENTIAL All' newline ...
' 		LEBEDEV_GRID 50' newline ...
' 		RADIAL_GRID ' Radial_Grid newline ...
'    &END KIND'];



end