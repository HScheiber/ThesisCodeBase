function [Geometry,Total_Charge] = Get_Complex_Structure(Mol)

% Required properties of Mol (input):
% Metal_Ion
% Metal_Charge
% Structure

% Output:
% Geometry is a XYZ map describing the system geometry
% Total charge is the charge of the systme

if ispc
    Template_Dir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\submission\Gaussian_scripts\Metal_Geometries';
else
    Template_Dir = '/home/scheiber/ThesisCodeBase/submission/Gaussian_scripts/Metal_Geometries';
end

Protons = regexp(Mol.Structure,'H[0-9]*','match','once');
if isempty(Protons)
    Ligand_Charge = -4;
elseif strcmp(Protons,'H')
    Ligand_Charge = -3;
else
    Ligand_Charge = -4 + str2double(Protons(2:end));
end

try
    Geometry = fileread( fullfile(Template_Dir,[Mol.Structure '.zmat']) );
catch
    error(['Unknown Structure: ' Mol.Structure])
end
Geometry = strrep(Geometry,'##METAL##',Mol.Metal_Ion);
Total_Charge = Ligand_Charge + Mol.Metal_Charge;
 
end