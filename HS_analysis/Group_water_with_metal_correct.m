fid = fopen('ammonium.gro','rt');

Input_data = textscan(fid,'%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n',...
        'Delimiter','','whitespace','','headerlines',2);
fclose(fid);

Atoms = length(Input_data{1});
Mols = Atoms/5;

Oxygen_coordinates = zeros(Mols,3);
Hydrogen_coordinates = zeros(Mols*2,3);
j=1;
k=1;
for i=1:Atoms
    if strcmp( strtrim(Input_data{3}{i}),'OW' )
        Oxygen_coordinates(j,:) = [Input_data{5}(i) Input_data{6}(i) Input_data{7}(i)];
        j=j+1;
    else
        Hydrogen_coordinates(k,:) = [Input_data{5}(i) Input_data{6}(i) Input_data{7}(i)];
        k=k+1;
    end
    
end

a = 0.13458335;
b = a;

output = ['h2 o; 1Created by VESTA' newline ...
    num2str(Mols*4) newline];

atom_index = 1;
mol_index = 1;
f = waitbar(0,'Calculating...');
for i=1:Mols
    Current_oxygen = Oxygen_coordinates(i,:);
    
    newtext = sprintf('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n',mol_index,'SOL','OW',...
        atom_index,Current_oxygen(1),Current_oxygen(2),Current_oxygen(3));
    
    output = [output newtext];
    
    atom_index = atom_index + 1;
    
    k = 1;
    for j=1:2*Mols
        
        Current_hydrogen = Hydrogen_coordinates(j,:);
        
        Distance = norm(Current_oxygen - Current_hydrogen);
        
        if Distance < 0.13
            
            if k == 1
                First_Hydrogen = Current_hydrogen;
            end
            
            HW_text = ['HW' num2str(k)];
            newtext = sprintf('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n',mol_index,'SOL',HW_text,...
                atom_index,Current_hydrogen(1),Current_hydrogen(2),Current_hydrogen(3));
            output = [output newtext];
            atom_index = atom_index + 1;
            
            
            k=k+1;
        end
    end
    
    mol_index = mol_index + 1;
    waitbar(i/Mols,f);
end
close(f)
fid = fopen('boxOutbig.gro','wt');
fprintf(fid, output);
fclose(fid);