% Open and read the content
Myfile = 'geom.txt';
fid = fopen(Myfile);
Text = textscan(fid,'%s','delimiter','\n');
Text = Text{1};
fclose(fid);


% Analyze
ResName_mol = 'UNK';
ResName_solv = 'TIP';
Zero_pres = '0';
L_H = ' (L|H)';
H_in_Mol = '';
H_in_Solv = '';
C_in_Mol = '';
O_in_Mol = '';
O_in_Solv = '';
All_in_Mol = '';
All_in_Solv = '';
for i = 1:length(Text)
    Current_Line = Text{i};
    if ~isempty(regexp(Current_Line,'PDBName=H,','ONCE'))
        H_in_Mol = [H_in_Mol num2str(i) ' '];
    end
    if ~isempty(regexp(Current_Line,'PDBName=H(1|2),','ONCE'))
        H_in_Solv = [H_in_Solv num2str(i) ' '];
        Line_Out = regexprep(Current_Line,['H\(.+?\) +' Zero_pres '(.+?)' L_H],[pad('H2(Iso=2,Spin=2,NMagM=0.857438228)',36) '$1']);
    end
    if ~isempty(regexp(Current_Line,'PDBName=C,','ONCE'))
        C_in_Mol = [C_in_Mol num2str(i) ' '];
    end
    if ~isempty(regexp(Current_Line,'PDBName=O,','ONCE'))
        O_in_Mol = [O_in_Mol num2str(i) ' '];
    end
    if ~isempty(regexp(Current_Line,'PDBName=OH2,','ONCE'))
        O_in_Solv = [O_in_Solv num2str(i) ' '];
        Line_Out = regexprep(Current_Line,['O\(.+?\) +' Zero_pres '(.+?)' L_H],[pad('O2',36) '$1']);
    end
    if ~isempty(regexp(Current_Line,['ResName=' ResName_mol ','],'ONCE'))
        All_in_Mol = [All_in_Mol num2str(i) ' '];
        Line_Out = regexprep(Current_Line,['([A-Z][a-z]*)\(.+?\) +' Zero_pres '(.+?)' L_H],[pad('$1$`1',39) '$2']);
    end
    if ~isempty(regexp(Current_Line,['ResName=' ResName_solv ','],'ONCE'))
        All_in_Solv = [All_in_Solv num2str(i) ' '];
    end
    
    disp(Line_Out)
    
end
% result
fprintf('Atoms in molecule: \t%s 0\n',strtrim(All_in_Mol))
fprintf('Atoms in solvent: \t%s 0\n',strtrim(All_in_Solv))
fprintf('H in molecule: \t%s 0\n',strtrim(H_in_Mol))
fprintf('H in solvent: \t%s 0\n',strtrim(H_in_Solv))
fprintf('C in molecule: \t%s 0\n',strtrim(C_in_Mol))
fprintf('O in molecule: \t%s 0\n',strtrim(O_in_Mol))
fprintf('O in solvent: \t%s 0\n',strtrim(O_in_Solv))
