function Sort_Atoms(filename)
fid = fopen(filename,'rt');

header = textscan(fid, '%s', 2,'Delimiter','\n');

% Header text
headertxt = [header{1}{1} newline ' ' header{1}{2}];
totmoles = str2double(header{1}{2});

Input_data = textscan(fid,'%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n',...
    'Delimiter','','whitespace','');

endidx = length(Input_data{1,1});
for i=1:length(Input_data)
    if length(Input_data{1,i}) == endidx
        Input_data{1,i}(endidx) = [];
    end
end

% Get end cell parameter text
frewind(fid);
finaltext = textscan(fid,'%s','delimiter','\n','headerlines',endidx+1,...
    'whitespace','');
boxcoords = finaltext{1}{1};

% Convert elements to alphabetical order
Input_data{2} = strrep(Input_data{2},'Li','A#');
Input_data{2} = strrep(Input_data{2},'Na','B#');
Input_data{2} = strrep(Input_data{2},'K','C#');
Input_data{2} = strrep(Input_data{2},'Rb','D#');
Input_data{2} = strrep(Input_data{2},'F','E#');
Input_data{2} = strrep(Input_data{2},'Cl','F#');
Input_data{2} = strrep(Input_data{2},'Br','G#');
Input_data{2} = strrep(Input_data{2},'I','H#');

% Sort by atom type
[B,Indx] = sortrows([cellstr(num2str(Input_data{1})) Input_data{2}],[2,1]);
Molnum = B(:,1);
Mols_sorted = B(:,2);
Atoms_sorted = Input_data{3}(Indx);
Atomnum = 1:1:totmoles;
X_Coords = Input_data{5}(Indx);
Y_Coords = Input_data{6}(Indx);
Z_Coords = Input_data{7}(Indx);
P1 = Input_data{8}(Indx);
P2 = Input_data{9}(Indx);
P3 = Input_data{10}(Indx);

% Convert elements to back
Mols_sorted = strrep(Mols_sorted,'A#','Li');
Mols_sorted = strrep(Mols_sorted,'B#','Na');
Mols_sorted = strrep(Mols_sorted,'C#','K');
Mols_sorted = strrep(Mols_sorted,'D#','Rb');
Mols_sorted = strrep(Mols_sorted,'E#','F');
Mols_sorted = strrep(Mols_sorted,'F#','Cl');
Mols_sorted = strrep(Mols_sorted,'G#','Br');
Mols_sorted = strrep(Mols_sorted,'H#','I');

% Create new labels
M = length(Mols_sorted);
NumLab = (1:M)';
Outstring = headertxt;
for i=1:M
    str = cell(9,1);
    str{1} = Molnum{i};
    str{2} = Mols_sorted{i};
    str{3} = Atoms_sorted{i};
    str{4} = pad(num2str(Atomnum(i)),5,'left');
    str{5} = pad(num2str(X_Coords(i),'%8.3f'),8,'left');
    str{6} = pad(num2str(Y_Coords(i),'%8.3f'),8,'left');
    str{7} = pad(num2str(Z_Coords(i),'%8.3f'),8,'left');
    str{8} = pad(num2str(P1(i),'%8.4f'),8,'left');
    str{9} = pad(num2str(P2(i),'%8.4f'),8,'left');
    str{10} = pad(num2str(P3(i),'%8.4f'),8,'left');
    
    current_line = [str{:}];
    
    Outstring = [Outstring newline current_line];
end

% Final line
Outstring = [Outstring newline boxcoords newline];

% Save file
fid2 = fopen(filename,'wt');
fwrite(fid2,Outstring);
fclose(fid2);
end

