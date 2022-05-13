% Short function to grab final converged energy
function Energy = GrabEnergyD4(Salt,OptDir,JobName,CryGeom,Theory,Dispersion)

[Metal,Halide] = Separate_Metal_Halide(Salt);

%% Generate transformation matrix from crystal basis to cartesian basis
abc = [CryGeom.a 0 0; 0 CryGeom.b 0; 0 0 CryGeom.c];
TM = abc*GenTransformMatrix(CryGeom);

%% Generate .vasp input file for dftd4
NF = CryGeom.NFc; % D4 works with the conventional cell as set up here  
geom_txt = ['New Structure' newline '1.0' newline];
for i = 1:3
    for j = 1:3
        geom_txt = [geom_txt pad(num2str(TM(i,j),'%10.10f'),20,'left')];
    end
    geom_txt = [geom_txt newline];
end
Nd2 = num2str(NF);

geom_txt = [geom_txt pad(Metal,5,'left') pad(Halide,5,'left') newline];
geom_txt = [geom_txt pad(Nd2,5,'left') pad(Nd2,5,'left') newline 'Direct' newline];

for i = 1:(NF)
    for j = 1:3
        geom_txt = [geom_txt pad(num2str(CryGeom.FC_Metal(i,j),'%10.10f'),16,'left')];
    end
    geom_txt = [geom_txt newline];
end
for i = 1:(NF)
    for j = 1:3
        geom_txt = [geom_txt pad(num2str(CryGeom.FC_Halide(i,j),'%10.10f'),16,'left')]; %#ok<*AGROW>
    end
    geom_txt = [geom_txt newline];
end

% save to file
filename1 = [OptDir filesep JobName '.vasp'];
fidPM = fopen(filename1,'wt');
fwrite(fidPM,regexprep(geom_txt,{'\r', '\n\n+'}',{'', '\n'}));
fclose(fidPM);

% run dftd4 on the structure
if ispc % for testing
    dftd4 = ['wsl source ~/.bashrc; dftd4 -f ' lower(Theory) ' '];
    fnunix = windows2unix(filename1);
elseif isunix
    srv = getenv('USERNAME');
    if strcmpi(srv,'user')
        dftd4 = ['source ~/.bashrc; dftd4_24 -f ' lower(Theory) ' '];
        fnunix = filename1;
    else
        dftd4 = ['dftd4 -f ' lower(Theory) ' '];
        fnunix = filename1;
    end
end

if ~strcmp(Dispersion,'D4TB') % two body terms only
    dftd4 = [dftd4 '-2 '];
end

[ercode,dftd4_out] = system([dftd4 fnunix]);
delete(filename1)

if ercode ~= 0 
    error(['DFTD4 Program failed. Printing output:' newline dftd4_out])
end

% Grab output energy file .EDISP
filename2 = [OptDir filesep '.EDISP'];
fidPM = fopen(filename2,'rt');
En = textscan(fidPM,'%f');
fclose(fidPM);
delete(filename2);

if isempty(En)
    error(['DFTD4 Program failed. Printing output:' newline dftd4_out])
end

Energy = En{1}/(NF); % a.u. per formula unit

end