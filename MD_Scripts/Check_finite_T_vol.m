
Salt = 'LiI';
Model = 'CM3';
startpoint = '12950'; % ps
endpoint = '13050'; % ps

DatDir = ['C:\Users\Hayden\Documents\Patey_Lab\Results\' Salt '\HiAcc_MeltNucCycle_R_JC_Model_' Model '_NPT'];


% Find energy file
Ener_info = dir([DatDir filesep '*.edr']);
Tpr_info = dir([DatDir filesep '*.tpr']);
Mat_info = dir([DatDir filesep '*.mat']);

tname = tempname;
mkdir(tname);

Ener_file = [DatDir filesep Ener_info.name]; % traj file
Tpr_file = [DatDir filesep Tpr_info.name]; % tpr file
Mat_file = [DatDir filesep Mat_info.name]; % mat file
Out_file = [tname filesep 'tempEnergy.xvg'];

% Check energy options
gmx_command = ['wsl source ~/.bashrc; echo 0 ^| gmx_d energy -f ' windows2unix(Ener_file)...
' -s ' windows2unix(Tpr_file)];
[err,outpt] = system(gmx_command);

if err ~= 1
    warndlg('Problem with energy check.')
    return
end

en_opts = regexp(outpt,'-+\n.+?-+\n','match','once');
En_set = '';
%En_set = char(regexp(en_opts,'([0-9]{1,2})  Potential','tokens','once'));
%En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Kinetic-En.','tokens','once'))];
%En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Total-Energy','tokens','once'))];
%En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Temperature','tokens','once'))];
%En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Pressure','tokens','once'))];
En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Volume','tokens','once'))];
%En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Density','tokens','once'))];
%En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Enthalpy','tokens','once'))];
En_set = [En_set ' 0'];
En_set = regexprep(En_set,' +',' ');

% Grab data from results
gmx_command = ['wsl source ~/.bashrc; echo ' En_set ' ^| gmx_d energy -f ' windows2unix(Ener_file)...
' -o ' windows2unix(Out_file) ' -s ' windows2unix(Tpr_file) ...
' -b ' startpoint ' -e ' endpoint];

err = system(gmx_command,'-echo');

if err ~= 0
    warndlg('Failed to collect data.')
    return
end

NSS = load(Mat_file);
NF = NSS.N_total/2; % number of ion pairs in supercell

Data = import_xvg(Out_file); % Gather RDF
delete(Out_file) % remove temp output file
            
Time = Data(:,1)./1000; % ns
Vol = mean(Data(:,2).*(10^3)./NF); % Angstrom^3/formula unit
disp(['Volume: ' num2str(Vol,'%.15f')])

Vol_UC = Vol*4; % Volume (Angstrom^3) per unit cell
a = Vol_UC^(1/3);
disp(['a: ' num2str(a,'%.15f')])
