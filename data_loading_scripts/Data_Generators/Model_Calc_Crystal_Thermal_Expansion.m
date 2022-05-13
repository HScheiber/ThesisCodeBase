% Script for generating linear thermal expansion coefficients for solid
% alkali halides given a simulation of a model.
% Simulation should be a solid that is heated up to it's melting point.

Salt = 'LiF';
Theory = 'JC';
Model = 'ER2';
startpoint = '0'; % ps
endpoint = '19100'; % ps

DatDir = ['C:\Users\Hayden\Documents\Patey_Lab\Results\' Salt '\RefBox_R_' Theory '_Model_' Model '_NPT'];
savedir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\ThermExp';

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
En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Temperature','tokens','once'))];
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
% remove temp output file
delete(Out_file) 
rmdir(tname);
            
Time = Data(:,1)./1000; % ns
T = Data(:,2); % K
Vol = Data(:,3).*(10^3)./NF; % Angstrom^3/formula unit
Vol_UC = Vol*4; % Volume (Angstrom^3) per unit cell
a = Vol_UC.^(1/3);

% Fit an exponential

% modelfun = @(b,x)b(1) + b(2)*exp(b(3)*x)+ b(4)*exp(b(5)*x);
% beta0 = [5 1e-2 1e-4 0.5 1e-4];

modelfun = @(b,x)b(1) + b(2)*exp(b(3)*x) + b(4)*exp(b(5)*x);
beta0 = [4 1e-2 1e-4 1e-2 1e-4];
opts = statset('Display','iter','TolFun',1e-20,'TolX',1e-20,'maxiter',1000);
mdl = fitnlm(T,a,modelfun,beta0,'Options',opts);

Tfit = (0:0.1:max(T))';
[afit,afitci] = predict(mdl,Tfit);

figure
hold on
plot(T,a,Tfit,afit)
plot(Tfit,afitci(:,1),':g')
plot(Tfit,afitci(:,2),':g')


filename = fullfile(savedir,[Salt '_' Theory '_' Model '.mat']);
save(filename,'mdl');
