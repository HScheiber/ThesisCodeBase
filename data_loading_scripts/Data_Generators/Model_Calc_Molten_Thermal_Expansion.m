% Script for generating liquid densities as a function of temperature for a
% given force-field model. Requires an MD simulation that runs over a range of
% temperatures.

Salt = 'LiF';
Theory = 'JC';
Model = 'EQ3';
startpoint = '500'; % ps
endpoint = '19500'; % ps

DatDir = ['C:\Users\Hayden\Documents\Patey_Lab\Results\' Salt '\NucTest_C_' Theory '_Model_' Model '_NPT'];
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

en_opts = regexp(outpt,'End your selection with an empty line or a zero.\n-+(.+?)\n\n','tokens','once');
en_opts = en_opts{1};
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

% Fit an exponential

% modelfun = @(b,x)b(1) + b(2)*exp(b(3)*x)+ b(4)*exp(b(5)*x);
% beta0 = [5 1e-2 1e-4 0.5 1e-4];

% modelfun = @(b,x)b(1) + b(2)*exp(b(3)*x) + b(4)*exp(b(5)*x);
% beta0 = [2.3351 15.702 0.00031681 1e-3 1e-4];

modelfun = @(b,x)b(1) + b(2)*exp(b(3)*x);
beta0 = [2.3351 15.702 0.00031681];
opts = statset('Display','iter','TolFun',1e-20,'TolX',1e-20,'maxiter',1000);
mdl = fitnlm(T,Vol,modelfun,beta0,'Options',opts);

Tfit = (0:0.1:max(T))';
[Vfit,Vfitci] = predict(mdl,Tfit);

figure
hold on
plot(T,Vol,Tfit,Vfit)
plot(Tfit,Vfitci(:,1),':g')
plot(Tfit,Vfitci(:,2),':g')


filename = fullfile(savedir,[Salt '_L_' Theory '_' Model '.mat']);
save(filename,'mdl');
