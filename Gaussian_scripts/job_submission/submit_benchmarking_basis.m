%% Inputs
version = 1; % 0 = g16; 1 = g09 ((BE SURE TO LOAD THE CORRECT MODULE))
hours = 3; % Max time for job (hours)
mins = 0; % Max time for job (minutes)
template_file = 'LiLi_Ion_Pair_B2GP-PLYP-D3BJ.template'; % Name of template file
label = 'LiLi'; % system label
pausetime = 2; % time to pause
nMols_per_Task = -1; % -1 to fix the number of cores
nCores = 32;
nTasks_per_Node = 32;
submit_type = 'submit_Gaussian_special.pl';

l1a = {'6-311+G*','CBSB7++','AUG-cc-pVDZ','AUG-cc-pVTZ','AUG-cc-pVQZ','Def2SVPP','Def2TZVPP','Def2QZVPP'}; % Basis sets
l1l = 'BASS';

%% Code
workDir = pwd;

for i = 1:length(l1a)

	l1s = l1a{i};
	l1d = [l1l '_' l1a{i}];
 
	% create/change dir
	if ~exist(l1d,'dir')
        mkdir(l1d);
	end
	dir1 = [workDir filesep l1d];
	cd(dir1)

	nTasks = num2str(nCores);
	hours_calc = num2str(hours);
	mins_calc = num2str(mins);
	nCores_per_Node = num2str(nTasks_per_Node);
	vers = num2str(version);

	% Open template
	fid = fopen(['/home/scheiber/templates' filesep template_file],'rt');
	X = fread(fid);
	fclose(fid);
	X = char(X.');

	% Replace strings
	Y = strrep(X,['##' l1l '##'], l1s);

	% Save file
	fid2 = fopen( [label '_' l1s '.com'],'wt');
	fwrite(fid2,Y);
	fclose(fid2);

	% Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe
	command = [submit_type ' -1 ' nTasks ' ' nCores_per_Node ...
	' -1 -1 ' hours_calc ' ' mins_calc ' ' label '_' l1s ...
	' ' label '_' l1s '.com ' label '_' l1s ' 1 ' vers];
	system(command);

	pause(pausetime);

	cd(workDir)
end
