%% Generates WURTZITE or ROCKSALT crystal structures
% NOTE: In Gaussian, empirical dispersion corrections have been disabled for PBC
%% Inputs and options
hours = 3; % Max time for job (hours)
mins = 0; % Max time for job (minutes)

salt_type = 'LiF';
crystal_type = 'Wurtzite'; 
template_file = 'LiHalide_Crystal.template'; % Name of template file
nMols_per_Task = -1; % -1 to fix the number of cores
nCores = 32;
nTasks_per_Node = 32;
mempernode = '-1'; % Memory request for server (default = '-1')
dispersion = ''; % Set to 'GD3BJ' / 'GD2' / 'PFD' / 'GD3' / '' (no dispersion);
conv = 8; % Convergence criteria
kpoint = 1000; % Number of K-points included in integration (0 FOR GAMMMA POINT ONLY, NEGATIVE FOR DEFAULT)
range = 100; % Integration range in bohr radii (negative number to switch to default)
memory = '64000MB'; % Memory request for Gaussian

% Theory. Note: B2PLYPD3, B97D3,'wB97XD' contain empirical dispersion and cannot be used
%theory = {'LSDA','PBEPBE','PBE1PBE','DSDPBEP86','B2GPPLYP','HSEH1PBE','wB97X','wB97','LC-B3LYP',...
%    'CAM-B3LYP','LC-wHPBE','B2PLYP','PW91PW91','mPW1PW91','PW1PW','TPSSTPSS','RevTPSSRevTPSS','HISSbPBE'}; 
%theory = {'HF','B3LYP','BLYP','BPW91','SLYP','SPW91','LSDA'};
theory = {'LSDA'};
% Basis set
%basis = {'pob-TZVP','7-311G-jansen','7-311G-crystal','7-311G-nada'};
basis = {'pob-TZVP'};

% Lattice parameter a
a = 2:0.1:8;

%% Wurtzite unit cell for LiF
if strcmp('Wurtzite',crystal_type)

    Coordinates = [0.000000000000000    0.577350269189626    0.612372435695794;...
                   0.000000000000000    0.577350269189626    0.000000000000000;...
                   0.500000000000000    0.288675134594813    1.428869016623521;...
                   0.500000000000000    0.288675134594813    0.816496580927726;...
                   1.000000000000000    0.000000000000000    0.000000000000000;...
                  -0.500000000000000    0.866025403784439    0.000000000000000;...
                   0.000000000000000    0.000000000000000    1.632993161855452];
               
    if strcmp('LiF',salt_type)
        Setting = {'Li'; 'F'; 'Li'; 'F'; 'Tv'; 'Tv'; 'Tv'};
        label = 'LiF_W';
    elseif strcmp('LiCl',salt_type)
        Setting = {'Li'; 'Cl'; 'Li'; 'Cl'; 'Tv'; 'Tv'; 'Tv'};
        label = 'LiCl_W';
    elseif strcmp('LiBr',salt_type)
        Setting = {'Li'; 'Br'; 'Li'; 'Br'; 'Tv'; 'Tv'; 'Tv'};
        label = 'LiBr_W';
    elseif strcmp('LiI',salt_type)
        Setting = {'Li'; 'I'; 'Li'; 'I'; 'Tv'; 'Tv'; 'Tv'};
        label = 'LiI_W';
    else
        error('Error: Unknown Salt Type');
    end    
               
elseif strcmp('Rocksalt',crystal_type)
    
    Coordinates = [0.00000000 0.00000000 0.00000000;...
                   0.50000000 0.00000000 0.00000000;...
                   0.00000000 0.50000000 0.50000000;...
                   0.50000000 0.00000000 0.50000000;...
                   0.50000000 0.50000000 0.00000000;...
                   0.00000000 0.50000000 0.00000000;...
                   0.00000000 0.00000000 0.50000000;...
                   0.50000000 0.50000000 0.50000000;...
                   1.00000000 0.00000000 0.00000000;...
                   0.00000000 1.00000000 0.00000000;...
                   0.00000000 0.00000000 1.00000000];

    if strcmp('LiF',salt_type)
        Setting = {'F'; 'Li'; 'F'; 'F'; 'F'; 'Li'; 'Li'; 'Li'; 'Tv'; 'Tv'; 'Tv'};
        label = 'LiF_R';
    elseif strcmp('LiCl',salt_type)
        Setting = {'Cl'; 'Li'; 'Cl'; 'Cl'; 'Cl'; 'Li'; 'Li'; 'Li'; 'Tv'; 'Tv'; 'Tv'};
        label = 'LiCl_R';
    elseif strcmp('LiBr',salt_type)
        Setting = {'Br'; 'Li'; 'Br'; 'Br'; 'Br'; 'Li'; 'Li'; 'Li'; 'Tv'; 'Tv'; 'Tv'};
        label = 'LiBr_R';
    elseif strcmp('LiI',salt_type)
        Setting = {'I'; 'Li'; 'I'; 'I'; 'I'; 'Li'; 'Li'; 'Li'; 'Tv'; 'Tv'; 'Tv'};
        label = 'LiI_R';
    else
        error('Error: Unknown Salt Type');
    end
    
else
    error('Error: Unknown Crystal Structure Type');
    
end
% For effective core potentials
if strcmp('LiI',salt_type)
    gen = 'GenECP';
else
    gen = 'Gen';
end

%% Code

% Autochoose home folder based on system type
if ispc
    home = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase'; % PC
elseif isunix
    home = '/home/scheiber/ThesisCodeBase'; % Cedar/Graham
else
    home = input('Please input home bin directory.\n','s');
end

% For multi-node processes, invoke Linda
if nCores/nTasks_per_Node > 1
    linda=['%LindaWorkers=##NODES##' newline '%DebugLinda' newline];
    nCores = nTasks_per_Node;
else
    linda='';
end

ECPlink = newline;
workDir = pwd;
nTasks = num2str(nCores);
hours_calc = num2str(hours);
mins_calc = num2str(mins);
nCores_per_Node = num2str(nTasks_per_Node);
submit_Gaussian = 'submit_Gaussian_array.pl';

for i = 1:length(theory)
    
    l1s = theory{i};  
    l1d = [label '_' regexprep(theory{i},'\((.?)\)','_$1')];
    
    % Set additional environment variables option for submission script
    if strcmp(l1s,'DSDPBEP86') && strcmp(dispersion,'GD3BJ')
        Extra = '0';
    elseif strcmp(l1s,'B2GPPLYP') && strcmp(dispersion,'GD3BJ')
        Extra = '1';
    else
        Extra = '-1';
    end
    
    % create/change dir
    if ~exist(l1d,'dir')
        mkdir(l1d);
    end
    
    dir1 = [workDir filesep l1d];
    cd(dir1)

    for j = 1:length(basis)

        l2s = basis{j};
        l2d = ['BASS_' basis{j}];

        % create/change dir
        if ~exist(l2d,'dir')
            mkdir(l2d);
        end
        dir2 = [workDir filesep l1d filesep l2d];
        cd(dir2)

        for k=1:length(a)
            Coord_Scaled = Coordinates.*a(k);

            coord_cell = cellfun(@(x) num2str(x,'%15.15f'),...
                num2cell(Coord_Scaled),'UniformOutput',false);
            coord_cell2 = [Setting coord_cell];
            l3s = strjoin(join(coord_cell2,'    '), '\n');
            
            l3d = ['APAR_' num2str(a(k),'%05.2f')];
            
            % create/change dir
            if ~exist(l3d,'dir')
                mkdir(l3d);
            end
            dir3 = [workDir filesep l1d filesep l2d filesep l3d];
            cd(dir3)
        
            % Open template
            fid = fopen([home filesep 'templates' filesep template_file],'rt');
            X = fread(fid);
            fclose(fid);
            newtext = char(X.');

            % Replace strings
            % Theory
            if strcmp(l1s,'B2GPPLYP')
                newtext = strrep(newtext,'##THERT##','B2GP-PLYP');
                newtext = strrep(newtext,'##THER##','b2plyp');
                newtext = strrep(newtext,'##IOPS##',...
                    [newline 'iop(3/125=0360003600,3/78=0640006400,3/76=0350006500,3/77=1000010000)']);
                if ~strcmp(dispersion,'GD3BJ') && ~strcmp(dispersion,'')
                    error(['Unknown parameters for dispersion type ' ...
                        dispersion ' for ' l1s]);
                end
            elseif strcmp(l1s,'DSDPBEP86')
                newtext = strrep(newtext,'##THERT##','DSD-PBEP86');
                newtext = strrep(newtext,'##THER##','b2plyp');
                newtext = strrep(newtext,'##IOPS##',...
                    [newline 'iop(3/125=0220005200,3/78=0440004400,3/76=0310006900,3/74=1004)']);
                if ~strcmp(dispersion,'GD3BJ') && ~strcmp(dispersion,'')
                    error(['Unknown parameters for dispersion type ' ...
                        dispersion ' for ' l1s]);
                end
            elseif strcmp(l1s,'PW1PW')
                newtext = strrep(newtext,'##THERT##','PW1PW');
                newtext = strrep(newtext,'##THER##','PW91PW91');
                newtext = strrep(newtext,'##IOPS##',...
                    [newline 'iop(3/76=1000002000,3/77=0800000000,3/78=1000000000)']);
                if ~strcmp(dispersion,'')
                    error(['Unknown parameters for dispersion type ' ...
                        dispersion ' for ' l1s]);
                end
            else
                newtext = strrep(newtext,'##THERT##',l1s);
                newtext = strrep(newtext,'##THER##',l1s);
                newtext = strrep(newtext,'##IOPS##','');
            end

            % Basis set
            if strcmp(l2s,'pob-TZVP')
                newtext = strrep(newtext,'##BASS##',gen);
                newtext = strrep(newtext,'##BASST##','pob-TZVP');
                % Custom basis sets
                newtext = strrep(newtext,'##GENBAS##',['@' home filesep ...
                    'basis_sets' filesep 'pob_TZVP.gbs' newline newline]);
                if strcmp('LiI',salt_type)
                    ECPlink = ['@' home filesep 'basis_sets' filesep ...
                        'I_ECP_pob_TZVP.ecp' newline newline];
                end
            elseif strcmp(l2s,'7-311G-crystal')
                newtext = strrep(newtext,'##BASS##',gen);
                newtext = strrep(newtext,'##BASST##','7-311G-crystal');
                % Custom basis sets
                newtext = strrep(newtext,'##GENBAS##',['@' home filesep ...
                    'basis_sets' filesep '7_311G_crystal.gbs' newline newline]);
            elseif strcmp(l2s,'7-311G-jansen')
                newtext = strrep(newtext,'##BASS##',gen);
                newtext = strrep(newtext,'##BASST##','7-311G-jansen');
                % Custom basis sets
                newtext = strrep(newtext,'##GENBAS##',['@' home filesep ...
                    'basis_sets' filesep '7_311G_jansen.gbs' newline newline]);
            elseif strcmp(l2s,'7-311G-nada')
                newtext = strrep(newtext,'##BASS##',gen);
                newtext = strrep(newtext,'##BASST##','7-311G-nada');
                % Custom basis sets
                newtext = strrep(newtext,'##GENBAS##',['@' home filesep ...
                    'basis_sets' filesep '7_311G_nada.gbs' newline newline]);
            else
                newtext = strrep(newtext,'##BASS##', l2s);
                newtext = strrep(newtext,'##BASST##', l2s);
                newtext = strrep(newtext,'##GENBAS##','');
            end
            
            newtext = strrep(newtext,'##APAR##', l3s);
            newtext = strrep(newtext,'##CONV##', num2str(conv));
            newtext = strrep(newtext,'##MEMORY##', num2str(memory));
            newtext = strrep(newtext,'##CRYS##', crystal_type);
            newtext = strrep(newtext,'##APART##',num2str(a(k),'%05.2f'));
			newtext = strrep(newtext,'##CORES##',nTasks);
            newtext = strrep(newtext,'##CTYPE##',salt_type);
            newtext = strrep(newtext,['##LINDA##' newline],linda);
            newtext = strrep(newtext,'##GENECP##',ECPlink);
            
            % Dispersion correction
            if isempty(dispersion)
                newtext = strrep(newtext,'##DISPT##','');
                newtext = strrep(newtext,'##DISP##','');
            else
                newtext = strrep(newtext,'##DISP##',...
                    ['EmpiricalDispersion=' dispersion]);
                newtext = strrep(newtext,'##DISPT##',...
                    ['and ' dispersion ' empirical dispersion correction']);
            end
            
			% Cell Range (negative number to switch off)
			if range < 0
                newtext = strrep(newtext,',##CELLR##','');
			else
				newtext = strrep(newtext,'##CELLR##',...
					['CellRange=' num2str(range)]);
			end

            % K points calculation options
            if kpoint == 0
                newtext = strrep(newtext,'##KPOINT##',...
                    'GammaOnly');
            elseif kpoint < 0
                newtext = strrep(newtext,'##KPOINT##,', '');
				newtext = strrep(newtext,'##KPOINT##', '');
            else
                newtext = strrep(newtext,'##KPOINT##', ['NKpoint=' num2str(kpoint)]);
            end
			
			% Catch empty case
			newtext = strrep(newtext,' PBC=()', '');

            % Save file
            fid2 = fopen( [l1d '_' l2s '-' num2str(a(k),'%05.2f') '.com'],'wt');
            fwrite(fid2,newtext);
            fclose(fid2);
            
            cd(dir2)
        end
        
        params = regexprep(num2str(a.*100),' +',',');
        
        % Submit Array job
        % Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe
        command = [submit_Gaussian ' -1 ' nTasks ' ' nCores_per_Node ...
        ' -1 ' mempernode ' ' hours_calc ' ' mins_calc ' ' l1d '_' l2s ...
        ' ' l1d '_' l2s ' ' l1d '_' l2s ' 1 '  params ' ' Extra];
        disp([l1d ' ' l2s ' ' num2str(min(a)) ' - ' num2str(max(a))])
        if ~ispc
            system(command);
        else
            disp(command); % for debugging
        end
        
        cd(dir1)
    end
    cd(workDir)
end
