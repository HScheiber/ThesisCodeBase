% Short function to grab final converged energy
% N is the number of formula units per unit cell
function Energy = GrabEnergyDFT(OptDir,JobName,NF)

% First find the log file
logfile_info = dir([OptDir filesep JobName '.out']);

if isempty(logfile_info)
    error(['No log file found for: ' JobName]);
    
% If multiple log files exist, choose most recent
elseif size(logfile_info,1) > 1
    Times = zeros(1,size(logfile_info,1));
    for y = 1:size(logfile_info,1)
        Times(y) = datenum(logfile_info(y).date,...
            'dd-mmm-yyyy HH:MM:SS');
    end
    [~,Idx] = max(Times);
    logfile_info = logfile_info(Idx);
end

% Otherwise, get filename
logfile = [OptDir filesep logfile_info.name];

% Open the file
fid = fopen(logfile,'rt');
logtext = fread(fid);
fclose(fid);
logtext = char(logtext.');

if isempty(logtext)
    error(['No log file found for: ' OptDir filesep JobName '.out']);
end

% Find energy outputs in hartree, abort if not found        
energy = regexp(logtext,...
    '(CONVERGENCE ON (ENERGY|TESTER) +E\(AU\) +)(-|\+|\.|[0-9]|E|)+','tokens');

if isempty(energy) 
    error(['No energy found in: ' logfile]);
end

energy_num = str2double(energy{1}{2});

% Get dispersion energy if available
dispersiontxt = regexp(logtext,...
    '(D3 DISPERSION ENERGY \(AU\) +)(-|\+|\.|[0-9]|E|)+','tokens');

% Get gCP energy if available
gCPtxt = regexp(logtext,...
    '(GCP ENERGY \(AU\) +)(-|\+|\.|[0-9]|E|)+',...
    'tokens');

if ~isempty(gCPtxt) && ~isempty(dispersiontxt)
    disp_E = str2double(dispersiontxt{1}{2});
    gCP_E = str2double(gCPtxt{1}{2});
elseif ~isempty(dispersiontxt)
    disp_E = str2double(dispersiontxt{1}{2});
    gCP_E = 0;
elseif ~isempty(gCPtxt)
    disp_E = 0;
    gCP_E = str2double(gCPtxt{1}{2});
else
    disp_E = 0;
    gCP_E = 0;
end
Energy = (energy_num + disp_E + gCP_E)/NF; % Units of Hartree

end