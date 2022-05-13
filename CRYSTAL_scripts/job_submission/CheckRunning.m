function Response = CheckRunning(Current_Directory,Stop_Running)

Response = false;

% Check for log files
Logfiles = dir([Current_Directory filesep '*.out']);

if isempty(Logfiles)
    return
end

% Check the logfile for running jobs
for i = 1:length(Logfiles)
    Logfile = [Current_Directory filesep Logfiles(i).name];
    
    fid = fopen(Logfile,'rt');
    Logtext = fread(fid);
    fclose(fid);
    Logtext = char(Logtext.');
    
    tempdir = regexp(Logtext,'temporary directory (.+?)\n','Tokens','ONCE');
    
    if isempty(tempdir)
        continue
    end
    
    jobID = regexp(tempdir{1},'([0-9]+)\.orca[1-2]\.ibb','Tokens','ONCE');
    
    if isempty(jobID)
        continue
    end
    
    % Check for job in queue
    [~,Output] = system(['qstat -u $USER | grep ' jobID{1}]);
    
    if isempty(Output)
        continue
    else
        if Stop_Running
            disp(['A Requested Job is already running with Job ID: ' jobID{1} ' - Ending this Job and starting a new one!'])
            [Ercode,DelOut] = system(['qdel ' jobID{1}]);
            if Ercode ~= 0
                error(['Unable to end job: ' jobID{1} '. Err: ' DelOut])
            end
            return
        else
            Response = true;
            disp(['A Requested Job is already running with Job ID: ' jobID{1} ' - Skipping!'])
            return
        end
    end
end
end