function is_complete = cp2k_check_complete(cdir,Is_Salt)

is_complete = false; % initialize as false

d_inp = dir([cdir filesep '*.inp']);
if isempty(d_inp)
    return
end
inptxt = fileread(fullfile(cdir,d_inp(1).name));

% Check for DKH error
if isempty(regexp(cdir,'-DKH','once'))
    if ~isempty(regexp(inptxt,'&RELATIVISTIC','once'))
        disp('DKH was accidently turned on!');
        return
    end
end

% Check for fixed reference cell
if ~isempty(regexp(inptxt,'&CELL_REF','once'))
    disp('The reference cell was held fixed!');
    return
end

if Is_Salt
    % Check for derivative-free opt complete
    d = dir([cdir filesep 'CALCULATION_CONVERGED.mat']);
    if ~isempty(d)
        is_complete = true;
        return
    end
    
    % Find all log files in current folder
    d = dir([cdir filesep '*.out']);
    % Search each log file for results of interest
    for ldx = 1:length(d)

        if ~isempty(regexp(d(ldx).name,'slurm.*','ONCE')) || ...
                ~isempty(regexp(d(ldx).name,'8\.out','ONCE'))
            continue
        end

        jobtxt = fileread(fullfile(cdir,d(ldx).name));
        comptxt = regexp(jobtxt,'(GEOMETRY OPTIMIZATION COMPLETED|run CONVERGED!).+DBCSR STATISTICS','match','once');

        if isempty(comptxt)
            continue
        else
            is_complete = true;
        end
    end

else % Ions

    % Find all log files in current folder
    d = dir([cdir filesep '*.out']);

    % Search each log file for results of interest
    for kdx = 1:length(d)

        if ~isempty(regexp(d(kdx).name,'slurm.*','ONCE')) || ...
                ~isempty(regexp(d(kdx).name,'8\.out','ONCE'))
            continue
        end

        jobtxt = fileread(fullfile(cdir,d(kdx).name));
        Entxt = regexp(jobtxt,'ENERGY\| Total FORCE_EVAL \( QS \) energy .a\.u\..: +(-|\.|[0-9]|E)+','tokens','once');
        warntxt = regexp(jobtxt,'SCF run NOT converged','match','once');

        if isempty(Entxt)
            continue
        elseif ~isempty(warntxt)
            continue
        else
            is_complete = true;
        end
    end
end


end