function HoursRemain = JobTimeRemaining
    err = 1;
    timechunks = {};
    while err ~= 0 || isempty(timechunks)
        [err,timeleft] = system('squeue -h -j $SLURM_JOBID -o %L');
        timeleft = strtrim(timeleft);
        timechunks = regexp(timeleft,'([0-9]{1,2}-)*([0-9]{1,2}:)*([0-9]{1,2}):([0-9]{2})','tokens','once');
    end
    
    if isempty(timechunks{1})
        Days = '00';
    else
        Days = sprintf('%02s', strrep(timechunks{1},'-',''));
    end
    if isempty(timechunks{2})
        Hours = '00';
    else
        Hours = sprintf('%02s', strrep(timechunks{2},':',''));
    end
    Minutes = timechunks{3};
    Seconds = timechunks{4};
    TimeString = [Days ':' Hours ':' Minutes ':' Seconds];
    HoursRemain = hours(duration(TimeString,'InputFormat','dd:hh:mm:ss'));
end
