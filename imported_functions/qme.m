function qme
    [~,~,computer,~] = find_home;
    switch lower(computer)
        case 'sockeye'
            [~,cmdout] = system('qstat -u $USER | tee /dev/stderr | grep haydensc | wc -l | sed -E ''s/^/Total Jobs: /''');
        case {'cedar' 'graham' 'narval'}
            [~,cmdout] = system('squeue -u $USER -o "%10i %35j %8T %10M %10l %5D %5C %10m %11S %20a %60R " | tee /dev/stderr | wc -l | xargs printf "%d-1\n" |  bc -l | sed -E "s/^/Total Jobs: /"');
        otherwise
            [~,cmdout] = system('qstat -u $USER');
    end
    disp(cmdout)
end