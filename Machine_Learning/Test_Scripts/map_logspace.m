function out = map_logspace(min_exp,y1,y2,N)
    logmap = logspace(min_exp,0,N);
    diff_log = abs(logmap(1) - logmap(end));
    logmap = logmap./diff_log;
    logmap = logmap - logmap(1);
    diff_minmax = y2 - y1;
    logmap = (logmap).*diff_minmax;
    out = logmap + y1;
end