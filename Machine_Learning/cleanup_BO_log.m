function cleanup_BO_log(logfile)

logtxt = fileread(logfile);

logtxt = regexprep(logtxt,'MATLAB is selecting.+?matlab:opengl.+?\n','');
logtxt = regexprep(logtxt,'Opening log file.+?\n','');
logtxt = regexprep(logtxt,'\n+','\n');
logtxt = regexprep(logtxt,'Starting parallel pool.+?\n','');
logtxt = regexprep(logtxt,'Connected to the parallel pool.+?\n','');
logtxt = regexprep(logtxt,'Copying objective function.+?\n','');
logtxt = regexprep(logtxt,'Done copying objective.+?\n','');
logtxt = regexprep(logtxt,'\s+< M A T L A B \(R\).+?mathworks\.com\.\n','');
logtxt = regexprep(logtxt,'MATLAB is selecting.+?\n','');
logtxt = regexprep(logtxt,'connected to [0-9]+ workers.+?\n','');
logtxt = regexprep(logtxt,'Your MATLAB session has timed out.+?\n','');
logtxt = regexprep(logtxt,'Parallel pool using the.+?\n','');
paramsummary = regexp(logtxt,'\*+\nModel Input Parameters.+?\*+\n\*+\n','match');
logtxt = regexprep(logtxt,'\*+\nModel Input Parameters.+?\*+\n\*+\n','');

if ~isempty(paramsummary)
    logtxt = [paramsummary{1} logtxt];
end

% Replace old file
fid = fopen(logfile,'wt');
fwrite(fid,logtxt);
fclose(fid);

end