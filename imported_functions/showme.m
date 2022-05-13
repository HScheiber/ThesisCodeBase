function showme
[~,uname_out] = system('uname -n | cut -c 1-3');
uname = strtrim(uname_out);
if strcmp(uname,'sea')
    [~,cmdout] = system('showq -u $USER');
else
    [~,cmdout] = system('squeue -lu $USER');
end
disp(cmdout)
end