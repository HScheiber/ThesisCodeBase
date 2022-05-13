% Function to convert windows paths into unix paths
function output = windows2unix(input)

if isunix
    output = input;
    return
end
output = strrep(regexprep(input,'([A-Z]):','/mnt/${lower($1)}'),'\','/');


end