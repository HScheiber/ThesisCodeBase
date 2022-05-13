function makedir(Directory)
if ispc
    system(['wsl mkdir ' windows2unix(Directory)]);
elseif isunix
    mkdir(Directory)
end
end