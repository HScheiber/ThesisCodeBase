function subFolders = find_directories(outer_files)

% Find directories in current folder, except . and ..
outer_files(ismember( {outer_files.name}, {'.', '..'})) = [];

% Get a logical vector that tells which is a directory.
dirFlags = [outer_files.isdir];

% Extract only those that are directories.
subFolders = { outer_files(dirFlags).name };

end