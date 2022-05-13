% Function wrapper for the system() command, allowing use of default system path
function [output,errorcode] = system_def(Input)
    % Save library paths
    MatlabPath = getenv('LD_LIBRARY_PATH');
    % Make Matlab use system libraries
    setenv('LD_LIBRARY_PATH',getenv('PATH'));
    [output,errorcode] = system(Input);
    % Reassign old library paths
    setenv('LD_LIBRARY_PATH',MatlabPath);
end