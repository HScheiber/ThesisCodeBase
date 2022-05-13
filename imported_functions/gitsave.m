function gitsave
disp('Saving ~/ThesisCodeBase to github...')
[~,Servertxt] = system('hostname -s | cut -c 1-3');
Server = strtrim(Servertxt);
workdir = pwd;
if ispc
    cd('C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase');
    [~,cmdout1] = system('wsl source ~/.bashrc; git add -A');
    disp(cmdout1)
    [~,cmdout2] = system('wsl source ~/.bashrc; git commit -m "Incrimental save point."');
    disp(cmdout2)
    [~,cmdout3] = system('wsl source ~/.bashrc; git pull');
    disp(cmdout3)
    [~,cmdout4] = system('wsl source ~/.bashrc; git push');
    disp(cmdout4)
elseif strcmp(Server,'pat')
    cd('/home/user/ThesisCodeBase');
    !git add -A
    !git commit -m "Incrimental save point."
    !git pull
    !git push
else
    cd('~/ThesisCodeBase');
    [~,cmdout1] = system('git add -A');
    disp(cmdout1)
    [~,cmdout2] = system('git commit -m "Incrimental save point."');
    disp(cmdout2)
    [~,cmdout3] = system('git pull');
    disp(cmdout3)
    [~,cmdout4] = system('git push');
    disp(cmdout4)
end
cd(workdir)
disp('Save Complete.')
end