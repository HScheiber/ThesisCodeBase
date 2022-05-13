function gitload
disp('Syncing ~/ThesisCodeBase with github...')
[~,Servertxt] = system('hostname -s | cut -c 1-3');
Server = strtrim(Servertxt);
workdir = pwd;
if ispc
    cd('C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase');
    !git pull
elseif strcmp(Server,'pat')
    cd('/home/user/ThesisCodeBase');
    !git pull
else
    cd('~/ThesisCodeBase');
    [~,cmdout] = system('git pull');
    disp(cmdout)
end
cd(workdir)
disp('Sync Complete.')
end