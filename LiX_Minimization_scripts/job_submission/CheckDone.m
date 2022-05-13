LatPars = 2:0.01:10;


N = length(LatPars);
LatPars_Formatted = cell(1,N);
NotDoneList = [];

for i = 1:N
    LatPars_Formatted{i} = ['LPAR_' pad(num2str(LatPars(i),'%5.3f'),6,'left','0')];
    
    if ~exist(LatPars_Formatted{i},'dir')
        NotDoneList(end+1) = LatPars(i);
    else
        if isempty(dir([LatPars_Formatted{i} filesep '*energies.xvg']))
            NotDoneList(end+1) = LatPars(i);
        end
    end
end

if ~isempty(NotDoneList)
    disp('Found Unfinished Data Points.')
    save('Unfinished.mat','NotDoneList')
else
    disp('Check Complete: All Data Points Located.')
end