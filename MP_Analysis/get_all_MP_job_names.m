function FullJobNames = get_all_MP_job_names(Data_Salt)

% loop through models
FullJobNames = cell(0,1);

Models = fieldnames(Data_Salt);

for idx = 1:length(Models)
    Model = Models{idx};
    
    Structures = fieldnames(Data_Salt.(Model));
    
    for jdx = 1:length(Structures)
        Structure = Structures{jdx};
        endtxt = ['_' StructureLabel(Structure) '_' Model '_NPT'];
        
        JobNames = fieldnames(Data_Salt.(Model).(Structure));
        
        FullJobNames = [FullJobNames; cellfun(@(x) [strrep(x,'MP_','') endtxt],JobNames,'UniformOutput',false)];
    end
end