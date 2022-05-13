% Script to rename models

outdir = '/home/scheiber/project/Model_Building';
Orig_Models = 'CT';
Rename_Models = 'CM';

dirs = dir(outdir);
for i=length(dirs):-1:1
	if strncmp( dirs(i).name,'.',1) || ...
            strcmpi( dirs(i).name,'completed')
        dirs(i) = [];
	end
end

for i=1:length(dirs)
    subdirs = dir(fullfile(outdir,dirs(i).name));
    
    for j=1:length(subdirs)
        
        % Find the files we want to change the names of
        if regexp(subdirs(j).name,Orig_Models)
            oldfile = fullfile(outdir,dirs(i).name,subdirs(j).name);
            
            % Create new file name that isn't already taken
            N = 0;
            while true
                N = N + 1;
                new_model_num = [Rename_Models num2str(N)];
                new_model_name = regexprep(subdirs(j).name,[Orig_Models '[0-9]{1,2}'],new_model_num,'once');
                
                newfile = fullfile(outdir,dirs(i).name,new_model_name);
                
                if ~exist(newfile, 'dir')
                    disp(oldfile)
                    disp(newfile)
                    break
                end
            end
            
            % Move the folder
            movefile(oldfile,newfile)
            
            % Find and rename all the files in the folder
            newdir = dir([newfile filesep '*' Orig_Models '*']);
            
            for k=1:length(newdir)
                oldobj = fullfile(newdir(k).folder,newdir(k).name);
                
                new_file_name = regexprep(newdir(k).name,[Orig_Models '[0-9]{1,2}'],new_model_num,'once');
                newobj = fullfile(newdir(k).folder,new_file_name);
                
                % Move the file
                movefile(oldobj,newobj)
                disp(oldobj)
                disp(newobj)
                
                % Rename the trial ID in the model input file
                if contains(newobj,'.inp')
                    dat = load(newobj,'-mat');
                    Model = dat.Model;
                    Model.Trial_ID = new_model_num;
                    save(newobj,'Model','-mat')
                end
                
                % Rename the trial ID in the log file
                if contains(newobj,'.log')
                    dat = fileread(newobj);
                    log = regexprep(dat,'Trial_ID:\s+''.{1,4}''',['Trial_ID: ''' new_model_num '''']);
                    fid = fopen(newobj,'wt');
                    fprintf(fid, log);
                    fclose(fid);
                end
            end
        end
    end
end
