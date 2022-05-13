Directory = '/home/scheiber/project/Molten_Salts_MD';
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl' 'CsCl'};

for idx = 1:length(Salts)
    idir = fullfile(Directory,Salts{idx});  
	files = dir(idir);
	isdir = logical(cat(1,files.isdir));
	dirs = files(isdir); % select only directory entries from the current listing
    dirs(cellfun(@(x) strncmp(x,'.',1),{dirs.name})) = [];
    
    for jdx = 1:length(dirs)
        jdir = fullfile(dirs(jdx).folder,dirs(jdx).name);
        JobName = dirs(jdx).name;
        if isfile(fullfile(jdir,'POSTPROCESS_COMPLETE'))
            disp(['Job Complete: ' Salts{idx} ' ' dirs(jdx).name])
            cd(jdir)
            system('find . -type f -iname "*-0*.subm" -delete','-echo');
            system('find . -type f -iname "*.std*" -delete','-echo');
            system('find . -type f -iname "#*#" -delete','-echo');
%             system('ls -tp | grep ''cpt'' | tail -n +2 | xargs rm; for i in *cpt; do mv "$i" "$(echo "$i" | sed -E ''s/-[0-9][0-9][0-9]\.cpt/.cpt/'')"; done');
            
        end
        
        
    end
    
end