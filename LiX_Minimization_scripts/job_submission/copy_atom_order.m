function output = copy_atom_order(filename)
[~,~,ext] = fileparts(filename);

fid = fopen(filename,'rt');

if strcmp(ext,'.gro')
    Input_data = textscan(fid,'%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n',...
        'Delimiter','','whitespace','','headerlines',2);
    fclose(fid);

    N = size(Input_data{1,1},1)-1;
    Outcell = cell(1,N);
    Previous = strtrim(Input_data{3}{1});
    cnt = 0;
    idx = 1;
    for i=1:N
        Current = strtrim(Input_data{3}{i});
        if ~strcmp(Current,Previous)
            if i == N % last one but different from previous
                Outcell{idx} = [Previous ' ' num2str(cnt)];
                idx = idx+1;
                cnt = 1;
                Outcell{idx} = [Current ' ' num2str(cnt)];
            else % different from previous but not last one
                Outcell{idx} = [Previous ' ' num2str(cnt)];
                idx = idx+1;
                cnt = 1;
            end
        elseif i == N % last one but same as previous
            cnt = cnt+1;
            Outcell{idx} = [Current ' ' num2str(cnt)];
        else % same as previous but not last one
            cnt = cnt+1;
        end        
        Previous = Current;
    end
elseif strcmp(ext,'.g96')
    curline = fgetl(fid);
    while ~strcmp(curline,'POSITION')
        curline = fgetl(fid);
    end
    Input_data = textscan(fid,'%6c%6c%6c%6c%15.9f%15.9f%15.9f\n',...
        'Delimiter','','whitespace','');
    fclose(fid);
    
    N = size(Input_data{3},1)-1;
    Outcell = cell(1,N);
    Previous = strtrim(Input_data{3}(1,:));
    cnt = 0;
    idx = 1;
    for i=1:N
        Current = strtrim(Input_data{3}(i,:));
        if ~strcmp(Current,Previous)
            if i == N % last one but different from previous
                Outcell{idx} = [Previous ' ' num2str(cnt)];
                idx = idx+1;
                cnt = 1;
                Outcell{idx} = [Current ' ' num2str(cnt)];
            else % different from previous but not last one
                Outcell{idx} = [Previous ' ' num2str(cnt)];
                cnt = 1;
                idx = idx+1;
            end
        elseif i == N % last one but same as previous
            cnt = cnt+1;
            Outcell{idx} = [Current ' ' num2str(cnt)];
        else % same as previous but not last one
            cnt = cnt+1;
        end        
        Previous = Current;
    end
end

Outcell = Outcell( cellfun(@(x) ~isempty(x),Outcell) );
W = [Outcell',[repmat({newline},numel(Outcell)-1,1);{[]}]]';
output = [W{:}];
end