function coords_out = Get_Fractional_Coords(logtext)
    % Get the fractional coordinates of asymmetric unit atoms
    num_UC = regexp(logtext,...
        '(ATOMS IN THE UNIT CELL: +)(-|\+|\.|[0-9]|E)+','tokens');
    UC_lines = str2double(num_UC{1}{2});

    coords_UC = regexp(logtext,...
        ['(ATOMS IN THE UNIT CELL.+?\n)(.+?\n){' num2str(UC_lines + 2) '}'],'tokens');

    % extract coordinates
    coords_UC = textscan(coords_UC{1}{2},'%*f %s %*f %s %f %f %f',...
        'MultipleDelimsAsOne',true,'Delimiter',' ','HeaderLines',2);
    coords_out = cell(UC_lines,4);
    indxx = 0;
    for indx = 1:UC_lines
        if strcmp(coords_UC{1}{indx},'T') % if T, then part of asymmetric unit
            indxx = indxx + 1;
            coords_out{indxx,1} = coords_UC{2}{indx};
            coords_out{indxx,2} = coords_UC{3}(indx);
            coords_out{indxx,3} = coords_UC{4}(indx);
            coords_out{indxx,4} = coords_UC{5}(indx);
        end
    end
    % Remove empty rows
    coords_out( all(cellfun(@isempty,coords_out),2),:) = [];
end