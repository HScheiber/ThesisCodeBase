function myCellSelectionCB_Old(hObject,eventdata,Data_Table,home)
    if isempty(eventdata.Indices)
        return
    end
    row = eventdata.Indices(1,1);
    column = eventdata.Indices(1,2);
    cellDat = strtrim(regexprep(hObject.Data{row,column}, '<.*?>',''));
    if strcmp(cellDat,'-') % Do nothing for empty cells
        return
    end
    
    % Get salt info
    salt = textscan(Data_Table{row,1}{1},'%s %*s');
    while iscell(salt)
        salt = salt{1};
    end
    [Met,Hal] = Separate_Metal_Halide(salt);
    
    if column == 11 % View structure no dispersion 
        a = sprintf('%8f',Data_Table{row,5});
        b = sprintf('%8f',Data_Table{row,6});
        c = sprintf('%8f',Data_Table{row,7});
        BL = sprintf('%8f',Data_Table{row,9});
        FC = textscan(Data_Table{row,10}{1},'(%f, %f, %f) (%f, %f, %f)');
        disptxt = '';
        disptxt2 = 'Without';
    elseif column == 19 % View structure dispersion
        a = sprintf('%8f',Data_Table{row,13});
        b = sprintf('%8f',Data_Table{row,14});
        c = sprintf('%8f',Data_Table{row,15});
        BL = sprintf('%8f',Data_Table{row,17});
        FC = textscan(Data_Table{row,18}{1},'(%f, %f, %f) (%f, %f, %f)');
        disptxt = '_D3';
        disptxt2 = 'With Empirical';
    else % Otherwise copy to clipboard if not empty
        clipboard('copy',cellDat)
        return
    end
        
    % Get fractional coordinates
    Metx = strrep(sprintf('%+9f',FC{1}),'+',' ');
    Mety = strrep(sprintf('%+9f',FC{2}),'+',' ');
    Metz = strrep(sprintf('%+9f',FC{3}),'+',' ');
    

    Halx = strrep(sprintf('%+9f',FC{4}),'+',' ');
    Haly = strrep(sprintf('%+9f',FC{5}),'+',' ');
    Halz = strrep(sprintf('%+9f',FC{6}),'+',' ');

    structure_cell = textscan(Data_Table{row,1}{1},'%s %s');
    structure = structure_cell{2}{1};
    theory = Data_Table{row,2}{1};
    if strcmp(theory,'JC (SPC/E)')
        theory = 'JC';
    elseif strcmp(theory,'JC (TIP4P/EW)')
        theory = 'JC4P';
    elseif strcmp(theory,'JC (TIP3P)')
        theory = 'JC3P';
    end

    filename = [salt '_' structure '_' theory '.vesta'];

    % Load template
    fid = fopen([home filesep 'templates' filesep structure '.template'],'rt');
    X = fread(fid);
    fclose(fid);
    template_text = char(X.');

    % Replace template strings
    template_text = strrep(template_text,'##INFO##',[salt ' ' theory ' ' disptxt2 ' Dispersion.']);
    template_text = strrep(template_text,'##APAR##',a);
    template_text = strrep(template_text,'##BPAR##',b);
    template_text = strrep(template_text,'##CPAR##',c);
    template_text = strrep(template_text,'##MET##',Met);
    template_text = strrep(template_text,'##METX##',Metx);
    template_text = strrep(template_text,'##METY##',Mety);
    template_text = strrep(template_text,'##METZ##',Metz);
    template_text = strrep(template_text,'##HAL##',Hal);
    template_text = strrep(template_text,'##HALX##',Halx);
    template_text = strrep(template_text,'##HALY##',Haly);
    template_text = strrep(template_text,'##HALZ##',Halz);
    template_text = strrep(template_text,'##BOND##',BL);

    % Save file
    tempfolder = fullfile(home,'analysis','AppData','temp');
    if ~exist(fullfile(home,'analysis','AppData','temp'),'dir')
        mkdir(tempfolder);
    end
    full_filename = fullfile(tempfolder,filename);

    fid2 = fopen(full_filename,'wt');
    fwrite(fid2,template_text);
    fclose(fid2);        

    winopen(full_filename)
end