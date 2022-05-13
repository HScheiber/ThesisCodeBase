function myCellSelectionCB(hObject,eventdata,Data_Table,home)
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
    salt = textscan(Data_Table{row,1}{1},'%s');
    while iscell(salt)
        salt = salt{1};
    end
    [Met,Hal] = Separate_Metal_Halide(salt);
    
    if column == 15 % View structure
        Optimization = Data_Table{row,4}{1};
        a = sprintf('%8f',Data_Table{row,6});
        b = sprintf('%8f',Data_Table{row,7});
        c = sprintf('%8f',Data_Table{row,8});
        alpha = sprintf('%8f',Data_Table{row,9});
        beta = sprintf('%8f',Data_Table{row,10});
        gamma = sprintf('%8f',Data_Table{row,11});
        BL = sprintf('%8f',Data_Table{row,13});
        FC = textscan(Data_Table{row,14}{1},'(%f, %f, %f)');
    else % Otherwise copy to clipboard if not empty
        clipboard('copy',cellDat)
        return
    end
    
    % If no structural data exists
    if isnan(Data_Table{row,6})
        return
    end
    
    N = length(FC{1});
    
    % Construct input text for P1
    if strcmp(Optimization,'Full (P1)')
        Met_Struc = '';
        Hal_Struc = '';
        Met_Theri = '';
        Hal_Theri = '';
        Met_Sitet = '';
        Hal_Sitet = '';
        for i=1:N/2
            Met_Struc = [Met_Struc '  ' num2str(i) ' ' pad(Met,2,'left') '         ' ...
                pad(Met,2,'left') '  1.0000   ' num2str(FC{1}(i),'%1.6f') ...
                '   ' num2str(FC{2}(i),'%1.6f') '   ' num2str(FC{3}(i),'%1.6f') '    1a       1' newline ...
                '                            0.000000   0.000000   0.000000  1.00' newline];

            Met_Theri = [Met_Theri '   ' num2str(i) '         ' pad(Met,2,'left') '  1.000000' newline];
            Met_Sitet = [Met_Sitet '   ' num2str(i) '         ' pad(Met,2,'left') '  1.0000 255   0   0 134 224 116 204  0' newline];

            Hal_Struc = [Hal_Struc '  ' num2str(N/2+i) ' ' pad(Hal,2,'left') '         ' ...
                pad(Hal,2,'left') '  1.0000   ' num2str(FC{1}(N/2+i),'%1.6f') ...
                '   ' num2str(FC{2}(N/2+i),'%1.6f') '   ' num2str(FC{3}(N/2+i),'%1.6f') '    1a       1' newline ...
                '                            0.000000   0.000000   0.000000 -1.00' newline];
            Hal_Theri = [Hal_Theri '   ' num2str(N/2+i) '         ' pad(Hal,2,'left') '  1.000000' newline];
            Hal_Sitet = [Hal_Sitet '   ' num2str(N/2+i) '         ' pad(Hal,2,'left') '  1.0000   0 255   0 176 185 230 204  0' newline];
        end
        Tot_Struc = [Met_Struc Hal_Struc];
        Tot_Theri = [Met_Theri Hal_Theri];
        Tot_Sitet = [Met_Sitet Hal_Sitet];
    else
        % Get fractional coordinates
        Metx = strrep(sprintf('%+9f',FC{1}(1:N/2)),'+',' ');
        Mety = strrep(sprintf('%+9f',FC{2}(1:N/2)),'+',' ');
        Metz = strrep(sprintf('%+9f',FC{3}(1:N/2)),'+',' ');

        Halx = strrep(sprintf('%+9f',FC{1}(N/2+1:N)),'+',' ');
        Haly = strrep(sprintf('%+9f',FC{2}(N/2+1:N)),'+',' ');
        Halz = strrep(sprintf('%+9f',FC{3}(N/2+1:N)),'+',' ');
    end

    structure_cell = textscan(Data_Table{row,2}{1},'%s');
    structure = structure_cell{1}{1};
    theory = Data_Table{row,3}{1};
    if strcmp(theory,'JC (SPC/E)')
        theory = 'JC';
    elseif strcmp(theory,'JC (TIP4P/EW)')
        theory = 'JC4P';
    elseif strcmp(theory,'JC (TIP3P)')
        theory = 'JC3P';
    end

    filename = [salt '_' structure '_' theory '_' strrep(Optimization,' ','-') '.vesta'];

    % Load template
    if strcmp(Optimization,'Full (P1)')
        fid = fopen([home filesep 'templates' filesep 'P1.template'],'rt');
        X = fread(fid);
        fclose(fid);
        template_text = char(X.');
        
        % Replace template strings
        template_text = strrep(template_text,'##INFO##',[salt ' ' structure ' ' theory '.']);
        template_text = strrep(template_text,'##APAR##',a);
        template_text = strrep(template_text,'##BPAR##',b);
        template_text = strrep(template_text,'##CPAR##',c);
        template_text = strrep(template_text,'##ALPHA##',alpha);
        template_text = strrep(template_text,'##BETA##',beta);
        template_text = strrep(template_text,'##GAMMA##',gamma);
        template_text = strrep(template_text,'##MET##',Met);
        template_text = strrep(template_text,'##HAL##',Hal);
        template_text = strrep(template_text,'##STRUC##',Tot_Struc);
        template_text = strrep(template_text,'##THERI##',Tot_Theri);
        template_text = strrep(template_text,'##SITET##',Tot_Sitet);
    else
        fid = fopen([home filesep 'templates' filesep structure '.template'],'rt');
        X = fread(fid);
        fclose(fid);
        template_text = char(X.');
        
        % Replace template strings
        template_text = strrep(template_text,'##INFO##',[salt ' ' structure ' ' theory '.']);
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
    end

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