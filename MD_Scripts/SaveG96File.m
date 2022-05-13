function SaveG96File(Filename,Crystal)
    
    txt = ['TITLE' newline '    ' Crystal.Salt newline 'END' newline 'POSITION' newline];
    for jj=1:Crystal.N
        txt = [txt pad(pad(num2str(Crystal.res_number(jj)),5,'left'),6,'right') ...
            pad(Crystal.res_name{jj},6,'right') ...
            pad(Crystal.atom_name{jj},6,'right') ...
            pad(num2str(Crystal.atom_number(jj)),6,'right') ...
            pad(num2str(Crystal.xyz(jj,1),'%15.9f'),15,'left')  ...
            pad(num2str(Crystal.xyz(jj,2),'%15.9f'),15,'left') ...
            pad(num2str(Crystal.xyz(jj,3),'%15.9f'),15,'left') newline];
    end
    txt = [txt 'END' newline 'BOX' newline];

    % Add in box
    txt = [txt '    ' num2str([Crystal.boxcoords{:}],'%15.9f') newline 'END' newline];
    
    % Save to file
    fid = fopen(Filename,'wt');
    fwrite(fid,regexprep(txt,'\r',''));
    fclose(fid);
end