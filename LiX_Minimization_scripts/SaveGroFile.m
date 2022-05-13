function txt = SaveGroFile(Filename,Geometry,saveflag)
    [~,~,ext] = fileparts(Filename);
    
    if strcmp(ext,'.gro')
        txt = [Geometry.Salt newline ' ' num2str(Geometry.N) newline];
        
        for jj=1:Geometry.N
            txt = [txt pad(num2str(Geometry.res_number(jj)),5,'left') ...
                pad(Geometry.res_name{jj},5,'right') ...
                pad(Geometry.atom_name{jj},5,'left') ...
                pad(num2str(mod(Geometry.atom_number(jj),1e5)),5,'left') ...
                pad(num2str(Geometry.xyz(jj,1),'%8.3f'),8,'left')  ...
                pad(num2str(Geometry.xyz(jj,2),'%8.3f'),8,'left') ...
                pad(num2str(Geometry.xyz(jj,3),'%8.3f'),8,'left') newline];
        end
        
        txt = [txt '   ' num2str([Geometry.boxcoords{:}],'%10.5f') newline];
    elseif strcmp(ext,'.g96')
        
        txt = ['TITLE' newline '    ' Geometry.Salt newline 'END' newline 'POSITION' newline];
        for jj=1:Geometry.N
            txt = [txt pad(pad(num2str(Geometry.res_number(jj)),5,'left'),6,'right') ...
                pad(Geometry.res_name{jj},6,'right') ...
                pad(Geometry.atom_name{jj},6,'right') ...
                pad(num2str(Geometry.atom_number(jj)),6,'right') ...
                pad(num2str(Geometry.xyz(jj,1),'%15.9f'),15,'left')  ...
                pad(num2str(Geometry.xyz(jj,2),'%15.9f'),15,'left') ...
                pad(num2str(Geometry.xyz(jj,3),'%15.9f'),15,'left') newline];
            
            
        end
        txt = [txt 'END' newline 'BOX' newline];
        % Add in box
        txt = [txt '    ' num2str([Geometry.boxcoords{:}],'%15.9f') newline 'END' newline];
    end
    
    
    % Save to file
    if saveflag
        fid = fopen(Filename,'wt');
        fwrite(fid,regexprep(txt,'\r',''));
        fclose(fid);
    end
end