% Short function to reset the geometry of a CRYSTAL17 input
function Crystal_Input_Template = Reset_Crystal_Input(NumCoord,MetalCAN,HalideCAN,OptDir,JobName,Cry,...
    Crystal_Input_text,Crystal_Geom_Template,Crystal_Input_Template,ResetMainGeom)

    % Generate lattice parameters
    a_param = num2str(Cry.a,'%15.14f');
    b_param = num2str(Cry.b,'%15.14f');
    c_param = num2str(Cry.c,'%15.14f');
    
    % Generate coordinates
    Coordinates = num2str(NumCoord);
    for i = 1:NumCoord/2
        Coordinates = [Coordinates newline num2str(MetalCAN) ' ' num2str(Cry.FC_Metal(i,:),' %15.14f')];
        Coordinates = [Coordinates newline num2str(HalideCAN) ' ' num2str(Cry.FC_Halide(i,:),' %15.14f')];
    end

    % Reset the main central geometry used to classify integrals
    if ResetMainGeom
        Crystal_Input_Template = strrep(Crystal_Input_Template,'##APAR##',a_param);
        Crystal_Input_Template = strrep(Crystal_Input_Template,'##BPAR##',b_param);
        Crystal_Input_Template = strrep(Crystal_Input_Template,'##CPAR##',c_param);
        Crystal_Input_Template = strrep(Crystal_Input_Template,'##GEOM##',Coordinates);
        % Add in GUESSP to use previous density matrix as initial guess if not already done
        if isempty(regexp(Crystal_Input_Template,'GUESSP','ONCE'))
            Crystal_Input_Template = regexprep(Crystal_Input_Template,'END$',['GUESSP' newline 'END']);
        end
        
        % Save file
        fid = fopen([OptDir filesep JobName '.d12'],'wt');
        fwrite(fid,Crystal_Input_Template);
        fclose(fid);
        
    % Reset the new geometry of interest
    else

        Geom_text = strrep(Crystal_Geom_Template,'##APAR##',a_param);
        Geom_text = strrep(Geom_text,'##BPAR##',b_param);
        Geom_text = strrep(Geom_text,'##CPAR##',c_param);
        Geom_text = strrep(Geom_text,'##GEOM##',Coordinates);
        
        
        % Add in GUESSP to use previous density matrix as initial guess if not already done
        if isempty(regexp(Crystal_Input_text,'GUESSP','ONCE'))
            Crystal_Input_text = regexprep(Crystal_Input_text,'END$',['GUESSP' newline 'END']);
        end
    
        % Add new geometry to old input file
        Crystal_Input_text = regexprep(Crystal_Input_text,'END$',...
            ['FIXINDEX' newline 'END' newline 'GEOM' newline Geom_text]);
        % Save file
        fid = fopen([OptDir filesep JobName '.d12'],'wt');
        fwrite(fid,Crystal_Input_text);
        fclose(fid);
    end
end