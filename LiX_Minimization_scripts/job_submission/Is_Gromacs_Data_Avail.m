function skip = Is_Gromacs_Data_Avail(Salt,Structure,Model,TF_Param,C6_Damp,CR_Damp,S,...
    Optimize_Position,Maintain_Symmetry,Workdir,...
    GAdjust_MX,GAdjust_MM,GAdjust_XX)

    G_Adj = '';
    if sum(GAdjust_MX(:,1)) ~= 0
        G_Adj = [G_Adj '_GMX'];
    end

    if sum(GAdjust_MM(:,1)) ~= 0
        G_Adj = [G_Adj '_GMM'];
    end
    
    if sum(GAdjust_XX(:,1)) ~= 0
        G_Adj = [G_Adj '_GXX'];
    end
    
    skip = false;
    
    if Optimize_Position && Maintain_Symmetry
        Opt = 'FULLOPT';
    elseif Optimize_Position && ~Maintain_Symmetry
        Opt = 'FULLOPT_SG1';
    elseif ~Optimize_Position && ~Maintain_Symmetry
        Opt = 'CELLOPT_SG1';
    elseif ~Optimize_Position && Maintain_Symmetry
        Opt = 'CELLOPT';
    end

    %% Generate name for model with current scaling parameters
    if TF_Param ~= 0 && strcmp(Model,'TF')
        Model = [Model num2str(TF_Param)];
    end
    
    if C6_Damp.MX == C6_Damp.MM && C6_Damp.MX == C6_Damp.XX
        if C6_Damp.MX ~= 0
            Model = [Model 'd' num2str(C6_Damp.MX)];
        end
    else
        for current_interaction = fieldnames(C6_Damp)'
            tag = ['d' current_interaction{1}];

            if C6_Damp.(current_interaction{1}) ~= 0
                Model = [Model tag num2str(C6_Damp.(current_interaction{1}))];
            end
        end
    end
    
    for current_scale = fieldnames(S)'
        for current_interaction = fieldnames(S.(current_scale{1}))'
            if strcmp(current_interaction{1},'All')
                tag = ['_' current_scale{1}];
            else
                tag = ['_' current_interaction{1} current_scale{1}];
            end

            if ~ismembertol(1.0,S.(current_scale{1}).(current_interaction{1}),1e-5)
                mtxt = strrep(num2str(S.(current_scale{1}).(current_interaction{1}),'%10.5f'),'-','N');
                Model = [Model tag mtxt];
            end
        end
    end

    for current_interaction = fieldnames(CR_Damp)'
        if strcmp(current_interaction{1},'All')
            tag = '_dCR';
        else
            tag = ['_dCR' current_interaction{1}];
        end

        if CR_Damp.(current_interaction{1}).r_d >= 0 && CR_Damp.(current_interaction{1}).b >= 0
            mtxt = ['rd' num2str(CR_Damp.(current_interaction{1}).r_d,'%10.3f') 'b' num2str(CR_Damp.(current_interaction{1}).b,'%10.3f')];
            Model = [Model tag mtxt];
        end
    end
    
    Model = [Model G_Adj];
    
    Struc_Label = Label_replace(Structure);
    
    Search_Dir = fullfile(Workdir,Salt,Structure,Model,Opt);
    
    logfilename = [Salt Struc_Label '_' Model '_' Opt '.optlog'];
    
    fullfilename = fullfile(Search_Dir,logfilename);
    
    % Check if file exists
    if ~isfile(fullfilename)
        return
    end
    
    % If file exists, open file and extract text
    logtext = fileread(fullfilename);
    
    % Look for "Stopping Here." in log file. One case where restart needed
    if ~isempty(regexp(logtext,'Stopping here\.','ONCE'))
        return
    end
    
    % Look for convergence not found due to unphysical energy
    if ~isempty(regexp(logtext,'Unphysical','ONCE'))
        return
    end
    
    % Look for convergence found text
    if ~isempty(regexp(logtext,'Final Optimized','ONCE'))
        skip = true;
        return
    end
    
    
end