function E = Calculate_Crystal_Energy(Param,Settings)

if ~isempty(Param)
    for idx = 1:length(Settings.DOF)
        Settings.Geometry.(Settings.DOF{idx}) = Param(idx);
    end
    
    if Settings.MinMDP.Maintain_Symmetry
        Settings.Geometry = SymmetryAdapt(Settings.Geometry,Settings.Structure); %#ok<*UNRCH>
    end
end

if strcmp(Settings.Theory,'TF') && Settings.TF_Paramset == 3
    [TF_U_MXnew, TF_U_MMnew, TF_U_XXnew,d4fail] = TF_Potential_Generator(Settings);

    % If D4 fails, do not update the tables
    if ~d4fail
        TF_U_MX = TF_U_MXnew;
        TF_U_MM = TF_U_MMnew;
        TF_U_XX = TF_U_XXnew;
    end

    % Save tables into current directory
    fidMX = fopen(Settings.TableFile_MX,'wt');
    fwrite(fidMX,regexprep(TF_U_MX,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidMX);

    fidMM = fopen(Settings.TableFile_MM,'wt');
    fwrite(fidMM,regexprep(TF_U_MM,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidMM);

    fidXX = fopen(Settings.TableFile_XX,'wt');
    fwrite(fidXX,regexprep(TF_U_XX,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidXX);
end

Settings.WorkDir = tempname(Settings.WorkDir);

OptimizationLoop(Settings);

E = GrabEnergy(Settings.WorkDir,Settings.FileBase);

% Clean up
rmdir(Settings.WorkDir,'s')

end


