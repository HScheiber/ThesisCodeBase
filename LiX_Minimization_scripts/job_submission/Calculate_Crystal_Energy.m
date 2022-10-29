function E = Calculate_Crystal_Energy(Param,Settings)

if ~isempty(Param)
    for idx = 1:length(Settings.DOF)
        Settings.Geometry.(Settings.DOF{idx}) = Param(idx);
    end
    
    if Settings.MinMDP.Maintain_Symmetry
        Settings.Geometry = SymmetryAdapt(Settings.Geometry,Settings.Structure); %#ok<*UNRCH>
    end
end

Settings.WorkDir = tempname(Settings.WorkDir);

OptimizationLoop(Settings);

E = GrabEnergy(Settings.WorkDir,Settings.FileBase);

% Clean up
rmdir(Settings.WorkDir,'s')

end


