% Script to run any post-processing to MD simulations
function Output = Alkali_Halide_Classifier_Tester(Settings)

PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(Settings.WorkDir, Settings.Salt,...
            pyargs('SystemName',Settings.JobName,...
            'RefStructure',Settings.Structure,...
            'FileType',Settings.CoordType,...
            'SavePredictionsImage',Settings.SavePredictionsImage,...
            'ML_TimeLength',Settings.ML_TimeLength,...
            'ML_TimeStep',Settings.ML_TimeStep,...
            'TimePerFrame',Settings.TimePerFrame,...
            'SaveTrajectory',Settings.SaveTrajectory,...
            'SaveFeatures', Settings.SaveFeatures,...
            'SavePredictions',Settings.SavePredictions,...
            'Qlm_Average',Settings.Qlm_Average,...
            'Voronoi',Settings.Voronoi));
    
    %Output.Froze = logical(PyOut{1});
    %Output.Melted = logical(PyOut{2});
    Output = double(PyOut{4}); % Ref fraction
    %Output.LiqFrac = double(PyOut{5});
    %Output.Froze_alt = logical(PyOut{6});
end