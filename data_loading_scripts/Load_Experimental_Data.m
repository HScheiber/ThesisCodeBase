function Experiment = Load_Experimental_Data
    home = find_home;
    Data_Directory = [home filesep 'data'];
    Exp = load(fullfile(Data_Directory,'Alkali_Halides_Exp.mat'),'Experiment');
    Experiment = Exp.Experiment;
end