WorkDir = 'C:\Users\Hayden\Documents\Patey_Lab\Melting_Point_Studies\NaCl\Test_R_TF_NPT\T_750.0000';
Salt = 'NaCl';



Test = py.LiXStructureDetector.Check_Structures(WorkDir, Salt,...
            pyargs('SystemName','Test_R_TF_NPT',...
            'FileType','gro',...
            'Verbose',true,...
            'SavePredictionsImage',true))

        KMP_DUPLICATE_LIB_OK
        
        from scipy.ndimage.filters import uniform_filter1d
        import seaborn as sns
        import matplotlib.pyplot as plt
        
        py.importlib.import_module('matplotlib.pyplot')
        
        py.seaborn.color_palette("muted",int64(7))
        
        py.matplotlib.pyplot.subplots
        fig, ax = plt.subplots()