def Check_Structures(WorkDir, Salt, SystemName=None,
                        SaveTrajectory=True, SaveFeatures=False, 
                        SavePredictions=False, SavePredictionsImage=True, 
                        ML_TimeLength=5, ML_TimeStep=1, TimePerFrame=1, 
                        FileType='gro', Verbose=False, StartPoint = None,
                        EndPoint = None, Version = 2, SaveDir = None,
                        InMemory = False, Temporal_Cutoff = 0):
    
    """
    Created on Fri Nov 20 18:54:17 2020
    
    This function determines whether a LiX simulation has melted or frozen, 
    and returns a list of outputs.
    Outputs:
        
    Inputs (Required):
        WorkDir                 The directory containing the trajectory, 
                                and will also contain any output files.
        Salt                    The LiX salt.
        
    Inputs (Optional):
        SystemName              The name of the files involved in this calculation.
        T                       The actual system temperature.
        SaveTrajectory          Saves an .xyz trajectory of the predicted 
                                results when true.
        SaveFeatures            A switch to save the calculated features to a 
                                .pkl file.
        SavePredictions         A switch to saved the calculated structure 
                                predictions to a .pkl file.
        SavePredictionsImage    A switch to save a png image of the calculated 
                                structure predictions.
        ML_TimeLength           Controls the particular structure predictor 
                                model to use.
        ML_TimeStep             Controls the particular structure predictor 
                                model to use.
        TimePerFrame            Sets the step size to use between observations
                                in ps
        FileType                Sets the structure file type, either gro 
                                or g96
        Verbose                 Sets the level of logging verbosity.
        Startpoint              Initial point to start calculation in the 
                                trajectory. Default behaviour is to start at 
                                the first available point in the trajectory.
        Endpoint                Last time point to end calculation. Default 
                                behaviour is to end at the last available time
                                point in the trajectory.
        InMemory                When true, loads the trajectory directly into 
                                memory, rather than operate on the trajectory 
                                on disk.
        Temporal_Cutoff         [Frames] Eliminates short-term fluctuations of solids 
                                in the liquid phase whose timescale is shorter 
                                than this value.
    @author: Hayden
    """
    
    # % Import libraries and define functions
    import os
    threads = os.cpu_count()
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
    import tensorflow as tf
    import MDAnalysis as md
    import numpy as np
    import freud
    freud.parallel.set_num_threads(nthreads=threads)
    import time
    import re
    import warnings
    import logging
    
    if SavePredictionsImage:
        #from scipy.ndimage.filters import uniform_filter1d
        import seaborn as sns
        import matplotlib.pyplot as plt
    
    if SaveFeatures or SavePredictions:
        import pickle
    
    physical_devices = tf.config.list_physical_devices('GPU')
    if len(physical_devices) > 0:
        tf.config.experimental.set_memory_growth(physical_devices[0], True)
    
    if Verbose:
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        logging.basicConfig(level=logging.INFO,format='%(message)s')
    else:
        logging.basicConfig(level=logging.WARNING,format='%(message)s')
    
    logging.info("\rModules successfully imported")
    
    def to_freud_box(dims):
        a = dims[0]/10
        b = dims[1]/10
        c = dims[2]/10
        alpha = np.deg2rad(dims[3])
        beta = np.deg2rad(dims[4])
        gamma = np.deg2rad(dims[5])
    
        Omega = a*b*c*np.sqrt(1 - (np.cos(alpha)**2) - (np.cos(beta)**2) -
                              (np.cos(gamma)**2) +
                              2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        cx = c*np.cos(beta)
        cy = c*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        cz = Omega/(a*b*np.sin(gamma))
    
        Transform_Matrix = [[a,             b*np.cos(gamma),      cx],
                            [0,             b*np.sin(gamma),      cy],
                            [0,             0,                    cz]]
        return freud.box.Box.from_matrix(Transform_Matrix)
    
    def hist_laxis(data, n_bins, range_limits):
        # Setup bins and determine the bin location for each element for the bins
        R = range_limits
        N = data.shape[-1]
        bins = np.linspace(R[0],R[1],n_bins+1)
        data2D = data.reshape(-1,N)
        idx = np.searchsorted(bins, data2D,'right')-1
    
        # Some elements would be off limits, so get a mask for those
        bad_mask = (idx==-1) | (idx==n_bins)
    
        # We need to use bincount to get bin based counts. To have unique IDs for
        # each row and not get confused by the ones from other rows, we need to 
        # offset each row by a scale (using row length for this).
        scaled_idx = n_bins*np.arange(data2D.shape[0])[:,None] + idx
    
        # Set the bad ones to be last possible index+1 : n_bins*data2D.shape[0]
        limit = n_bins*data2D.shape[0]
        scaled_idx[bad_mask] = limit
    
        # Get the counts and reshape to multi-dim
        counts = np.bincount(scaled_idx.ravel(),minlength=limit+1)[:-1]
        counts.shape = data.shape[:-1] + (n_bins,)
        return counts
    
    def Separate_Metal_Halide(Salt):
        Metal = ''
        Halide = ''
        mylist = ['M', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Fe', 'Al', 'F', 'Cl', 'Br', 'I', 'At', 'X', 'O']
        for idx, s in enumerate(mylist):
            if s in Salt and idx < 8:
                Metal = s
            elif s in Salt:
                Halide = s
                return [Metal,Halide]
        
        return [Metal,Halide]
    
    
    # Calculation setup and loading the ML model
    
    # Check that the input directory is valid
    if not os.path.isdir(WorkDir):
        raise NotADirectoryError(WorkDir + ' is not a valid directory.')
    os.chdir(WorkDir)
    
    # Load some model parameters
    Trial_ID = str(int(ML_TimeLength)) + '-' + str(int(ML_TimeStep))
    if Trial_ID == '0-1':
        N_neighbour_list = [4,6,18] # List of nearest neighbour numbers to calculate order parameters
        L_list = [[2,4],[2,4,6],[2,4,6,8]] # List of lists: for each neighbour number, a list of l values for the order parameters
        timeless = True
    elif Trial_ID == '0-0':
        N_neighbour_list = [18] # List of nearest neighbour numbers to calculate order parameters
        L_list = [[2,4,6,8]] # List of lists: for each neighbour number, a list of l values for the order parameters
        timeless = True
    else:
        N_neighbour_list = [18] # List of nearest neighbour numbers to calculate order parameters
        L_list = [[2,4,6,8]] # List of lists: for each neighbour number, a list of l values for the order parameters
        timeless = False
    
    if Version == 1:
        ML_Name = 'V1NN-' + Trial_ID
    elif Version == 2:
        ML_Name = 'V2NN-' + Trial_ID
    Include_Wl = True # includes third order steinhart order parameters (Wl) as features when true
    Include_ID = True # includes atom identity (metal vs halide) as a feature when true
    
    # Load the ML model
    pwd = os.path.abspath(os.path.dirname(__file__))
    MLModelDir = os.path.join(pwd, "ML_Models/")
    
    if Version == 1:
        model_loc = os.path.join(MLModelDir,'LiX_Sturcture_Selector_CNN_Model_' + Trial_ID + '.tf')
    elif Version == 2:
        model_loc = os.path.join(MLModelDir,'MX_Sturcture_Classifier_Model_' + Trial_ID + '.tf')
    
    if not os.path.isdir(model_loc):
        raise FileNotFoundError('Unable to load Structure Selector model from: ' + model_loc + ' (model not found).')
    model = tf.keras.models.load_model(model_loc)
    
    # Separate the metal and halide 
    [Metal,Halide] = Separate_Metal_Halide(Salt)
    
    # If note give, attempt to find system name
    # Subdirectory after the salt name should be the system name
    if SystemName is None:
        p = re.compile(Salt + '/(.+?)/{0,1}')
        SystemName = p.findall(WorkDir)[0]
    
    if Version == 1:
        labels = ["Liquid","Rocksalt","Wurtzite","5-5","NiAs","Sphalerite","$\\beta$-BeO"]
        map_dict = {0: "Liquid",
                    1: "Rocksalt",
                    2: "Wurtzite",
                    3: "5-5",
                    4: "NiAs",
                    5: "Sphalerite",
                    6: "b-BeO"}
        map_label = {0: "L",
                    1: "R",
                    2: "W",
                    3: "F",
                    4: "N",
                    5: "S",
                    6: "B"}
    else:
        labels = ["Liquid","Rocksalt","Wurtzite","5-5","NiAs","Sphalerite","$\\beta$-BeO","AntiNiAs","CsCl"]
        map_dict = {0: "Liquid",
                    1: "Rocksalt",
                    2: "Wurtzite",
                    3: "5-5",
                    4: "NiAs",
                    5: "Sphalerite",
                    6: "b-BeO",
                    7: "AntiNiAs",
                    8: "CsCl"}
        map_label = {0: "L",
                    1: "R",
                    2: "W",
                    3: "F",
                    4: "N",
                    5: "S",
                    6: "B",
                    7: "A",
                    8: "C"}
    
    n_classes = len(map_dict)
    
    
    trajfile = os.path.join(WorkDir, SystemName + "." + 'trr')
    grofile = os.path.join(WorkDir, SystemName + "." + FileType)
    mdpfile = os.path.join(WorkDir, SystemName + "." + 'mdp')
    
    # Check that the required files exist
    if not os.path.isfile(trajfile):
        raise FileNotFoundError('Unable to find trajectory output file: ' + trajfile + ' (file not found).')
    if not os.path.isfile(grofile):
        # If no gro file found with default name, search for others
        import glob
        list_of_gro_files = glob.glob(os.path.join(WorkDir,'*.' + FileType))
        list_of_gro_files = [x for x in list_of_gro_files if not 'UnitCell' in x]
        list_of_gro_files = [x for x in list_of_gro_files if not 'OutConf' in x]
        if len(list_of_gro_files) == 0:
            raise FileNotFoundError('Unable to find initial system geometry file: ' + grofile + ' (file not found).')
        else:
            grofile = min(list_of_gro_files, key=os.path.getctime)
    if not os.path.isfile(mdpfile):
        raise FileNotFoundError('Unable to find molecular dynamics parameters file: ' + mdpfile + ' (file not found).')
    
    # Check if pbc is on
    textfile = open(mdpfile, 'r')
    filetext = textfile.read()
    textfile.close()
    if re.findall("pbc\s*=\s*xyz", filetext):
        pbc_on = True
    else:
        pbc_on = False
        
    # Load the trajectory and grab some basic info
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t = md.Universe(grofile, trajfile, in_memory=InMemory)
    
    traj_timestep = (t.trajectory[-1].time - t.trajectory[0].time)/(t.trajectory.n_frames-1) # ps, time per trajectory frame
    if StartPoint is None:
        min_time = t.trajectory[0].time  # ps
    else:
        min_time = StartPoint   # ps
    if EndPoint is None:
        max_time = t.trajectory[-1].time  # ps
    else:
        max_time = EndPoint
    min_step = int(min_time/traj_timestep)
    max_step = int(max_time/traj_timestep)
    
    min_traj_time = t.trajectory[0].time  # ps
    max_traj_time = t.trajectory[-1].time  # ps
    min_traj_step = int(min_traj_time/traj_timestep)
    max_traj_step = int(max_traj_time/traj_timestep)
    
    steps_per_frame = int(max(ML_TimeStep/traj_timestep,1))
    
    # output file name of processed 
    if SaveDir == None:
        outfile = ML_Name + '_' + Salt + '_' + SystemName + '.xyz'
    else:
        outfile = os.path.join(SaveDir, ML_Name + '_' + Salt + '_' + SystemName + '.xyz')
    
    # List of index points to examine
    steps_per_init_frame = int(TimePerFrame/traj_timestep)
    Traj_starts = list(range(min_step, max_step+1, steps_per_init_frame)) # Steps
    ML_TimeLength_in_steps = int(np.ceil(ML_TimeLength/traj_timestep)) # number of trajectory steps required to traverse the CNN time slice
    Half_ML_TimeLength_in_steps = int(np.floor(ML_TimeLength_in_steps/2))
    
    
    # What is the number of features used in this model (useful for preallocation)
    n_features = 0
    for Ls in L_list:
        n_features = n_features + len(Ls)*(1+int(Include_Wl))
    if Include_ID:
        n_features = n_features + 1
    
    t_slice_len = int(np.ceil(max(ML_TimeLength_in_steps,1)/steps_per_frame))
    num_traj_starts = len(Traj_starts)
    
    #% Calculate features
    # Loop through the number of time batches
    t_tot = time.time()
    Ql_result = np.empty([n_features,t.atoms.n_atoms])
    Ql_result_traj = np.empty([t_slice_len,t.atoms.n_atoms,n_features])
    Ql_result_all = np.empty([num_traj_starts,t_slice_len,t.atoms.n_atoms,n_features])
    
    logging.info('Generating features...')
    if Include_ID:
        namelist = t.atoms.names
        Metal_Index = np.array([x not in Metal for x in namelist]).astype(int) # 0 = metal, 1 = halide
    
    for traj_start_idx, Init_step in enumerate(Traj_starts):
        t_cur = time.time()
        
        # Select out a slice of the trajectory of the correct ML_TimeLength [steps]
        traj_slice = np.array(range(Init_step-Half_ML_TimeLength_in_steps, Init_step+Half_ML_TimeLength_in_steps+1, steps_per_frame))
        
        # Make sure there are not step indeces outside of the trajectory range by shifting back extremes
        if any(traj_slice < min_traj_step):
            ds = min_traj_step - np.min(traj_slice)
            traj_slice = traj_slice + ds
        elif any(traj_slice > max_traj_step):
            ds = max_traj_step - np.max(traj_slice)
            traj_slice = traj_slice + ds
        
        # Calculate Ql values at each time point in the time slice
        prev_time = np.nan;
        for t_idx, ts in enumerate(t.trajectory[traj_slice]):
            
            # Check if previous time slice is the same
            if ts.time == prev_time:
                Ql_result_traj[t_idx] = np.column_stack(Ql_result)
                continue
            
            # update previous time
            prev_time = ts.time
            
            if pbc_on:
                box_data = to_freud_box(ts.dimensions)
            else:
                # If pbc is off, place particles in ~infinite box
                box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
            
            # Select out the atoms
            point_data = ts.positions/10 # in nm
            system = [box_data,ts.positions/10]
            
            # Loop through the N_neighbour_list to built up a set of features
            fidx = 0
            if Include_ID:
                Ql_result[fidx] = Metal_Index # Append the atom idenities: metal vs non-metal
                fidx += 1
            for Neighbour_idx, N_neighbour  in enumerate(N_neighbour_list):
                
                # Construct the neighbour filter
                query_args = dict(mode='nearest', num_neighbors=N_neighbour, exclude_ii=True)
                
                # Build the neighbour list
                nlist = freud.locality.AABBQuery(box_data, point_data).query(point_data, query_args).toNeighborList()
                
                for L in L_list[Neighbour_idx]:
                    ql = freud.order.Steinhardt(L,wl=Include_Wl,wl_normalize=True)
                    ql_calc = ql.compute(system, neighbors=nlist)
                    
                    # Append order parameters to output
                    Ql_result[fidx] = ql_calc.ql
                    if Include_Wl:
                        Ql_result[fidx+1] = ql_calc.particle_order
                        fidx += 2
                    else:
                        fidx += 1
            Ql_result_traj[t_idx] = np.column_stack(Ql_result)
            
        Ql_result_all[traj_start_idx] = Ql_result_traj
        logging.info("\rTime Point: {:.2f} ps. Time Elapsed: {:.1f} s. ({:.2f}%, {:.2f} time points/s)".format(
            Init_step*traj_timestep,
            time.time() - t_tot,
            (traj_start_idx+1)*100/num_traj_starts,
            (1)/(time.time() - t_cur)))
    
    if np.shape(Ql_result_all)[1] == 1:
        Ql_result_all = np.squeeze(Ql_result_all,axis=1)
    
    
    # Optional: Save the calculated features
    if SaveFeatures:
        with open(ML_Name + '_' + Salt + '_' + SystemName + '_Features.pkl', 'wb') as f:
            pickle.dump([Ql_result_all], f)
    
    # Calculate the prediction for each atom at each time step
    logging.info('\nInferring Predictions...')
    t_tot = time.time()
    
    if timeless:
        Ql_result_all_t = Ql_result_all
    else:
        Ql_result_all_t = Ql_result_all.transpose(0, 2, 1, 3)
    
    # Concatenate all time slices into one
    Ql_result_all_cat = np.concatenate(Ql_result_all_t,axis=0)
    
    # Prediction probability
    predicted_prob = model.predict(Ql_result_all_cat, verbose=int(Verbose), batch_size = 10000)
    
    # Convert to predictions and convert back to time slices
    predicted_classes = np.array(np.split(np.argmax(predicted_prob, axis=1), num_traj_starts, axis=0))
    
    logging.info('\nFinished Inferring Predictions.')
    
    # Filter out transient predictions
    if Temporal_Cutoff > 0:
        t_tot = time.time()
        logging.info('\nRemoving short-timescale liquid fluctuations...')
        predicted_classes_grouped = predicted_classes.copy()
        for atm in range(0,np.shape(predicted_classes)[1]):
            c_atom_trj = (predicted_classes[:,atm].copy() == 0)
            
            i = 0
            for k, g in groupby(c_atom_trj):
                grp = len(list(g))
                if grp < Temporal_Cutoff:
                    predicted_classes_grouped[range(i, i+grp),atm] = 0
                i += grp
        
        predicted_classes = predicted_classes_grouped
        logging.info('\nFinished removing liquid fluctuations. Time elapsed: {:.1f} s'.format(time.time() - t_tot))
    
    
    if SaveTrajectory:
        t_tot = time.time()
        # Check if processed outfle already exists and delete
        if os.path.isfile(outfile):
            os.remove(outfile)
        
        # Convert predictions to their labels
        predicted_labels = np.vectorize(map_label.get)(predicted_classes)
        predicted_names = np.array([t.atoms.names + '_' + x for x in predicted_labels])
        
        W = md.coordinates.XYZ.XYZWriter(outfile, n_atoms=t.atoms.n_atoms, remark = '')
        for idx,ts in enumerate(t.trajectory[Traj_starts]):
            t_elap = time.time() - t_tot
            logging.info("\rWriting Time Point to file: {:.2f} ps. Time Elapsed: {:.1f} s. ({:.2f}%, {:.2f} time points/s)".format(
                ts.frame*traj_timestep,
                t_elap,
                (idx+1)*100/num_traj_starts,
                (idx+1)/t_elap))
            
            t.atoms.names = predicted_names[idx]
            
            dims = np.ndarray.flatten(t.coord.triclinic_dimensions)
            comment_txt = 'Lattice=' + '"' + str(dims[0]) + ' '\
                    + str(dims[1]) + ' ' + str(dims[2]) + ' ' + str(dims[3])\
                    + ' ' + str(dims[4]) + ' ' + str(dims[5]) + ' ' + str(dims[6])\
                    + ' ' + str(dims[7]) + ' ' + str(dims[8]) + '"'\
                    + ' Properties=species:S:1:pos:R:3 Time=' + str(ts.time)
            W.remark = comment_txt
            W.write(t.atoms)
        W.close()
    
    # Optional: save the predictions
    if SavePredictions:
        if SaveDir == None:
            with open(ML_Name + '_' + Salt + '_' + SystemName + '_Predictions.pkl', 'wb') as f:
                pickle.dump([predicted_prob,predicted_classes], f)
        else:
            with open(ML_Name + '_' + Salt + '_' + SystemName + '_Predictions.pkl', 'wb') as f:
                pickle.dump(os.path.join(SaveDir,[predicted_prob,predicted_classes]), f)
    
    # Finally, generate a series of line plots: partitions vs time
    if SavePredictionsImage:
        y_max = 100
        
        if Temporal_Cutoff > 0:
            init = Temporal_Cutoff-1
            finit = np.shape(predicted_classes_grouped)[0]-(Temporal_Cutoff-1)
        else:
            init = 0
            finit = np.shape(predicted_classes_grouped)[0]
        
        # Separate into metal and halide
        Metal_Ind = (Metal_Index == 0) # 0 = metal, 1 = halide
        Halide_Ind = (Metal_Index == 1) # 0 = metal, 1 = halide
        
        # Generate the histogram along one dimension
        d_metal = hist_laxis(predicted_classes[init:finit,Metal_Ind], n_classes, [0,n_classes])
        d_metal_norm = d_metal/sum(Metal_Ind)
        d_halide = hist_laxis(predicted_classes[init:finit,Halide_Ind], n_classes, [0,n_classes])
        d_halide_norm = d_halide/sum(Halide_Ind)
        
        #d = hist_laxis(predicted_classes, n_classes, [0,n_classes])
        #d_norm = d/np.shape(predicted_classes)[1]
        
        x = np.array([i * traj_timestep + min_traj_time for i in Traj_starts[init:finit,Halide_Ind]])
        #x = x - x[0]
        
        # colors
        cols = sns.color_palette("muted",n_classes)
        
        # Set up the matplotlib figure
        fig, ax = plt.subplots()
        for idx,y_dat_metal in enumerate(d_metal_norm.T):
            #yhat = uniform_filter1d(y_dat,size=int(np.ceil(len(y_dat)/30)))
            #plt.plot(x,yhat*100, label=labels[idx], color=cols[idx], linewidth=2, linestyle='dashed')
            y_dat_halide = d_halide_norm.T[idx,:]
            plt.plot(x,y_dat_metal*100, label=labels[idx], color=cols[idx], alpha=0.8, linewidth=2, linestyle='solid')
            plt.plot(x,y_dat_halide*100, label=None, color=cols[idx], alpha=0.8, linewidth=2, linestyle='dashed')
            #y_dat = d_norm[idx,:]
            #plt.plot(x,y_dat_metal*100, label=labels[idx], color=cols[idx], alpha=0.8, linewidth=3, linestyle='solid')
        
        
        plt.title('[' + ML_Name + ']: ' + Salt + ' ' + SystemName)
        plt.grid()
        plt.xlabel('Trajectory Time [ps]')
        plt.ylabel('Fraction in Structure [%]')
        ax.set_yticks(np.arange(0, 110, step=10))
        plt.legend()
        if y_max < 100:
            plt.ylim([0,y_max])
            
        if SaveDir == None:
            fig.savefig(ML_Name + '_' + Salt + '_' + SystemName + '.png', dpi=300, format='png')
        else:
            fig.savefig(os.path.join(SaveDir,ML_Name + '_' + Salt + '_' + SystemName + '.png'), dpi=300, format='png')
        
    return