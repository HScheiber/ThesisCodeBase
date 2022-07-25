def Calculate_Liquid_Fraction(WorkDir, Salt, SystemName=None, T=None,
                              T_Ref=None, RefStructure='Liquid', CheckFullTrajectory=False, 
                              SaveTrajectory=False, SaveFeatures=False, 
                              SavePredictions=False, SavePredictionsImage=False,
                              InitialRefFrac=None, RefChangeThreshold=0.2, 
                              SlopeThreshold=0.0125, SlopeCheckBegin=0.1,
                              ML_TimeLength=20, ML_TimeStep=5, TimePerFrame=5, 
                              FileType='gro', Verbose=False, Version = 2,
                              Temporal_Cutoff = 0,Voronoi = False, Qlm_Average = True,
                              Prob_Interfacial = None,Spatial_Reassignment = False,
                              Spatial_Interfacial = None,
                              SaveTrajectoryAux = 0):
    
    """
    Created on Fri Nov 20 18:54:17 2020
    
    Calculate_Liquid_Fraction(WorkDir, Salt, SystemName=None, T=None,
                                  T_Ref=None, RefStructure='Liquid', CheckFullTrajectory=False, 
                                  SaveTrajectory=False, SaveFeatures=False, 
                                  SavePredictions=False, SavePredictionsImage=False,
                                  InitialRefFrac=None, RefChangeThreshold=0.2, 
                                  SlopeThreshold=0.0125, SlopeCheckBegin=0.1,
                                  ML_TimeLength=20, ML_TimeStep=5, TimePerFrame=5, 
                                  FileType='gro', Verbose=False, Version = 2,
                                  Temporal_Cutoff = 0,Voronoi = False, Qlm_Average = False,
                                  Prob_Interfacial = None,Spatial_Reassignment = False,
                                  Spatial_Interfacial = None)
    This function determines whether a LiX simulation has melted or frozen, 
    and returns a list of outputs.
    Outputs:
        system_froze            is true if the algorithm detects that the 
                                system freezes.
        system_melted           is true if the algorithm detects that the 
                                system melts.
        time_to_phase_change    is number indicating the time (in ps) at which 
                                the phase change occurs.
        
    Inputs (Required):
        WorkDir                 The directory containing the trajectory, 
                                and will also contain any output files.
        MLModelDir              The directory containing the ML structure 
                                predictor to be used.
        Salt                    The LiX salt.
        
    Inputs (Optional):
        SystemName              The name of the files involved in this calculation.
        T_Ref                   The reference temperature used to generate the
                                system geometry.
        T                       The actual system temperature.
        RefStructure            Tells the algorithm which structure to check 
                                for melting/freezing. If a solid structure is selected,
                                the algorithm checks both the liquid and the solid.
        CheckFullTrajectory     Checks the entire trajectory when true, or 
                                only the final time point when false. Note that 
                                when false, the time_to_phase_change output is 
                                set to NaN.
        SaveTrajectory          Saves an .xyz trajectory of the predicted 
                                results when true.
        SaveFeatures            A switch to save the calculated features to a 
                                .pkl file.
        SavePredictions         A switch to saved the calculated structure 
                                predictions to a .pkl file.
        SavePredictionsImage    A switch to save a png image of the calculated 
                                structure predictions.
        InitialRefFrac          What fraction of the system is initiall in the 
                                reference structure?
        RefChangeThreshold      If RefChangeThrehold is a number between 0 and 1, then
                                it corresponds to the minimum change in the 
                                fraction of RefStructure required to indicate 
                                a phase change has completed. If RefChangeThreshold
                                is a number greater than 1, then RefChangeThreshold
                                indicates the absolute number of atoms that must change
                                phase in order to indicate a phase change has completed.
        SlopeThreshold          A second requirement: the change in the asbolute number or
                                fraction per unit time must be smaller than the 
                                absolute value of this threshold for the system 
                                to be considered at the melting point. Units of
                                [Structure Fraction/ps] or [atoms/ps]
        SlopeCheckBegin         At what fraction [0,1) of the total trajectory to 
                                start using data for the slope. This should be 
                                slightly greater than 0 to avoid incorporating
                                changes due to the establishment of equilibrium.
        ML_TimeLength           Controls the particular structure predictor 
                                model to use.
        ML_TimeStep             Controls the particular structure predictor 
                                model to use.
        TimePerFrame            Sets the step size to use between observations
                                in ps, only meaningful when CheckFullTrajectory 
                                is True
        FileType                Sets the structure file type, either gro 
                                or g96
        Verbose                 Sets the level of logging verbosity.    
        Temporal_Cutoff         [Frames] Eliminates short-term fluctuations of solids 
                                in the liquid phase whose timescale is shorter 
                                than this value.
        Voronoi                 Use NN classification models based on Neighbourhoods 
                                determined by Voronoi tesselation
        Qlm_Average             Use NN classification models based on averaged Qlm 
                                parameters
        Prob_Interfacial        Use the probability of structure to classify interfacial 
                                atoms
        Spatial_Reassignment    Re-assign atom classes based on the class of the surrounding
                                neighbourhood
        Spatial_Interfacial     Re-assign atoms to interfacial if Spatial_Interfacial 
                                fraction of surrounding atoms are not all at least one class
        SaveTrajectoryAux       Only applicable when SaveTrajectory is enabled. 
                                Tells the program how much auxiliary data to export.
                                0 = No auxiliary data, 1 = max probability of any one class,
                                2 = also include liquid probability, ... Up to 10.
    
    @author: Hayden Scheiber
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
    import collections

    if CheckFullTrajectory:
        from scipy.stats import linregress
        if SavePredictionsImage:
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
    
    def GetRefStructure(RefStructure):
        rfs = RefStructure.lower()
        if rfs in set(['liquid','liq','l']):
            return 'Liquid'
        elif rfs in set(['rocksalt','rs','r']):
            return 'Rocksalt'
        elif rfs in set(['wurtzite','wz','w']):
            return 'Wurtzite'
        elif rfs in set(['fivefive','5-5','f']):
            return '5-5'
        elif rfs in set(['nias','ns','n']):
            return 'NiAs'
        elif rfs in set(['sphalerite','sphal','s']):
            return 'Sphalerite'
        elif rfs in set(['betabeo','bbeo','b-beo','beta-beo','b']):
            return 'b-BeO'
        elif rfs in set(['antinias','anias','a']):
            return 'AntiNiAs'
        elif rfs in set(['cscl','cc','c']):
            return 'CsCl'
        else:
            return None
    
    # Calculation setup and loading the ML model
    
    # Check that the input directory is valid
    if not os.path.isdir(WorkDir):
        raise NotADirectoryError(WorkDir + ' is not a valid directory.')
    os.chdir(WorkDir)
    
    RefStructure = GetRefStructure(RefStructure)
    
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
    
    if Voronoi:
        Trial_ID = Trial_ID + 'v'
    
    if Qlm_Average:
        Trial_ID = Trial_ID + 'a'
    
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
        model_loc = os.path.join(MLModelDir,'LiX_Structure_Selector_CNN_Model_' + Trial_ID + '.tf')
    elif Version == 2:
        model_loc = os.path.join(MLModelDir,'MX_Structure_Classifier_Model_' + Trial_ID + '.tf')
    
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
    
    if Spatial_Interfacial is not None or Prob_Interfacial is not None:
        labels.append("Interfacial")
        map_dict[n_classes] = "Interfacial"
        map_label[n_classes] = "I"
        
        n_classes += 1
    
    ref_idx = list(map_dict.keys())[list(map_dict.values()).index(RefStructure)]
    liq_idx = list(map_dict.keys())[list(map_dict.values()).index('Liquid')]
    
    if T_Ref is None:
        T_txt = ""
    else:
        T_txt = f"_{T_Ref:.4f}"
    
    trajfile = os.path.join(WorkDir, SystemName + "." + 'trr')
    grofile = os.path.join(WorkDir, SystemName + T_txt + "." + FileType)
    mdpfile = os.path.join(WorkDir, SystemName + "." + 'mdp')
    
    # Check that the required files exist
    if not os.path.isfile(trajfile):
        raise FileNotFoundError('Unable to find trajectory output file: ' + trajfile + ' (file not found).')
    if not os.path.isfile(grofile):
        raise FileNotFoundError('Unable to find initial system geometry file: ' + grofile + ' (file not found).')
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
        t = md.Universe(grofile, trajfile)
    
    traj_timestep = (t.trajectory[-1].time - t.trajectory[0].time)/(t.trajectory.n_frames-1) # ps, time per trajectory frame
    if CheckFullTrajectory:
        min_time = t.trajectory[0].time  # ps
        max_time = t.trajectory[-1].time  # ps
    else:
        min_time = t.trajectory[-1].time  - np.floor(ML_TimeLength/2) # ps
        max_time = t.trajectory[-1].time  - np.floor(ML_TimeLength/2) # ps
    min_step = int(min_time/traj_timestep)
    max_step = int(max_time/traj_timestep)
    
    min_traj_time = t.trajectory[0].time  # ps
    max_traj_time = t.trajectory[-1].time  # ps
    min_traj_step = int(min_traj_time/traj_timestep)
    max_traj_step = int(max_traj_time/traj_timestep)
    
    steps_per_frame = int(max(ML_TimeStep/traj_timestep,1))
    
    # output file name of processed trajectory
    outfile = ML_Name + '_' + Salt + '_' + SystemName + '.xyz'
    
    # List of index points to examine
    steps_per_init_frame = int(TimePerFrame/traj_timestep)
    Traj_starts = list(range(min_step, max_step+1, steps_per_init_frame)) # Steps
    ML_TimeLength_in_steps = int(np.ceil((ML_TimeLength)/traj_timestep)) # number of trajectory steps required to traverse the CNN time slice
    Half_ML_TimeLength_in_steps = int(np.floor(ML_TimeLength_in_steps/2))
    
    # Ensure the last time step is checked when CheckFullTrajectory is active
    if CheckFullTrajectory and Traj_starts[-1] != max_step:
        Traj_starts.append(max_step)
    
    # What is the number of features used in this model (useful for preallocation)
    n_features = 0
    for Ls in L_list:
        n_features = n_features + len(Ls)*(1+int(Include_Wl))
    if Include_ID:
        n_features = n_features + 1
    
    t_slice_len = int(np.ceil((ML_TimeLength_in_steps+1)/steps_per_frame))
    num_traj_starts = len(Traj_starts)
    
    #% Calculate features
    # Loop through the number of time batches
    t_tot = time.time()
    Ql_result = np.empty([n_features,t.atoms.n_atoms])
    Ql_result_traj = np.empty([t_slice_len,t.atoms.n_atoms,n_features])
    Ql_result_all = np.empty([num_traj_starts,t_slice_len,t.atoms.n_atoms,n_features])
    
    if Spatial_Reassignment or Spatial_Interfacial is not None:
        Save_Neighbour = True
        Neighbourlist_list = []
        Neighbourlist_slice = [None] * t_slice_len
        Neighbourlist_slice_prev = [None] * t_slice_len
    else:
        Save_Neighbour = False
    
    logging.info('Generating features...')
    central_idx = int((t_slice_len - 1)/2)
    if Include_ID:
        namelist = t.atoms.names
        Metal_Index = np.array([x not in Metal for x in namelist]).astype(int) # 0 = metal, 1 = halide
        
    if Voronoi:
        voro = freud.locality.Voronoi()
    
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
            
        prev_idx = int(traj_start_idx - np.ceil(steps_per_frame/steps_per_init_frame))
        if prev_idx < 0 or timeless:
            prev_traj_slice = np.array([])
        else:
            prev_init_step = Traj_starts[prev_idx]
            prev_traj_slice = np.array(range(prev_init_step-Half_ML_TimeLength_in_steps, prev_init_step+Half_ML_TimeLength_in_steps+1, steps_per_frame))
            
            if any(prev_traj_slice < min_traj_step):
                ds = min_traj_step - np.min(prev_traj_slice)
                prev_traj_slice = prev_traj_slice + ds
            elif any(prev_traj_slice > max_traj_step):
                ds = max_traj_step - np.max(prev_traj_slice)
                prev_traj_slice = prev_traj_slice + ds
        
        # Calculate Ql values at each time point in the time slice
        for t_idx, ts in enumerate(t.trajectory[traj_slice]):
            
            if traj_slice[t_idx] in prev_traj_slice:
                pridx = list(prev_traj_slice).index(traj_slice[t_idx])
                Ql_result_traj[t_idx] = Ql_result_all[prev_idx,pridx,:,:]
                if Save_Neighbour:
                    Neighbourlist_slice[t_idx] = Neighbourlist_slice_prev[pridx]
                    if t_idx == central_idx:
                        Neighbourlist_list.append(Neighbourlist_slice[t_idx])
                continue
            
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
                
                # Build the neighbour list
                if Voronoi:
                    nlist = voro.compute((box_data, point_data)).nlist
                else:
                    # Construct the neighbour filter and build
                    query_args = dict(mode='nearest', num_neighbors=N_neighbour, exclude_ii=True)
                    nlist = freud.locality.AABBQuery(box_data, point_data).query(point_data, query_args).toNeighborList()
                
                if Save_Neighbour and (Neighbour_idx+1 == len(N_neighbour_list)):
                    Neighbourlist_slice[t_idx] = nlist
                    if t_idx == central_idx:
                        Neighbourlist_list.append(nlist)
                
                for L in L_list[Neighbour_idx]:
                    ql = freud.order.Steinhardt(L,wl=Include_Wl,
                                                wl_normalize=True,
                                                average=Qlm_Average,
                                                weighted=Voronoi)
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
        prev_traj_slice = traj_slice
        if Save_Neighbour:
            Neighbourlist_slice_prev = Neighbourlist_slice.copy()
        logging.info("\rTime Point: {:.2f} ps. Time Elapsed: {:.1f} s. ({:.2f}%, {:.2f} time points/s)".format(
            Init_step*traj_timestep,
            time.time() - t_tot,
            (traj_start_idx+1)*100/num_traj_starts,
            (1)/(time.time() - t_cur)))
    
    if np.shape(Ql_result_all)[1] == 1:
        Ql_result_all = np.squeeze(Ql_result_all,axis=1)
    
    # Optional: Save the calculated features
    if SaveFeatures:
        with open(ML_Name + '_' + SystemName + '_Features.pkl', 'wb') as f:
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
    with tf.device('/cpu:0'):
        predicted_prob = model.predict(Ql_result_all_cat, verbose=int(Verbose), batch_size = 10000)
    
    # Convert to predictions and convert back to time slices
    pred_class = np.array(np.split(np.argmax(predicted_prob, axis=1), num_traj_starts, axis=0))
    
    logging.info('\nFinished Inferring Predictions.')
    
    # Applying various post-processing elements
    predicted_classes = pred_class.copy()
    if Prob_Interfacial is not None:
        max_prob = np.amax(predicted_prob, axis=1)
        max_prob_below_cutoff = np.array(np.split((max_prob <= Prob_Interfacial), num_traj_starts, axis=0))
        predicted_classes[max_prob_below_cutoff] = 9
    
    # Filter out transient predictions
    if Temporal_Cutoff > 0:
        t_tot = time.time()
        logging.info('\nRemoving short-timescale liquid fluctuations...')
        for atm in range(0,np.shape(pred_class)[1]):
            c_atom_trj = (pred_class[:,atm].copy() == 0)
            
            i = 0
            for k, g in groupby(c_atom_trj):
                grp = len(list(g))
                if grp < Temporal_Cutoff:
                    predicted_classes[range(i, i+grp),atm] = 0
                i += grp
        logging.info('\nFinished removing liquid fluctuations. Time elapsed: {:.1f} s'.format(time.time() - t_tot))        
    
    if Spatial_Reassignment or Spatial_Interfacial is not None:
        logging.info('\nReclassifying Atoms...')
        t_tot = time.time()
        for t_idx, traj_start in enumerate(Traj_starts):
            t_cur = time.time()
            nlist = Neighbourlist_list[t_idx]
            cut_idx = np.array([nlist.find_first_index(x) for x in range(1,t.atoms.n_atoms)])
            
            nlist_grps_id = np.split(nlist[:,1],cut_idx)
            nlist_grps = [pred_class[t_idx,x] for x in nlist_grps_id]
            
            mcl = np.array([collections.Counter(x).most_common(1)[0] for x in nlist_grps])
            unmatched_mcl = pred_class[t_idx] != mcl[:,0]
            if Spatial_Reassignment:
                predicted_classes[t_idx,unmatched_mcl] =  mcl[unmatched_mcl,0]
            
            if Spatial_Interfacial is not None:
                predicted_classes[t_idx,mcl[:,1] <= nlist.neighbor_counts*Spatial_Interfacial] = 9
            
            logging.info("\rTime Point: {:.2f} ps. Time Elapsed: {:.1f} s. ({:.2f}%, {:.2f} time points/s)".format(
            Traj_starts[t_idx]*traj_timestep,
            time.time() - t_tot,
            (t_idx+1)*100/num_traj_starts,
            (1)/(time.time() - t_cur)))
        logging.info('\nFinished Reclassifying Atoms. Time elapsed: {:.1f} s'.format(time.time() - t_tot))
    
    if SaveTrajectory:
        if SaveTrajectoryAux > 0:
            Certainty_ts = np.split(np.amax(predicted_prob, axis=1,keepdims=True), num_traj_starts, axis=0)
            Probs_ts = np.split(predicted_prob, num_traj_starts, axis=0)
            aux_dat = np.concatenate((Certainty_ts,Probs_ts), axis=2)[:,:,0:SaveTrajectoryAux]
        else:
            aux_dat = [None] * num_traj_starts
        
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
            W.aux_data = aux_dat[idx]
            W.write(t.atoms)
        W.close()
        
    # Optional: save the predictions
    if SavePredictions:
        with open(ML_Name + '_' + SystemName + '_Predictions.pkl', 'wb') as f:
            pickle.dump([predicted_prob,predicted_classes], f)
    
    # Find the initial fraction of reference structure if not given
    if InitialRefFrac is None and not CheckFullTrajectory:
        t_if = time.time()
        logging.info('\nCalculating Initial Fraction of ' + map_dict[ref_idx] + '...')
        
        traj_slice = np.array(range(min_traj_step, min_traj_step+2*Half_ML_TimeLength_in_steps+1, steps_per_frame))
        
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
                
                # Build the neighbour list
                if Voronoi:
                    nlist = voro.compute((box_data, point_data)).nlist
                else:
                    # Construct the neighbour filter
                    query_args = dict(mode='nearest', num_neighbors=N_neighbour, exclude_ii=True)
                    nlist = freud.locality.AABBQuery(box_data, point_data).query(point_data, query_args).toNeighborList()
                
                if Save_Neighbour and (t_idx == 0 and Neighbour_idx+1 == len(N_neighbour_list)):
                    nlist_init = nlist
                
                for L in L_list[Neighbour_idx]:
                    ql = freud.order.Steinhardt(L,wl=Include_Wl,
                                                wl_normalize=True,
                                                average=Qlm_Average,
                                                weighted=Voronoi)
                    ql_calc = ql.compute(system, neighbors=nlist)
                    
                    # Append order parameters to output
                    Ql_result[fidx] = ql_calc.ql
                    if Include_Wl:
                        Ql_result[fidx+1] = ql_calc.particle_order
                        fidx += 2
                    else:
                        fidx += 1
            
            Ql_result_traj[t_idx] = np.column_stack(Ql_result)
        
        if timeless:
            Ql_initial = Ql_result_traj
        else:
            Ql_initial = Ql_result_traj.transpose(1, 0, 2)
        
        predicted_prob_init = model.predict(Ql_initial, verbose=int(Verbose), batch_size = 10000)
        pred_class_init = np.argmax(predicted_prob_init, axis=1)
        
        # Applying various post-processing elements
        predicted_classes_init = pred_class_init.copy()
        if Prob_Interfacial is not None:
            max_prob = np.amax(predicted_prob_init, axis=1)
            max_prob_below_cutoff = max_prob <= Prob_Interfacial
            predicted_classes_init[max_prob_below_cutoff] = 9
            
        if Spatial_Reassignment or Spatial_Interfacial is not None:
            logging.info('\nReclassifying Atoms...')
            t_tot = time.time()
            
            cut_idx = np.array([nlist_init.find_first_index(x) for x in range(1,t.atoms.n_atoms)])
            
            nlist_grps_id = np.split(nlist_init[:,1],cut_idx)
            nlist_grps = [pred_class_init[x] for x in nlist_grps_id]
            
            mcl = np.array([collections.Counter(x).most_common(1)[0] for x in nlist_grps])
            unmatched_mcl = pred_class_init != mcl[:,0]
            
            if Spatial_Reassignment:
                predicted_classes_init[unmatched_mcl] =  mcl[unmatched_mcl,0]
            
            if Spatial_Interfacial is not None:
                predicted_classes_init[mcl[:,1] <= nlist_init.neighbor_counts*Spatial_Interfacial] = 9
            logging.info('\nFinished Reclassifying Atoms. Time elapsed: {:.1f} s'.format(time.time() - t_tot))
        
        d = hist_laxis(predicted_classes_init, n_classes, [0,n_classes])
        d_norm = d/t.atoms.n_atoms
        InitialRefFrac = d_norm[ref_idx]
        InitialLiqFrac = d_norm[liq_idx]
        if liq_idx != ref_idx:
            logging.info('\nInitial Fraction of ' + map_dict[liq_idx] + 
                         ' Phase is: {:.2f} %.'.format(
                             InitialLiqFrac*100))
        logging.info('\nInitial Fraction of ' + map_dict[ref_idx] + 
                     ' Phase is: {:.2f} %. Time to calculate: {:.1f} s.'.format(
                         InitialRefFrac*100,
                         time.time() - t_if))
    elif InitialRefFrac is None and CheckFullTrajectory:
        t_if = time.time()
        logging.info('\nCalculating Initial Fraction of ' + map_dict[ref_idx] + '...')
        Ql_initial = Ql_result_all_t[0]
        
        predicted_prob_init = model.predict(Ql_initial, verbose=int(Verbose), batch_size = 10000)
        pred_class_init = np.argmax(predicted_prob_init, axis=1)
        
        # Applying various post-processing elements
        predicted_classes_init = pred_class_init.copy()
        if Prob_Interfacial is not None:
            max_prob = np.amax(predicted_prob_init, axis=1)
            max_prob_below_cutoff = max_prob <= Prob_Interfacial
            predicted_classes_init[max_prob_below_cutoff] = 9
            
        if Spatial_Reassignment or Spatial_Interfacial is not None:
            logging.info('\nReclassifying Atoms...')
            t_tot = time.time()
            
            nlist_init = Neighbourlist_list[0]
            cut_idx = np.array([nlist_init.find_first_index(x) for x in range(1,t.atoms.n_atoms)])
            
            nlist_grps_id = np.split(nlist_init[:,1],cut_idx)
            nlist_grps = [pred_class_init[x] for x in nlist_grps_id]
            
            mcl = np.array([collections.Counter(x).most_common(1)[0] for x in nlist_grps])
            unmatched_mcl = pred_class_init != mcl[:,0]
            
            if Spatial_Reassignment:
                predicted_classes_init[unmatched_mcl] =  mcl[unmatched_mcl,0]
            
            if Spatial_Interfacial is not None:
                predicted_classes_init[mcl[:,1] <= nlist_init.neighbor_counts*Spatial_Interfacial] = 9
            logging.info('\nFinished Reclassifying Atoms. Time elapsed: {:.1f} s'.format(time.time() - t_tot))
        
        d = hist_laxis(predicted_classes_init, n_classes, [0,n_classes])
        d_norm = d/t.atoms.n_atoms
        InitialRefFrac = d_norm[ref_idx]
        InitialLiqFrac = d_norm[liq_idx]
        if liq_idx != ref_idx:
            logging.info('\nInitial Fraction of ' + map_dict[liq_idx] + 
                         ' Phase is: {:.2f} %.'.format(
                             InitialLiqFrac*100))
        logging.info('\nInitial Fraction of ' + map_dict[ref_idx] + 
                     ' Phase is: {:.2f} %. Time to calculate: {:.1f} s.'.format(
                         InitialRefFrac*100,
                         time.time() - t_if))
    else:
        InitialLiqFrac = 1 - InitialRefFrac
    
    # Finally, generate a series of line plots: partitions vs time
    y_max = 100
    
    # Separate into metal and halide
    Metal_Ind = (Metal_Index == 0) # 0 = metal, 1 = halide
    Halide_Ind = (Metal_Index == 1) # 0 = metal, 1 = halide
    
    # Generate the histogram along one dimension: combined metal/halide
    d = hist_laxis(predicted_classes, n_classes, [0,n_classes])
    d_norm = d/t.atoms.n_atoms
    final_ref_frac = d_norm[-1,ref_idx]
    final_liq_frac = d_norm[-1,liq_idx]
    x = np.array([i * traj_timestep + min_traj_time for i in Traj_starts])
    
    if CheckFullTrajectory and SavePredictionsImage:
        
        # Generate the histogram along one dimension: separate for metal and halide
        d_metal = hist_laxis(predicted_classes[:,Metal_Ind], n_classes, [0,n_classes])
        d_metal_norm = d_metal/sum(Metal_Ind)
        d_halide = hist_laxis(predicted_classes[:,Halide_Ind], n_classes, [0,n_classes])
        d_halide_norm = d_halide/sum(Halide_Ind)
        
        # colors
        cols = sns.color_palette("muted",n_classes)
        
        # Set up the matplotlib figure
        
        fig, ax = plt.subplots()
        
        for idx,y_dat_metal in enumerate(d_metal_norm.T):
            y_dat_halide = d_halide_norm.T[idx,:]
            plt.plot(x,y_dat_metal*100, label=labels[idx], color=cols[idx], alpha=0.8, linewidth=2, linestyle='solid')
            plt.plot(x,y_dat_halide*100, label=None, color=cols[idx], alpha=0.8, linewidth=2, linestyle='dashed')
        
        if T is None:
            T_txt1 = ""
            T_txt2 = ""
        else:
            T_txt1 = f" T = {T:.4f} K"
            T_txt2 = f"T_{T:.4f}"
        
        plt.title('[' + ML_Name + ']: ' + SystemName + T_txt1)
        plt.grid()
        plt.xlabel('Trajectory Time [ps]')
        plt.ylabel('Fraction in Structure [%]')
        ax.set_yticks(np.arange(0, 110, step=10))
        plt.legend()
        if y_max < 100:
            plt.ylim([0,y_max])
        
        fig.savefig(ML_Name + '_' + Salt + '_' + SystemName + T_txt2 + '.png', dpi=300, format='png')
        plt.close(fig)
    
    # Determine the fraction of liquid, find first time step where liquid becomes less than or greater than the given threshold 
    system_melting = False
    system_freezing = False
    system_frzing_alt = False
            
    # Leave SlopeThreshold and InitialRefFrac as fractions or convert them to absolute numbers
    if RefChangeThreshold <= 1:
        Threshold_Upper_sol = min(InitialRefFrac + RefChangeThreshold,0.95) # Above this number indicates melting/freezing [solid]
        Threshold_Lower_sol = max(InitialRefFrac - RefChangeThreshold,0.05) # Below this number indicates freezing/melting [solid]
        Threshold_Upper_liq = min(InitialLiqFrac + RefChangeThreshold,0.95) # Above this number indicates melting/freezing [liquid]
        Threshold_Lower_liq = max(InitialLiqFrac - RefChangeThreshold,0.05) # Below this number indicates freezing/melting [liquid]
        ref_fraction = d_norm[:,ref_idx]
        liq_fraction = d_norm[:,liq_idx]
    else: # Convert to absolute units
        ref_fraction = d[:,ref_idx]
        liq_fraction = d[:,liq_idx]
        SlopeThreshold = SlopeThreshold*t.atoms.n_atoms # [Structure Fraction/ps] -> [Atoms/ps]
        InitialRefFrac = InitialRefFrac*t.atoms.n_atoms # Initial fraction -> Initial Num atoms
        InitialLiqFrac = InitialLiqFrac*t.atoms.n_atoms # Initial fraction -> Initial Num atoms
        Threshold_Upper_sol = min(InitialRefFrac + RefChangeThreshold,t.atoms.n_atoms*0.95) # Above this number indicates melting/freezing [liquid/solid]
        Threshold_Lower_sol = max(InitialRefFrac - RefChangeThreshold,t.atoms.n_atoms*0.05) # Below this number indicates freezing/melting [liquid/solid]
        Threshold_Upper_liq = min(InitialLiqFrac + RefChangeThreshold,t.atoms.n_atoms*0.95) # Above this number indicates melting/freezing [liquid/solid]
        Threshold_Lower_liq = max(InitialLiqFrac - RefChangeThreshold,t.atoms.n_atoms*0.05) # Below this number indicates freezing/melting [liquid/solid]
        
    if ref_idx == 0: # Reference Structure is Liquid
        is_frozen = ref_fraction <= Threshold_Lower_liq
        is_melted = ref_fraction >= Threshold_Upper_liq
        is_frozen_alt = False # unable to check if only looking at liquid
        # Check slope of last 90% of trajectory (fit linear regression)
        if CheckFullTrajectory:
            x1 = x[x >= max_time*SlopeCheckBegin] # [ps] time points
            y1 = ref_fraction[x >= max_time*SlopeCheckBegin] # fraction (or absolute number) of reference structure
            regmodel = linregress(x1, y1)
            slope = regmodel.slope
            if slope > SlopeThreshold:
                system_melting = True
                est_time_to_phase_change = (Threshold_Upper - y1[0])/slope + x1[0]
            elif slope < -SlopeThreshold:
                system_freezing = True
                est_time_to_phase_change = (Threshold_Lower - y1[0])/slope + x1[0]

    else: # Reference structure is a solid structure
        is_frozen_sol = ref_fraction >= Threshold_Upper_sol
        is_melted_sol = ref_fraction <= Threshold_Lower_sol
        
        is_frozen_liq = liq_fraction <= Threshold_Lower_liq
        is_melted_liq = liq_fraction >= Threshold_Upper_liq
        
        is_frozen = is_frozen_sol | is_frozen_liq
        is_melted = is_melted_sol | is_melted_liq
        
        alt_fraction = 1 - (liq_fraction + ref_fraction)
        is_frozen_alt = alt_fraction >= RefChangeThreshold
        
        # Check slope of last $SlopeCheckBegin% of trajectory
        if CheckFullTrajectory:
            x1 = x[x >= max_time*SlopeCheckBegin] # [ps] time points
            y1 = ref_fraction[x >= max_time*SlopeCheckBegin] # fraction (or absolute number) of reference structure
            y2 = liq_fraction[x >= max_time*SlopeCheckBegin]
            y3 = alt_fraction[x >= max_time*SlopeCheckBegin]
            regmodel_sol = linregress(x1, y1)
            regmodel_liq = linregress(x1, y2)
            regmodel_alt = linregress(x1, y3)
            slope_sol = regmodel_sol.slope
            slope_liq = regmodel_liq.slope
            slope_alt = regmodel_alt.slope
            if slope_sol > SlopeThreshold:
                system_freezing = True
                est_time_to_phase_change = (Threshold_Upper_sol - y1[0])/slope_sol + x1[0]
            elif slope_liq < -SlopeThreshold:
                system_freezing = True
                est_time_to_phase_change = (Threshold_Upper_liq - y2[0])/slope_liq + x1[0]
            elif slope_sol < -SlopeThreshold:
                system_melting = True
                est_time_to_phase_change = (Threshold_Lower_sol - y1[0])/slope_sol + x1[0]
            elif slope_liq > SlopeThreshold:
                system_melting = True
                est_time_to_phase_change = (Threshold_Lower_liq - y2[0])/slope_liq + x1[0]
            elif slope_alt > SlopeThreshold:
                system_frzing_alt = True
                est_time_to_phase_change = (RefChangeThreshold - y3[0])/slope_sol + x1[0]

    # If system is passed the threshold into frozen or melting, use that
    system_froze = any(is_frozen)
    system_melted = any(is_melted)
    system_frz_alt = any(is_frozen_alt) | system_frzing_alt
    
    # Otherwise, if the system is not yet frozen or melting, check if it is freezing/melting based on slope
    if not system_froze and not system_melted:
        system_froze = system_freezing
        system_melted = system_melting
    
    # If both freezing and melting is detected, check which occurs later
    if system_froze and system_melted:
        if max(x[is_frozen]) > max(x[is_melted]):
            system_froze = True
            system_melted = False
        else:
            system_froze = False
            system_melted = True
    
    # If liquid has neither fully melted nor fully frozen, report back
    if not system_froze and not system_melted and not system_frz_alt:
        time_to_phase_change = np.nan
    elif system_froze:
        if any(is_frozen):
            idx_phase_change = np.where(is_frozen)[0][0]
            time_to_phase_change = x[idx_phase_change]
        else:
            time_to_phase_change = est_time_to_phase_change
        
    elif system_melted:
        if any(is_melted):
            idx_phase_change = np.where(is_melted)[0][0]
            time_to_phase_change = x[idx_phase_change]
        else:
            time_to_phase_change = est_time_to_phase_change
    elif system_frz_alt:
        if any(is_frozen_alt):
            idx_phase_change = np.where(is_frozen_alt)[0][0]
            time_to_phase_change = x[idx_phase_change]
        else:
            time_to_phase_change = est_time_to_phase_change
        
    if not CheckFullTrajectory:
        time_to_phase_change = np.nan
    
    return [system_froze,system_melted,time_to_phase_change,final_ref_frac,final_liq_frac,system_frz_alt]