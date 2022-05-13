# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 18:54:17 2020

QlEdit should be an number
All lengths should be nm
All times should be ps

trajfile = 'C:\\Users\\Hayden\\Documents\\Patey_Lab\\Results\\LiI\\Step_W_JC_NPT\\Step_W_JC_NPT.trr'
grofile = 'C:\\Users\\Hayden\\Documents\\Patey_Lab\\Results\\LiI\\Step_W_JC_NPT\\Step_W_JC_NPT.gro'
Metal = 'Li'
Halide = 'I'
Ref_Type = 'All-All'
Start_Points = np.array([1490, 3490, 5490, 7490, 9490, 11490, 13490, 15490, 17490, 19490, 21490, 23490, 25490, 27490, 29490])
End_Points = np.array([1510,  3510,  5510,  7510,  9510, 11510, 13510, 15510, 17510, 19510, 21510, 23510, 25510, 27510, 29510])
TimePerFrame = 1
First_Shell_Switch = False
Second_Shell_Switch = True
SelRadius = False
NearestN = 0
QlEditMin = 0
QlEdit = 4
pbc_on = True
L_list = np.array([2, 4, 6, 8, 10])

@author: Hayden
"""
def Calculate_Ql(trajfile,grofile,Metal,Halide,Ref_Type,Start_Points,
                    End_Points,TimePerFrame,First_Shell_Switch,Second_Shell_Switch,
                    SelRadius,NearestN,QlEditMin,QlEdit,pbc_on,L_list):
    
    import MDAnalysis as md
    import numpy as np
    import PySimpleGUI as sg
    import freud
    import warnings
    
    def to_freud_box(dims):
        a = dims[0]/10
        b = dims[1]/10
        c = dims[2]/10
        alpha = np.deg2rad(dims[3])
        beta = np.deg2rad(dims[4])
        gamma = np.deg2rad(dims[5])
    
        Omega = a*b*c*np.sqrt(1 - (np.cos(alpha)**2) - (np.cos(beta)**2) - (np.cos(gamma)**2) + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        cx = c*np.cos(beta)
        cy = c*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        cz = Omega/(a*b*np.sin(gamma))
    
        Transform_Matrix = [[a,             b*np.cos(gamma),      cx],
                            [0,             b*np.sin(gamma),      cy],
                            [0,             0,                    cz]]
        
        return freud.box.Box.from_matrix(Transform_Matrix)
    
    def baseline_als(y, lam, p, niter=10):
        L = len(y)
        D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
        w = np.ones(L)
        for i in range(niter):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.transpose())
            z = spsolve(Z, w*y)
            w = p * (y > z) + (1-p) * (y < z)
        return z
    
    N_L = len(L_list)
    
    if SelRadius:
        Max_R = QlEdit/10; # Ql/Wl radius of nearest neighbours to search in nm
        Min_R = QlEditMin/10;
    elif First_Shell_Switch:
        Min_R = 0;
    elif not Second_Shell_Switch:
        Max_R = np.nan
        Min_R = np.nan
        N_Neigh = int(NearestN)
        
    SelNeigh = not (SelRadius | First_Shell_Switch | Second_Shell_Switch)
    
    RDF_Req = First_Shell_Switch | Second_Shell_Switch
    
    [ref,sel] = Ref_Type.split('-', 1)
    
    if np.isscalar(Start_Points):
        Start_Points = np.array([Start_Points],'float32')
        End_Points = np.array([End_Points],'float32')
    
    N_SP = len(Start_Points)
    
    # Load the trajectory
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t = md.Universe(grofile, trajfile)
    
    # Pre-allocate array for later
    Ql_result = np.empty([N_L,t.atoms.n_atoms])*np.nan
    
    # trajectory time step in ps
    timestep = (t.trajectory[-1].time - t.trajectory[0].time)/(t.trajectory.n_frames-1)
    
    # number of trajectory steps between each frame that we keep
    steps_per_frame = int(TimePerFrame/timestep)
    max_time = t.trajectory[-1].time
    min_time = t.trajectory[0].time
    
    # Pre calculate total number of order parameters
    time_slices = []
    tot_calcs = 0
    for idx, Spoints in enumerate(Start_Points):
        # index of the start point within trajectory
        Start_Step = int(np.max([Spoints,min_time])/timestep)
        
        # index of the end point within trajectory
        End_Step = int(np.min([End_Points[idx],max_time])/timestep)
        
        time_slices.append(list(range(Start_Step, End_Step+1, steps_per_frame)))
        tot_calcs += len(time_slices[idx])
        
    sg.one_line_progress_meter('', 0.1, tot_calcs, 'key',
                           'Calculating Order Parameters',no_titlebar=False,
                           orientation="h")
    sbidx = 0;
    
    Ql_result_TimePoints = []
    locmin_prev = [0,0]
    for idx in range(0,N_SP):

        time_slice = time_slices[idx]
        n_time_slice = len(time_slice)
                
        if RDF_Req:
            from scipy.ndimage import gaussian_filter1d
            from scipy.signal import argrelextrema
            from scipy import sparse
            from scipy.sparse.linalg import spsolve
            
            # initialize rdf
            rdf = freud.density.RDF(200, 0.8)
            
            sg.one_line_progress_meter('', sbidx+0.5, tot_calcs, 'key',
                                   'Calculating Shells',no_titlebar=False,
                                   orientation="h")
            
            # compute rdf with Freud
            for ts in t.trajectory[time_slice]:
                if pbc_on:
                    box_data = to_freud_box(ts.dimensions)
                else:
                    # If pbc is off, place particles in ~infinite box
                    box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
                    
                pos_data = t.atoms.positions/10
                rdf.compute(system=(box_data, pos_data), reset=False)
            
            rdf_x = rdf.bin_centers
            rdf_y = rdf.rdf
            
            # Smooth rdf
            rdf_yp = gaussian_filter1d(rdf_y, 3)
            # Remove baseline
            rdf_ypp = baseline_als(rdf_yp, 1e5, 0.5, niter=10)
    
            # import matplotlib.pyplot as plt
            # fig, ax = plt.subplots(1, 1, figsize=(12, 8))
            # ax.plot(rdf_x, rdf_y, label='RDF')
            # ax.plot(rdf_x, rdf_yp, label='RDF')
            # ax.plot(rdf_x, rdf_yp - rdf_ypp, label='RDF')
    
            # get local minima
            locmin_idx = argrelextrema(rdf_yp - rdf_ypp, np.less)
            locmin = rdf_x[locmin_idx[0]]
            locmin = np.delete(locmin,0)
            
            # get inflection points
            smooth_d2 = np.gradient(np.gradient(rdf_yp))
            infls_idx = np.where(np.diff(np.sign(smooth_d2)))[0]
            infls = rdf_x[infls_idx]
            
            if len(locmin) < 2 & Second_Shell_Switch:
                if len(locmin) < 1:
                    locmin = [0, 0]
                
                if len(infls) < 4:
                    if idx == 0:
                        import warnings
                        warnings.warn("Unable to locate second shell")
                        return np.nan
                    else:
                        locmin[1] = locmin_prev[1]
                else:
                    locmin[1] = infls[3]
                    
            elif len(locmin) < 1 & (First_Shell_Switch | Second_Shell_Switch):
                if len(infls) < 2:
                    if idx == 0:
                        import warnings
                        warnings.warn("Unable to locate first shell")
                        return np.nan
                    else:
                        locmin[0] = locmin_prev[0]
                else:
                    locmin[0] = infls[1]
                    
                    if First_Shell_Switch:
                        First_Shell = locmin[1] # nm
                    
                    if Second_Shell_Switch:
                        First_Shell = locmin[0] # nm
                        Second_Shell = locmin[1] # nm
            
            if First_Shell_Switch:
                First_Shell = locmin[0] # nm
            
            if Second_Shell_Switch:
                First_Shell = locmin[0] # nm
                Second_Shell = locmin[1] # nm
    
            if First_Shell_Switch & Second_Shell_Switch:
                Max_R = Second_Shell
                
            elif First_Shell_Switch:
                Max_R = First_Shell
                
            elif Second_Shell_Switch:
                Max_R = Second_Shell
                Min_R = First_Shell
            
            locmin_prev = locmin
    
        sg.one_line_progress_meter('', sbidx+0.5, tot_calcs, 'key',
                               'Calculating Order Parameters',no_titlebar=False,
                               orientation="h")
        if SelNeigh:
            query_args = dict(mode='nearest', num_neighbors=N_Neigh, exclude_ii=True)
        else:
            query_args = dict(mode='ball', r_min=Min_R, r_max=Max_R, exclude_ii=True)
            
        
        Ql_result_TimeWindow = np.empty([n_time_slice,N_L,t.atoms.n_atoms])*np.nan
        for t_idx, ts in enumerate(t.trajectory[time_slice]):
            
            if pbc_on:
                box_data = to_freud_box(t.atoms.dimensions)
            else:
                # If pbc is off, place particles in ~infinite box
                box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
                        
            query_data = t.atoms.positions/10 # in nm
            
            if ref == Metal:
                query_data[t.select_atoms("name " + Halide).indices] = np.nan
            elif ref == Halide:
                query_data[t.select_atoms("name " + Metal).indices] = np.nan
        
            #filter_idx = np.logical_and(query_data[:,2] > 0.1,query_data[:,2] < 0.7)
            #query_data_filtered = query_data[filter_idx,:]
        
        
            point_data = t.atoms.positions/10 # in nm
            system = [box_data,point_data]
            if sel == Metal:
                point_data[t.select_atoms("name " + Halide).indices] = np.nan
            elif sel == Halide:
                point_data[t.select_atoms("name " + Metal).indices] = np.nan
            
            nlist = freud.locality.AABBQuery(box_data, point_data).query(query_data, query_args).toNeighborList()
            
            for L_idx,L in enumerate(L_list):
                
                ql = freud.order.Steinhardt(L)
                ql_calc = ql.compute(system, neighbors=nlist)
                
                # Append to output
                Ql_result[L_idx] = ql_calc.ql
            
            Ql_result_TimeWindow[t_idx] = Ql_result
            
            sbidx += 1
            sg.one_line_progress_meter('', sbidx, tot_calcs, 'key',
                                   'Calculating Order Parameters',no_titlebar=False,
                                   orientation="h")
            
        Ql_result_TimePoints.append(np.concatenate(Ql_result_TimeWindow,axis=1))
        
    return Ql_result_TimePoints

def Calculate_Wl(trajfile,grofile,Metal,Halide,Ref_Type,Start_Points,
                    End_Points,TimePerFrame,First_Shell_Switch,Second_Shell_Switch,
                    SelRadius,NearestN,QlEditMin,QlEdit,pbc_on,L_list):
    
    import MDAnalysis as md
    import numpy as np
    import PySimpleGUI as sg
    import freud
    
    def to_freud_box(dims):
        a = dims[0]/10
        b = dims[1]/10
        c = dims[2]/10
        alpha = np.deg2rad(dims[3])
        beta = np.deg2rad(dims[4])
        gamma = np.deg2rad(dims[5])
    
        Omega = a*b*c*np.sqrt(1 - (np.cos(alpha)**2) - (np.cos(beta)**2) - (np.cos(gamma)**2) + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        cx = c*np.cos(beta)
        cy = c*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        cz = Omega/(a*b*np.sin(gamma))
    
        Transform_Matrix = [[a,             b*np.cos(gamma),      cx],
                            [0,             b*np.sin(gamma),      cy],
                            [0,             0,                    cz]]
        
        return freud.box.Box.from_matrix(Transform_Matrix)
    
    def baseline_als(y, lam, p, niter=10):
        L = len(y)
        D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
        w = np.ones(L)
        for i in range(niter):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.transpose())
            z = spsolve(Z, w*y)
            w = p * (y > z) + (1-p) * (y < z)
        return z
    
    if SelRadius:
        Max_R = QlEdit/10; # Ql/Wl radius of nearest neighbours to search in nm
        Min_R = QlEditMin/10;
    elif First_Shell_Switch:
        Min_R = 0;
    elif not Second_Shell_Switch:
        Max_R = np.nan
        Min_R = np.nan
        N_Neigh = int(NearestN)
        
    SelNeigh = not (SelRadius | First_Shell_Switch | Second_Shell_Switch)
    
    RDF_Req = First_Shell_Switch | Second_Shell_Switch
    
    [ref,sel] = Ref_Type.split('-', 1)
    
    if np.isscalar(Start_Points):
        Start_Points = np.array([Start_Points],'float32')
        End_Points = np.array([End_Points],'float32')
    
    N_SP = len(Start_Points)
    
    # Load the trajectory
    t = md.Universe(grofile, trajfile)
    timestep = t.trajectory.dt # ps
    steps_per_frame = int(TimePerFrame/timestep)
    max_time = t.trajectory[-1].time
    min_time = t.trajectory[0].time
    
    # Pre calculate total number of order parameters for the progress bar
    time_slices = []
    tot_calcs = 0
    for idx, Spoints in enumerate(Start_Points):
        # index of the start point within trajectory
        Start_Step = int(np.max([Spoints,min_time])/timestep)
        
        # index of the end point within trajectory
        End_Step = int(np.min([End_Points[idx],max_time])/timestep)
        
        time_slices.append(list(range(Start_Step, End_Step+1, steps_per_frame)))
        tot_calcs += len(time_slices[idx])
        
    sg.one_line_progress_meter('', 0.1, tot_calcs, 'key',
                           'Calculating Order Parameters',no_titlebar=False,
                           orientation="h")
    sbidx = 0;
    
    Ql_result_TimePoints = []
    Wl_result_TimePoints = []
    locmin_prev=[0,0]
    for idx in range(0,N_SP):
        
        time_slice = time_slices[idx]
        
        if RDF_Req:
            from scipy.ndimage import gaussian_filter1d
            from scipy.signal import argrelextrema
            from scipy import sparse
            from scipy.sparse.linalg import spsolve
            
            # initialize rdf
            rdf = freud.density.RDF(200, 0.8)
            
            # compute rdf with Freud
            for ts in t.trajectory[time_slice]:
                if pbc_on:
                    box_data = to_freud_box(t.atoms.dimensions)
                else:
                    # If pbc is off, place particles in ~infinite box
                    box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
                    
                pos_data = t.atoms.positions/10
                rdf.compute(system=(box_data, pos_data), reset=False)
            
            rdf_x = rdf.bin_centers
            rdf_y = getattr(rdf, 'rdf')
            
            # Smooth rdf
            rdf_yp = gaussian_filter1d(rdf_y, 6)
            # Remove baseline
            rdf_ypp = baseline_als(rdf_yp, 1e5, 0.5, niter=10)
    
            # get local minima
            locmin_idx = argrelextrema(rdf_yp - rdf_ypp, np.less)
            locmin = rdf_x[locmin_idx[0]]
            locmin = np.delete(locmin,0)
            
            # get inflection points
            smooth_d2 = np.gradient(np.gradient(rdf_yp))
            infls_idx = np.where(np.diff(np.sign(smooth_d2)))[0]
            infls = rdf_x[infls_idx]
            
            if len(locmin) < 2 & Second_Shell_Switch:
                if len(locmin) < 1:
                    locmin = [0, 0]
                
                if len(infls) < 4:
                    
                    if idx == 0:
                        import warnings
                        warnings.warn("Unable to locate second shell")
                        return np.nan
                    else:
                        locmin[1] = locmin_prev[1]
                else:
                    locmin[1] = infls[3]
                    
            elif len(locmin) < 1 & (First_Shell_Switch | Second_Shell_Switch):
                if len(infls) < 2:
                    if idx == 0:
                        import warnings
                        warnings.warn("Unable to locate second shell")
                        return np.nan
                    else:
                        locmin[0] = locmin_prev[0]
                else:
                    locmin[0] = infls[1]
                    
                    if First_Shell_Switch:
                        First_Shell = locmin[1] # nm
                    
                    if Second_Shell_Switch:
                        First_Shell = locmin[0] # nm
                        Second_Shell = locmin[1] # nm
            
            if First_Shell_Switch:
                First_Shell = locmin[0] # nm
            
            if Second_Shell_Switch:
                First_Shell = locmin[0] # nm
                Second_Shell = locmin[1] # nm
    
            if First_Shell_Switch & Second_Shell_Switch:
                Max_R = Second_Shell
                
            elif First_Shell_Switch:
                Max_R = First_Shell
                
            elif Second_Shell_Switch:
                Max_R = Second_Shell
                Min_R = First_Shell
            
            locmin_prev = locmin
    
        if SelNeigh:
            query_args = dict(mode='nearest', num_neighbors=N_Neigh, exclude_ii=True)
        else:
            query_args = dict(mode='ball', r_min=Min_R, r_max=Max_R, exclude_ii=True)
        
        Ql_result_TimeWindow = []
        Wl_result_TimeWindow = []
        for ts in t.trajectory[time_slice]:
            
            if pbc_on:
                box_data = to_freud_box(t.atoms.dimensions)
            else:
                # If pbc is off, place particles in ~infinite box
                box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
                        
            query_data = t.atoms.positions/10 # in nm
            if ref == Metal:
                query_data[t.select_atoms("name " + Halide).indices] = np.nan
            elif ref == Halide:
                query_data[t.select_atoms("name " + Metal).indices] = np.nan
        
            point_data = t.atoms.positions/10 # in nm
            if sel == Metal:
                point_data[t.select_atoms("name " + Halide).indices] = np.nan
            elif sel == Halide:
                point_data[t.select_atoms("name " + Metal).indices] = np.nan
            
            nlist = freud.locality.AABBQuery(box_data, point_data).query(query_data, query_args).toNeighborList()
            
            system = [box_data,t.atoms.positions/10]
            
            Ql_result = []
            Wl_result = []
            for L_idx,L in enumerate(L_list):
                
                ql = freud.order.Steinhardt(L,wl=True,wl_normalize=True)
                ql_calc = ql.compute(system, neighbors=nlist)
                
                # Append to output
                Ql_result.append(ql_calc.ql)
                Wl_result.append(ql_calc.particle_order)

            Ql_result_TimeWindow.append(Ql_result)
            Wl_result_TimeWindow.append(Wl_result)
            sbidx += 1
            sg.one_line_progress_meter('', sbidx, tot_calcs, 'key',
                                   'Calculating Order Parameters',no_titlebar=False,
                                   orientation="h")
            
        Ql_result_TimePoints.append(np.concatenate(Ql_result_TimeWindow,axis=1))
        Wl_result_TimePoints.append(np.concatenate(Wl_result_TimeWindow,axis=1))
        
    return [Wl_result_TimePoints,Ql_result_TimePoints]

def Calculate_RDF(trajfile,grofile,Metal,Halide,Ref_Type,Start_Points,
                    End_Points,TimePerFrame,pbc_on):
    
    import MDAnalysis as md
    import numpy as np
    import PySimpleGUI as sg
    import freud
    
    def to_freud_box(dims):
        a = dims[0]/10
        b = dims[1]/10
        c = dims[2]/10
        alpha = np.deg2rad(dims[3])
        beta = np.deg2rad(dims[4])
        gamma = np.deg2rad(dims[5])
    
        Omega = a*b*c*np.sqrt(1 - (np.cos(alpha)**2) - (np.cos(beta)**2) - (np.cos(gamma)**2) + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        cx = c*np.cos(beta)
        cy = c*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        cz = Omega/(a*b*np.sin(gamma))
    
        Transform_Matrix = [[a,             b*np.cos(gamma),      cx],
                            [0,             b*np.sin(gamma),      cy],
                            [0,             0,                    cz]]
        
        return freud.box.Box.from_matrix(Transform_Matrix)
            
    [ref,sel] = Ref_Type.split('-', 1)
    
    if np.isscalar(Start_Points):
        Start_Points = np.array([Start_Points],'float32')
        End_Points = np.array([End_Points],'float32')
    
    N_SP = len(Start_Points)
    
    # Load the trajectory
    t = md.Universe(grofile, trajfile)
    timestep = t.trajectory.dt # ps
    steps_per_frame = int(TimePerFrame/timestep)
    max_time = t.trajectory[-1].time
    min_time = t.trajectory[0].time
    
    Max_R = np.min([1.2,np.min(t.atoms.dimensions[0:3])/10])
    
    # Calculate number of bins
    rdf_bins = Max_R/0.002
    
    query_args = dict(mode='ball', r_min=0, r_max=Max_R, exclude_ii=True)
    
    # Pre calculate total number of RDF frames for the progress bar
    time_slices = []
    tot_calcs = 0
    for idx, Spoints in enumerate(Start_Points):
        # index of the start point within trajectory
        Start_Step = int(np.max([Spoints,min_time])/timestep)
        
        # index of the end point within trajectory
        End_Step = int(np.min([End_Points[idx],max_time])/timestep)
        
        time_slices.append(list(range(Start_Step, End_Step+1, steps_per_frame)))
        tot_calcs += len(time_slices[idx])
        
    sg.one_line_progress_meter('', 1, tot_calcs, 'key',
                           'Calculating RDF(s)',no_titlebar=False,
                           orientation="h")
    sbidx = 0;
    
    RDF_result_TimePoints = []
    CDF_result_TimePoints = []
    for idx in range(0,N_SP):

        time_slice = time_slices[idx]
        
        # initialize rdf
        rdf = freud.density.RDF(bins=rdf_bins, r_max=Max_R, r_min=0)
        
        # Loop through time steps of the current time slice
        for ts in t.trajectory[time_slice]:
            
            if pbc_on:
                box_data = to_freud_box(t.atoms.dimensions)
            else:
                # If pbc is off, place particles in ~infinite box
                box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
            
            if ref == Metal:
                query_data = t.select_atoms("name " + Metal).atoms.positions/10
            elif ref == Halide:
                query_data = t.select_atoms("name " + Halide).atoms.positions/10
            else:
                query_data = t.atoms.positions/10
            
            if sel == Metal:
                pos_data = t.select_atoms("name " + Metal).atoms.positions/10
            elif sel == Halide:
                pos_data = t.select_atoms("name " + Halide).atoms.positions/10
            else:
                pos_data = t.atoms.positions/10
            
            # Create RDF neighbourlist
            nlist = freud.locality.AABBQuery(box_data, pos_data).query(query_data, query_args).toNeighborList()
            
            # Compute the RDF using the neighbourlist
            rdf.compute(system=(box_data, pos_data), neighbors=nlist, reset=False)
            
            sbidx += 1
            sg.one_line_progress_meter('', sbidx, tot_calcs, 'key',
                                   'Calculating RDF(s)',no_titlebar=False,
                                   orientation="h")

        RDF_result_TimePoints.append(np.array([rdf.bin_centers,rdf.rdf]))
        CDF_result_TimePoints.append(np.array([rdf.bin_centers,rdf.n_r]))
        
    return [np.array(RDF_result_TimePoints),np.array(CDF_result_TimePoints)]

def Calculate_NN(trajfile,grofile,Metal,Halide,Ref_Type,Start_Points,
                    End_Points,TimePerFrame,First_Shell_Switch,Second_Shell_Switch,
                    QlEditMin,QlEdit,pbc_on):
    
    import MDAnalysis as md
    import numpy as np
    import PySimpleGUI as sg
    import freud
    
    def to_freud_box(dims):
        a = dims[0]/10
        b = dims[1]/10
        c = dims[2]/10
        alpha = np.deg2rad(dims[3])
        beta = np.deg2rad(dims[4])
        gamma = np.deg2rad(dims[5])
    
        Omega = a*b*c*np.sqrt(1 - (np.cos(alpha)**2) - (np.cos(beta)**2) - (np.cos(gamma)**2) + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        cx = c*np.cos(beta)
        cy = c*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        cz = Omega/(a*b*np.sin(gamma))
    
        Transform_Matrix = [[a,             b*np.cos(gamma),      cx],
                            [0,             b*np.sin(gamma),      cy],
                            [0,             0,                    cz]]
        
        return freud.box.Box.from_matrix(Transform_Matrix)
    
    def baseline_als(y, lam, p, niter=10):
        L = len(y)
        D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
        w = np.ones(L)
        for i in range(niter):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.transpose())
            z = spsolve(Z, w*y)
            w = p * (y > z) + (1-p) * (y < z)
        return z
    
    if First_Shell_Switch:
        Min_R = 0;
            
    RDF_Req = First_Shell_Switch | Second_Shell_Switch
    
    if not RDF_Req:
        Max_R = QlEdit/10; # Radius of nearest neighbours to search in nm
        Min_R = QlEditMin/10;
    elif First_Shell_Switch:
        Min_R = 0;
    
    [ref,sel] = Ref_Type.split('-', 1)
    
    if np.isscalar(Start_Points):
        Start_Points = np.array([Start_Points],'float32')
        End_Points = np.array([End_Points],'float32')
    
    N_SP = len(Start_Points)
    
    # Load the trajectory
    t = md.Universe(grofile, trajfile)
    timestep = t.trajectory.dt # ps
    steps_per_frame = int(TimePerFrame/timestep)
    max_time = t.trajectory[-1].time
    min_time = t.trajectory[0].time
    
    # Pre calculate total number of order parameters for the progress bar
    time_slices = []
    tot_calcs = 0
    for idx, Spoints in enumerate(Start_Points):
        # index of the start point within trajectory
        Start_Step = int(np.max([Spoints,min_time])/timestep)
        
        # index of the end point within trajectory
        End_Step = int(np.min([End_Points[idx],max_time])/timestep)
        
        time_slices.append(list(range(Start_Step, End_Step+1, steps_per_frame)))
        tot_calcs += len(time_slices[idx])
        
    sg.one_line_progress_meter('', 0.1, tot_calcs, 'key',
                           'Calculating Nearest Neighbours',no_titlebar=False,
                           orientation="h")
    sbidx = 0;
    
    NN_result_TimePoints = []
    locmin_prev = [0,0]
    for idx in range(0,N_SP):

        time_slice = time_slices[idx]
                
        if RDF_Req:
            from scipy.ndimage import gaussian_filter1d
            from scipy.signal import argrelextrema
            from scipy import sparse
            from scipy.sparse.linalg import spsolve
            
            # initialize rdf
            rdf = freud.density.RDF(200, 0.8)
            
            sg.one_line_progress_meter('', sbidx+0.5, tot_calcs, 'key',
                                   'Calculating Shells',no_titlebar=False,
                                   orientation="h")
            
            # compute rdf with Freud
            for ts in t.trajectory[time_slice]:
                if pbc_on:
                    box_data = to_freud_box(t.atoms.dimensions)
                else:
                    # If pbc is off, place particles in ~infinite box
                    box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
                    
                pos_data = t.atoms.positions/10
                rdf.compute(system=(box_data, pos_data), reset=False)
            
            rdf_x = rdf.bin_centers
            rdf_y = getattr(rdf, 'rdf')
            
            # Smooth rdf
            rdf_yp = gaussian_filter1d(rdf_y, 6)
            # Remove baseline
            rdf_ypp = baseline_als(rdf_yp, 1e5, 0.5, niter=10)
    
            # import matplotlib.pyplot as plt
            # fig, ax = plt.subplots(1, 1, figsize=(12, 8))
            # ax.plot(rdf_x, rdf_y, label='RDF')
            # ax.plot(rdf_x, rdf_yp, label='RDF')
            # ax.plot(rdf_x, rdf_yp - rdf_ypp, label='RDF')
    
            # get local minima
            locmin_idx = argrelextrema(rdf_yp - rdf_ypp, np.less)
            locmin = rdf_x[locmin_idx[0]]
            locmin = np.delete(locmin,0)
            
            # get inflection points
            smooth_d2 = np.gradient(np.gradient(rdf_yp))
            infls_idx = np.where(np.diff(np.sign(smooth_d2)))[0]
            infls = rdf_x[infls_idx]
            
            if len(locmin) < 2 & Second_Shell_Switch:
                if len(locmin) < 1:
                    locmin = [0, 0]
                
                if len(infls) < 4:
                    if idx == 0:
                        import warnings
                        warnings.warn("Unable to locate second shell")
                        return np.nan
                    else:
                        locmin[1] = locmin_prev[1]
                else:
                    locmin[1] = infls[3]
                    
            elif len(locmin) < 1 & (First_Shell_Switch | Second_Shell_Switch):
                if len(infls) < 2:
                    if idx == 0:
                        import warnings
                        warnings.warn("Unable to locate first shell")
                        return np.nan
                    else:
                        locmin[0] = locmin_prev[0]
                else:
                    locmin[0] = infls[1]
                    
                    if First_Shell_Switch:
                        First_Shell = locmin[1] # nm
                    
                    if Second_Shell_Switch:
                        First_Shell = locmin[0] # nm
                        Second_Shell = locmin[1] # nm
            
            if First_Shell_Switch:
                First_Shell = locmin[0] # nm
            
            if Second_Shell_Switch:
                First_Shell = locmin[0] # nm
                Second_Shell = locmin[1] # nm
    
            if First_Shell_Switch & Second_Shell_Switch:
                Max_R = Second_Shell
                
            elif First_Shell_Switch:
                Max_R = First_Shell
                
            elif Second_Shell_Switch:
                Max_R = Second_Shell
                Min_R = First_Shell
            
            locmin_prev = locmin
    
        sg.one_line_progress_meter('', sbidx+0.5, tot_calcs, 'key',
                               'Calculating Nearest Neighbours',no_titlebar=False,
                               orientation="h")
        
        query_args = dict(mode='ball', r_min=Min_R, r_max=Max_R, exclude_ii=True)
    
        NN_result_TimeWindow = np.empty_like([],dtype=int)
        for ts in t.trajectory[time_slice]: #jdx, pos_data in enumerate(current_traj.xyz):
            
            if pbc_on:
                box_data = to_freud_box(t.atoms.dimensions)
            else:
                # If pbc is off, place particles in ~infinite box
                box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
                        
            query_data = t.atoms.positions/10 # in nm
            if ref == Metal:
                query_data[t.select_atoms("name " + Halide).indices] = np.nan
            elif ref == Halide:
                query_data[t.select_atoms("name " + Metal).indices] = np.nan
        
            point_data = t.atoms.positions/10 # in nm
            if sel == Metal:
                point_data[t.select_atoms("name " + Halide).indices] = np.nan
            elif sel == Halide:
                point_data[t.select_atoms("name " + Metal).indices] = np.nan
            
            nlist = freud.locality.AABBQuery(box_data, point_data).query(query_data, query_args).toNeighborList()
            NN_list = nlist.neighbor_counts;
            
            # append to data.
            NN_result_TimeWindow = np.append(NN_result_TimeWindow,NN_list)
            
            sbidx += 1
            sg.one_line_progress_meter('', sbidx, tot_calcs, 'key',
                                   'Calculating Nearest Neighbours',no_titlebar=False,
                                   orientation="h")
            
        NN_result_TimePoints.append(NN_result_TimeWindow)
        
    return NN_result_TimePoints

def Slice_Traj(trajfile,grofile,outfile,Start_Point,End_Point,TimePerFrame,Filters,pbc_on,Viewer):
    
    import MDAnalysis as md
    import numpy as np
    import PySimpleGUI as sg
    import freud
    
    def to_freud_box(dims):
        a = dims[0]/10
        b = dims[1]/10
        c = dims[2]/10
        alpha = np.deg2rad(dims[3])
        beta = np.deg2rad(dims[4])
        gamma = np.deg2rad(dims[5])
    
        Omega = a*b*c*np.sqrt(1 - (np.cos(alpha)**2) - (np.cos(beta)**2) - (np.cos(gamma)**2) + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        cx = c*np.cos(beta)
        cy = c*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        cz = Omega/(a*b*np.sin(gamma))
    
        Transform_Matrix = [[a,             b*np.cos(gamma),      cx],
                            [0,             b*np.sin(gamma),      cy],
                            [0,             0,                    cz]]
        
        return freud.box.Box.from_matrix(Transform_Matrix)
    
    def baseline_als(y, lam, p, niter=10):
        L = len(y)
        D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
        w = np.ones(L)
        for i in range(niter):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.transpose())
            z = spsolve(Z, w*y)
            w = p * (y > z) + (1-p) * (y < z)
        return z
    
    num_filters = len(Filters)
    
    # Load the trajectory
    t = md.Universe(grofile, trajfile)
    timestep = t.trajectory.dt # ps
    steps_per_frame = int(TimePerFrame/timestep)
    max_time = t.trajectory[-1].time
    min_time = t.trajectory[0].time
    
    # index of the start point within trajectory
    Start_Step = int(np.max([Start_Point,min_time])/timestep)
    
    # index of the end point within trajectory
    End_Step = int(np.min([End_Point,max_time])/timestep)
    
    time_slice = list(range(Start_Step, End_Step+steps_per_frame, steps_per_frame))
    
    # Check if an RDF is needed
    RDF_Req = False
    First_Shell_Switch = False
    Second_Shell_Switch = False
    for Filter in Filters:
        if Filter['Type'] == 'order' or Filter['Type'] == 'neighbour':
            if Filter['Second_Shell_Sel']:
                RDF_Req = True
                First_Shell_Switch = True
                Second_Shell_Switch = True
            elif Filter['First_Shell_Sel']:
                First_Shell_Switch = True
                RDF_Req = True
    
    # If RDF is needed to get the first and second shells, perform one
    if RDF_Req:
        from scipy.ndimage import gaussian_filter1d
        from scipy.signal import argrelextrema
        from scipy import sparse
        from scipy.sparse.linalg import spsolve
        
        # initialize rdf
        num_ts = len(time_slice)
        id_ts = 0
        rdf = freud.density.RDF(200, 0.8)
        
        sg.one_line_progress_meter('', id_ts, num_ts, 'key',
                               'Calculating Shells',no_titlebar=False,
                                orientation="h")
        
        # compute rdf with Freud
        for ts in t.trajectory[time_slice]:
            if pbc_on:
                box_data = to_freud_box(t.atoms.dimensions)
            else:
                # If pbc is off, place particles in ~infinite box
                box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
                
            pos_data = t.atoms.positions/10
            rdf.compute(system=(box_data, pos_data), reset=False)
            id_ts += 1
            sg.one_line_progress_meter('', id_ts, num_ts, 'key',
                                   'Calculating Shells',no_titlebar=False,
                                   orientation="h")
        
        rdf_x = rdf.bin_centers
        rdf_y = getattr(rdf, 'rdf')
        
        # Smooth rdf
        rdf_yp = gaussian_filter1d(rdf_y, 6)
        # Remove baseline
        rdf_ypp = baseline_als(rdf_yp, 1e5, 0.5, niter=10)

        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        # ax.plot(rdf_x, rdf_y, label='RDF')
        # ax.plot(rdf_x, rdf_yp, label='RDF')
        # ax.plot(rdf_x, rdf_yp - rdf_ypp, label='RDF')

        # get local minima
        locmin_idx = argrelextrema(rdf_yp - rdf_ypp, np.less)
        locmin = rdf_x[locmin_idx[0]]
        locmin = np.delete(locmin,0)
        
        # get inflection points
        smooth_d2 = np.gradient(np.gradient(rdf_yp))
        infls_idx = np.where(np.diff(np.sign(smooth_d2)))[0]
        infls = rdf_x[infls_idx]
        
        if len(locmin) < 2 & Second_Shell_Switch:
            if len(locmin) < 1:
                locmin = [0, 0]
            
            if len(infls) < 4:
                import warnings
                warnings.warn("Unable to locate second shell")
                return np.nan
            else:
                locmin[1] = infls[3]
                
        elif len(locmin) < 1 & First_Shell_Switch:
            if len(infls) < 2:
                import warnings
                warnings.warn("Unable to locate first shell")
                return np.nan
            else:
                locmin[0] = infls[1]
        
        if First_Shell_Switch:
            First_Shell = locmin[0] # nm
        
        if Second_Shell_Switch:
            First_Shell = locmin[0] # nm
            Second_Shell = locmin[1] # nm
    
    
    # Generate global filters
    filtdx = 0;
    Global_Filter_Index = set(t.atoms.indices);
    Local_Filters = [];
    for Filter in Filters:
        if Filter['Time_Point_Sel']:
            if Filter['Type'] == 'order':

                Time_Point = Filter['Time_Point'] # Reference time in ps
                Time_Frame = int(Time_Point/timestep)
                ts = t.trajectory[Time_Frame]

                if Filter['N_Neighbours_Sel']:
                    query_args = dict(mode='nearest', num_neighbors=Filter['N_Neighbours'], exclude_ii=True)
                elif Filter['Search_Radius_Sel']:
                    query_args = dict(mode='ball', r_min=(Filter['Search_Radius_Min']/10), r_max=(Filter['Search_Radius_Max']/10), exclude_ii=True)
                elif Filter['First_Shell_Sel'] and Filter['Second_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=0, r_max=Second_Shell, exclude_ii=True)
                elif Filter['First_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=0, r_max=First_Shell, exclude_ii=True)
                elif Filter['Second_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=First_Shell, r_max=Second_Shell, exclude_ii=True)
                
                if pbc_on:
                    box_data = to_freud_box(ts.dimensions)
                else:
                    # If pbc is off, place particles in ~infinite box
                    box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
                            
                query_data = ts.positions/10 # in nm
                nlist = freud.locality.AABBQuery(box_data, query_data).query(query_data, query_args).toNeighborList()
                system = [box_data,query_data]

                if Filter['OrderParam'] == 'Ql':
                    ql = freud.order.Steinhardt(Filter['L'],wl=False,wl_normalize=True)
                elif Filter['OrderParam'] == 'Wl':
                    ql = freud.order.Steinhardt(Filter['L'],wl=True,wl_normalize=True)
                
                ql_calc = ql.compute(system, neighbors=nlist)
                
                idx_max = set(np.where(ql_calc.particle_order <= Filter['Ql_Max'])[0])
                idx_min = set(np.where(ql_calc.particle_order >= Filter['Ql_Min'])[0])
                
                # Combine the current filter with the running global filter
                Global_Filter_Index = Global_Filter_Index.intersection(idx_max.intersection(idx_min))
                
            elif Filter['Type'] == 'neighbour':
                
                Time_Point = Filter['Time_Point'] # Reference time in ps
                Time_Frame = int(Time_Point/timestep)
                ts = t.trajectory[Time_Frame]
                
                if Filter['Search_Radius_Sel']:
                    query_args = dict(mode='ball', r_min=(Filter['Search_Radius_Min']/10), r_max=(Filter['Search_Radius_Max']/10), exclude_ii=True)
                elif Filter['First_Shell_Sel'] and Filter['Second_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=0, r_max=Second_Shell, exclude_ii=True)
                elif Filter['First_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=0, r_max=First_Shell, exclude_ii=True)
                elif Filter['Second_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=First_Shell, r_max=Second_Shell, exclude_ii=True)
                
                if pbc_on:
                    box_data = to_freud_box(ts.dimensions)
                else:
                    # If pbc is off, place particles in ~infinite box
                    box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
                            
                query_data = ts.positions/10 # in nm
                nlist = freud.locality.AABBQuery(box_data, query_data).query(query_data, query_args).toNeighborList()
                
                neighbours = nlist.neighbor_counts
                
                idx_max = set(np.where(neighbours <= Filter['N_Neighbours_Max'])[0])
                idx_min = set(np.where(neighbours >= Filter['N_Neighbours_Min'])[0])
                
                # Combine the current filter with the running global filter
                Global_Filter_Index = Global_Filter_Index.intersection(idx_max.intersection(idx_min))
                
            elif Filter['Type'] == 'type':
                if Filter['Type_sel'] != 'All':
                    type_idx = set(np.where(t.atoms.names == Filter['Type_sel'])[0])
                    Global_Filter_Index = Global_Filter_Index.intersection(type_idx)
                
                
        else:
            
            if Filter['Type'] == 'order':

                if Filter['N_Neighbours_Sel']:
                    query_args = dict(mode='nearest', num_neighbors=Filter['N_Neighbours'], exclude_ii=True)
                elif Filter['Search_Radius_Sel']:
                    query_args = dict(mode='ball', r_min=(Filter['Search_Radius_Min']/10), r_max=(Filter['Search_Radius_Max']/10), exclude_ii=True)
                elif Filter['First_Shell_Sel'] and Filter['Second_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=0, r_max=Second_Shell, exclude_ii=True)
                elif Filter['First_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=0, r_max=First_Shell, exclude_ii=True)
                elif Filter['Second_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=First_Shell, r_max=Second_Shell, exclude_ii=True)
                
            elif Filter['Type'] == 'neighbour':
                
                if Filter['Search_Radius_Sel']:
                    query_args = dict(mode='ball', r_min=(Filter['Search_Radius_Min']/10), r_max=(Filter['Search_Radius_Max']/10), exclude_ii=True)
                elif Filter['First_Shell_Sel'] and Filter['Second_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=0, r_max=Second_Shell, exclude_ii=True)
                elif Filter['First_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=0, r_max=First_Shell, exclude_ii=True)
                elif Filter['Second_Shell_Sel']:
                    query_args = dict(mode='ball', r_min=First_Shell, r_max=Second_Shell, exclude_ii=True)
                    
            Filter['query_args'] = query_args # Save these so we don't have to re-create them at each time step
            Local_Filters.append(Filter)
            
        filtdx += 1; # Increment the filter counter
        sg.one_line_progress_meter('', filtdx, num_filters, 'key',
                               'Applying Global Filters',no_titlebar=False,
                               orientation="h")
            
    num_filters = len(time_slice)*len(Local_Filters)
    filtdx = 0 # filter counter
    #atoms = t.atoms
    
    #newres = t.add_Residue(segment=t.segments[0], resid=3, resname='SEL',resnum=3)
    
    #core_segment = t.add_Segment(segid='SEL')
    #core_segment.atoms
    
    
    
    with md.Writer(outfile, t.atoms.n_atoms) as W:
        for idx, ts in enumerate(t.trajectory[time_slice]):
            Local_Filter_Index = Global_Filter_Index
            
            for Filter in Local_Filters:
                if Filter['Type'] == 'order':
    
                    if pbc_on:
                        box_data = to_freud_box(ts.dimensions)
                    else:
                        # If pbc is off, place particles in ~infinite box
                        box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
                                
                    query_data = ts.positions/10 # in nm
                    nlist = freud.locality.AABBQuery(box_data, query_data).query(query_data, Filter['query_args']).toNeighborList()
                    system = [box_data,query_data]
    
                    if Filter['OrderParam'] == 'Ql':
                        ql = freud.order.Steinhardt(Filter['L'],wl=False,wl_normalize=True)
                    elif Filter['OrderParam'] == 'Wl':
                        ql = freud.order.Steinhardt(Filter['L'],wl=True,wl_normalize=True)
                    
                    ql_calc = ql.compute(system, neighbors=nlist)
                    
                    idx_max = set(np.where(ql_calc.particle_order <= Filter['Ql_Max'])[0])
                    idx_min = set(np.where(ql_calc.particle_order >= Filter['Ql_Min'])[0])
                    
                    # Combine the current filter with the running global filter and other local filters
                    Local_Filter_Index = Local_Filter_Index.intersection(idx_max.intersection(idx_min))
                    
                elif Filter['Type'] == 'neighbour':
                    
                    if pbc_on:
                        box_data = to_freud_box(ts.dimensions)
                    else:
                        # If pbc is off, place particles in ~infinite box
                        box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
                                
                    query_data = ts.positions/10 # in nm
                    nlist = freud.locality.AABBQuery(box_data, query_data).query(query_data, Filter['query_args']).toNeighborList()
                    
                    neighbours = nlist.neighbor_counts
                    
                    idx_max = set(np.where(neighbours <= Filter['N_Neighbours_Max'])[0])
                    idx_min = set(np.where(neighbours >= Filter['N_Neighbours_Min'])[0])
                    
                    # Combine the current filter with the running global filter
                    Local_Filter_Index = Local_Filter_Index.intersection(idx_max.intersection(idx_min))
            
                filtdx += 1; # Increment the filter counter
                sg.one_line_progress_meter('', filtdx, num_filters, 'key',
                       'Applying Time-Based Filters',no_titlebar=False,
                       orientation="h")
            
            #newseg = t.add_Segment(segid='SEL')
            #t.atoms[sorted(Local_Filter_Index)].segments = newseg
            
            if Viewer == 'VMD':
                t.atoms.residues.resids = 0
                t.atoms.residues.resnames = 'NOT'
                t.atoms.residues.resnums = 0
                t.atoms[sorted(Local_Filter_Index)].residues.resids = 1
                t.atoms[sorted(Local_Filter_Index)].residues.resnames = 'SEL'
                t.atoms[sorted(Local_Filter_Index)].residues.resnums = 1
                W.write(t.atoms)
            else:
                W.write(t.atoms[sorted(Local_Filter_Index)])

            
            #t.atoms[sorted(Local_Filter_Index)].residues = newres
            #sel_atoms = t.atoms[sorted(Local_Filter_Index)]
            #test.resids[sorted(Local_Filter_Index)] = 0
            #t.atoms
            
            #mask=np.full(len(t.atoms),True,dtype=bool)
            #mask[sorted(Local_Filter_Index)]=False
            
            #W.write(t.atoms)
            #W.write(t.atoms[sorted(Local_Filter_Index)])
    return

def Calculate_Liquid_Fraction(WorkDir, MLModelDir, Salt, SystemName=None, T=None,
                              T_Ref=None, RefStructure='Liquid', CheckFullTrajectory=False, 
                              SaveTrajectory=False, SaveFeatures=False, 
                              SavePredictions=False, SavePredictionsImage=False, 
                              Threshold=0.75, SlopeThreshold=0.0125, ML_TimeLength=5, 
                              ML_TimeStep=1, TimePerFrame=1, FileType='gro', 
                              Verbose=False):
    
    """
    Created on Fri Nov 20 18:54:17 2020
    
    Calculate_Liquid_Fraction(WorkDir, MLModelDir, Salt, RefStructure='Liquid', 
                              CheckFullTrajectory=True, SaveFeatures=False, 
                              SavePredictions=False, SavePredictionsImage=True,
                              Threshold=0.95, ML_TimeLength = 0, ML_TimeStep = 1)
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
                                for melting/freezing.
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
        Threshold               The minimum fraction of RefStructure required 
                                to indicate a phase change has completed
        SlopeThreshold          A second requirement: the change in the % 
                                fraction per unit time must be smaller than the 
                                absolute value of this threshold for the system 
                                to be considered at the melting point. Units of
                                [% Structure Fraction/ps]
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
    
    @author: Hayden
    """
    
    # % Import libraries and define functions
    import os
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
    import tensorflow as tf
    import MDAnalysis as md
    import numpy as np
    import freud
    import sys
    import time
    import re
    import warnings
    import logging

    if CheckFullTrajectory:
        from scipy.stats import linregress
        if SavePredictionsImage:
            from scipy.ndimage.filters import uniform_filter1d
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
        ML_Name = 'DNN-' + Trial_ID
        timeless = True
    elif Trial_ID == '0-0':
        N_neighbour_list = [18] # List of nearest neighbour numbers to calculate order parameters
        L_list = [[2,4,6,8]] # List of lists: for each neighbour number, a list of l values for the order parameters
        ML_Name = 'DNN-' + Trial_ID
        timeless = True
    else:
        N_neighbour_list = [18] # List of nearest neighbour numbers to calculate order parameters
        L_list = [[2,4,6,8]] # List of lists: for each neighbour number, a list of l values for the order parameters
        ML_Name = 'CNN-' + Trial_ID
        timeless = False
    Include_Wl = True # includes third order steinhart order parameters (Wl) as features when true
    Include_ID = True # includes atom identity (metal vs halide) as a feature when true
    
    # Load the ML model
    model_loc = os.path.join(MLModelDir,'LiX_Sturcture_Selector_CNN_Model_' + Trial_ID + '.tf')
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
    
    classes = [0,1,2,3,4,5,6]
    map_dict = {0: "Liquid",
                1: "Rocksalt",
                2: "Wurtzite",
                3: "5-5",
                4: "NiAs",
                5: "Sphalerite",
                6: "b-BeO"}
    n_classes = len(map_dict)
    
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
    outfile = ML_Name + '_' + SystemName + '.xyz'
    
    # List of index points to examine
    steps_per_init_frame = int(TimePerFrame/traj_timestep)
    Traj_starts = list(range(min_step, max_step+1, steps_per_init_frame)) # Steps
    ML_TimeLength_in_steps = int(ML_TimeLength/traj_timestep) # number of trajectory steps required to traverse the CNN time slice
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
    
    t_slice_len = int(np.ceil(max(ML_TimeLength_in_steps,1)/steps_per_frame))
    num_traj_starts = len(Traj_starts)
    
    #% Calculate features
    # Loop through the number of time batches
    t_tot = time.time()
    Ql_result = np.empty([n_features,t.atoms.n_atoms])
    Ql_result_traj = np.empty([t_slice_len,t.atoms.n_atoms,n_features])
    Ql_result_all = np.empty([num_traj_starts,t_slice_len,t.atoms.n_atoms,n_features])
    
    logging.info('Generating features...')
    for traj_start_idx, Init_step in enumerate(Traj_starts):
        t_cur = time.time()
        if Include_ID:
            namelist = t.atoms.names
            Metal_Index = np.array([x not in Metal for x in namelist]).astype(int) # 0 = metal, 1 = halide
        
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
            query_data = ts.positions/10 # in nm
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
        with open(ML_Name + '_' + SystemName + '_Features.pkl', 'wb') as f:
            pickle.dump([Ql_result_all], f)
    
    # Calculate the prediction for each atom at each time step
    logging.info('\nMaking Predictions...')
    t_tot = time.time()
    predicted_prob = np.empty([num_traj_starts,t.atoms.n_atoms,n_classes])
    predicted_classes = np.empty([num_traj_starts,t.atoms.n_atoms])
    
    if SaveTrajectory:
        # Check if processed outfle already exists and delete
        if os.path.isfile(outfile):
            os.remove(outfile)
        
        W = md.Writer(outfile, t.atoms.n_atoms)
        for idx,ts in enumerate(t.trajectory[Traj_starts]):
            t_elap = time.time() - t_tot
            logging.info("\rInferring Time Point: {:.2f} ps. Time Elapsed: {:.1f} s. ({:.2f}%, {:.2f} time points/s)".format(
                ts.frame*traj_timestep,
                t_elap,
                (idx+1)*100/num_traj_starts,
                (idx+1)/t_elap))
            
            if timeless:
                Current_data = Ql_result_all[idx]
            else:
                Current_data = Ql_result_all[idx].transpose(1, 0, 2)
            
            # Run a prediction test using the probability model on the test set
            predicted_prob[idx] = model.predict(Current_data)
        
            # Gather the predicted class for each data point and place into an array
            predicted_classes[idx] = np.argmax(predicted_prob[idx], axis=1)
            
            # Pick out the indeces of each class of interest and rename them
            index_of_interest = np.empty([len(classes),t.atoms.n_atoms],dtype=bool)
            old_atom_names = []
            for cidx, class_of_interest in enumerate(classes):
                class_sym = '_' + map_dict[class_of_interest][0]
                
                index_of_interest[cidx] = [int(x) == int(class_of_interest) for x in predicted_classes[idx]]
                
                # Save names
                old_atom_names.append(t.atoms[index_of_interest[cidx]].names)
                
                # Rename atoms of interest names
                t.atoms[index_of_interest[cidx]].names = t.atoms[index_of_interest[cidx]].names + class_sym
            
            # Combine the list of all the atoms of interest
            index_of_interest_all = np.any(index_of_interest, axis=0)
            
            
            # Picking out atoms surrounding the class of interest
            # Select out the atoms
            query_data = ts.positions[index_of_interest_all]/10 # in nm
            point_data = ts.positions/10 # in nm
                
            if pbc_on:
                box_data = to_freud_box(ts.dimensions)
            else:
                # If pbc is off, place particles in ~infinite box
                box_data = freud.box.Box(1e5,1e5,1e5,0,0,0)
            
            # Construct the neighbour filter
            query_args = dict(mode='nearest', num_neighbors=18, exclude_ii=True)
            
            # Build a neighbour list for the 18 atoms surrounding selected atoms
            nlist = freud.locality.AABBQuery(box_data, point_data).query(query_data, query_args).toNeighborList()
            
            # exclude any newly selected atoms that are already in the query data
            sel_points = np.unique(nlist.point_indices)
            query_points = np.where(index_of_interest_all)[0]
            
            # List of atom indeces for neighbour atoms that are not part of the selected set
            sel_list = np.setdiff1d(sel_points,query_points)
            
            # Modify the identity of the newly selected atoms surrounding the selected atoms
            names_old = t.atoms[sel_list].names
            t.atoms[sel_list].names = t.atoms[sel_list].names + '_Surr'
            
            # Combine the original atom index list with the nearby atom list
            index_of_interest_all = np.sort(np.concatenate((query_points,sel_list)))
                
            # Set the number of atoms at this step
            W.n_atoms = np.size(t.atoms[index_of_interest_all])
            
            # Write everything to file
            W.write(t.atoms[index_of_interest_all])
            
            # Reset Names
            for cidx, class_of_interest in enumerate(classes):
                t.atoms[index_of_interest[cidx]].names = old_atom_names[cidx]
            t.atoms[sel_list].names = names_old
        W.close()
        
        # Add in box shape at each time step
        t_tot = time.time()
        logging.info('\nModifying trajectory...')
        with open(outfile, 'r+') as f:
            
            text = f.read()
            comments = re.finditer('frame.+?Written by MDAnalysis XYZWriter.+?\n', text)
            
        
            for idx,ts in enumerate(t.trajectory[Traj_starts]):
                t_elap = time.time() - t_tot
                logging.info("\rEditing time point: {:.2f} ps. Time Elapsed: {:.1f} s. ({:.2f}%, {:.2f} time points/s)".format(
                    ts.frame*traj_timestep,
                    t_elap,
                    (idx+1)*100/num_traj_starts,
                    (idx+1)/t_elap))
                
                # Get the current box shape
                dims = np.ndarray.flatten(t.coord.triclinic_dimensions)
            
                comment_txt = 'Lattice=' + '"' + str(dims[0]) + ' '\
                        + str(dims[1]) + ' ' + str(dims[2]) + ' ' + str(dims[3])\
                        + ' ' + str(dims[4]) + ' ' + str(dims[5]) + ' ' + str(dims[6])\
                        + ' ' + str(dims[7]) + ' ' + str(dims[8]) + '"'
                comment_txt = comment_txt + ' Properties=species:S:1:pos:R:3 Time=' + str(ts.time) + '\n'
                text = re.sub('frame.+?Written by MDAnalysis XYZWriter.+?\n', r'' + comment_txt, text, count=1)
        
        with open(outfile, 'w+') as f:
            f.write(text)
            f.close()
    else:
        if timeless:
            for idx,ts in enumerate(t.trajectory[Traj_starts]):
                t_elap = time.time() - t_tot
                logging.info("\rInferring Time Point: {:.2f} ps. Time Elapsed: {:.1f} s. ({:.2f}%, {:.2f} time points/s)".format(
                    ts.frame*traj_timestep,
                    t_elap,
                    (idx+1)*100/num_traj_starts,
                    (idx+1)/t_elap))
                
                # Run a prediction test using the probability model on the test set
                predicted_prob[idx] = model.predict(Ql_result_all[idx])
                
                # Gather the predicted class for each data point and place into an array
                predicted_classes[idx] = np.argmax(predicted_prob[idx], axis=1)
            
        else:
            for idx,ts in enumerate(t.trajectory[Traj_starts]):
                t_elap = time.time() - t_tot
                logging.info("\rInferring Time Point: {:.2f} ps. Time Elapsed: {:.1f} s. ({:.2f}%, {:.2f} time points/s)".format(
                    ts.frame*traj_timestep,
                    t_elap,
                    (idx+1)*100/num_traj_starts,
                    (idx+1)/t_elap))
                
                Current_data = Ql_result_all[idx].transpose(1, 0, 2)
                
                # Run a prediction test using the probability model on the test set           
                predicted_prob[idx] = model.predict(Current_data)
                
                # Gather the predicted class for each data point and place into an array            
                predicted_classes[idx] = np.argmax(predicted_prob[idx], axis=1)
    
    # Optional: save the predictions
    if SavePredictions:
        with open(ML_Name + '_' + SystemName + '_Predictions.pkl', 'wb') as f:
            pickle.dump([predicted_prob,predicted_classes], f)
    
    
    # Finally, generate a series of line plots: partitions vs time
    y_max = 100
    labels = ["Liquid","Rocksalt","Wurtzite","5-5","NiAs","Sphalerite","$\\beta$-BeO"]
    
    n_classes = len(labels)
    # Generate the histogram along one dimension
    d = hist_laxis(predicted_classes, n_classes, [0,n_classes])
    d_norm = d/t.atoms.n_atoms
    x = np.array([i * traj_timestep + min_time for i in Traj_starts])
    
    if CheckFullTrajectory and SavePredictionsImage:
    
        # colors
        cols = sns.color_palette("muted",n_classes)
        
        # Set up the matplotlib figure
        
        fig, ax = plt.subplots()
        
        for idx,y_dat in enumerate(d_norm.T):
            yhat = uniform_filter1d(y_dat,size=int(np.ceil(len(y_dat)/30)))
            plt.plot(x,y_dat*100, label=None, color=cols[idx], alpha=0.6, linewidth=3)
            plt.plot(x,yhat*100, label=labels[idx], color=cols[idx], linewidth=2, linestyle='dashed')
        
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
        
        fig.savefig(ML_Name + '_' + SystemName + T_txt2 + '.png', dpi=300, format='png')
    
    # Determine the fraction of liquid, find first time step where liquid becomes less than or greater than the given threshold 
    ref_idx = list(map_dict.keys())[list(map_dict.values()).index(RefStructure)]
    system_melting = False
    system_freezing = False
    if ref_idx == 0: # Liquid
        liquid_fraction = d_norm[:,ref_idx]
        is_frozen = liquid_fraction < 1 - Threshold
        is_melted = liquid_fraction > Threshold
        # Check slope of last 50% of trajectory (fit linear regression)
        if CheckFullTrajectory:
            x1 = x[x >= max_time*0.5] # [ps] time points
            y1 = liquid_fraction[x >= max_time*0.5]*100 # % fraction of reference
            regmodel = linregress(x1, y1)
            slope = regmodel.slope
            if slope > SlopeThreshold:
                system_melting = True
                est_time_to_phase_change = (Threshold*100 - y1[0])/slope + x1[0]
            elif slope < -SlopeThreshold:
                system_freezing = True
                est_time_to_phase_change = (100-Threshold*100 - y1[0])/slope + x1[0]

    else:
        solid_fraction = d_norm[:,ref_idx]
        is_frozen = solid_fraction > Threshold
        is_melted = solid_fraction < 1 - Threshold
        # Check slope of last 50% of trajectory
        if CheckFullTrajectory:
            x1 = x[x >= max_time*0.5] # [ps] time points
            y1 = solid_fraction[x >= max_time*0.5]*100 # % fraction of reference
            regmodel = linregress(x1, y1)
            slope = regmodel.slope
            if slope > SlopeThreshold:
                system_freezing = True
                est_time_to_phase_change = (Threshold*100 - y1[0])/slope + x1[0]
            elif slope < -SlopeThreshold:
                system_melting = True
                est_time_to_phase_change = (100-Threshold*100 - y1[0])/slope + x1[0]

    # If system is passed the threshold into frozen or melting, use that
    system_froze = any(is_frozen)
    system_melted = any(is_melted)
    
    # Otherwise, if the system is not yet frozen or melting, check if it is freezing/melting based on slope
    if not system_froze and not system_melted:
        system_froze = system_freezing
        system_melted = system_melting
    
    # If liquid has neither fully melted nor fully frozen, report back
    if not system_froze and not system_melted:
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
    if not CheckFullTrajectory:
        time_to_phase_change = np.nan
    
    return [system_froze,system_melted,time_to_phase_change]


def Check_Trr_Status(grofile,trajfile):
    import MDAnalysis as md
    import warnings
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t = md.Universe(grofile, trajfile)
        
    dt = t.trajectory[-1].time - t.trajectory[0].time # ps
    return dt


# import numpy as np
# filter1 = { "Type": "order",
#             "OrderParam": "Ql",
#             "L": 2,
#             "Ql_Min": 0,
#             "Ql_Max": 0.3,
#             "First_Shell_Sel": True,
#             "Second_Shell_Sel": False,
#             "N_Neighbours_Sel": False,
#             "Search_Radius_Sel": False,
#             "Search_Radius_Max": 4,
#             "Search_Radius_Min": 0,
#             "N_Neighbours": 8,
#             "Time_Point_Sel": False,
#             "Time_Point": 0
# }

# filter2 = { "Type": "neighbour",
# 			"N_Neighbours_Min": 4,
# 			"N_Neighbours_Max": 7,
# 			"First_Shell_Sel": True,
# 			"Second_Shell_Sel": False,
# 			"Search_Radius_Sel": False,
# 			"Search_Radius_Max": 4,
# 			"Search_Radius_Min": 0,
# 			"Time_Point_Sel": False,
# 			"Time_Point": 10000
# }

# filter3 = { "Type": "type",
# 			"Type_sel": "All",
# 			"Time_Point_Sel": True
# }

# trajfile = 'C:\\Users\Hayden\\Documents\\Patey_Lab\\Results\\LiCl\\Liq2RN_L_JC_MMD375.00000_dCRMMrd0.210b75.000_NPT\\Liq2RN_L_JC_MMD375.00000_dCRMMrd0.210b75.000_NPT.trr'
# grofile = 'C:\\Users\\Hayden\\Documents\\Patey_Lab\\Results\\LiCl\\Liq2RN_L_JC_MMD375.00000_dCRMMrd0.210b75.000_NPT\\Liq2RN_L_JC_MMD375.00000_dCRMMrd0.210b75.000_NPT.gro'
# outfile = 'C:\\Users\\Hayden\\AppData\\Local\\Temp\\RDF_App\\tp9574c9d4_b4a3_472a_99ad_2ea78701622a\\tempTRJ.pdb'
# Start_Point = 0
# End_Point = 1000
# TimePerFrame = 100
# Filters = ()
# pbc_on = True
# Viewer='OVITO'
# Slice_Traj(trajfile,grofile,outfile,Start_Point,End_Point,TimePerFrame,Filters,pbc_on,Viewer)

# import numpy as np
# trajfile = 'C:\\Users\\Hayden\\Documents\\Patey_Lab\\Results\\LiI\\StepWurtz_W_JC_NPT\\StepWurtz_W_JC_NPT.trr'
# grofile = 'C:\\Users\\Hayden\\Documents\\Patey_Lab\\Results\\LiI\\StepWurtz_W_JC_NPT\\StepWurtz_W_JC_NPT.gro'
# Metal = 'Li'
# Halide = 'I'
# Ref_Type = 'All-All'
# Start_Points = np.array([0, 10000, 35000])
# End_Points = np.array([10, 10010, 35010])
# TimePerFrame = 1
# First_Shell_Switch = False
# Second_Shell_Switch = True
# SelRadius = False
# QlEditMin = 0
# QlEdit = 4
# pbc_on = True
# L_list = np.array([2, 4, 6, 8, 10])

# test = Calculate_NN(trajfile,grofile,Metal,Halide,Ref_Type,Start_Points,
#     End_Points,TimePerFrame,First_Shell_Switch,Second_Shell_Switch,
#     QlEditMin,QlEdit,pbc_on)


# WorkDir = 'C:/Users/Hayden/Documents/Patey_Lab/Melting_Point_Results/MP-M-3_R_JC_Model_CM3_NPT'
# MLModelDir = 'C:/Users/Hayden/Documents/Patey_Lab/Machine_Learning_Structure_Finder/LiX_Structure_Selector_Optimization'
# Salt = 'LiCl'

# Test = Calculate_Liquid_Fraction(WorkDir, MLModelDir, Salt, SystemName='MP-M-3_R_JC_Model_CM3_NPT', T_Ref=1201.45,
#                               RefStructure='Liquid', CheckFullTrajectory=True, 
#                               SaveTrajectory=False, SaveFeatures=False, 
#                               SavePredictions=False, SavePredictionsImage=True, 
#                               Threshold=0.75, ML_TimeLength = 5, ML_TimeStep = 1, 
#                               TimePerFrame=10, FileType='gro', Verbose=True,
#                               SlopeThreshold=0.0125)
