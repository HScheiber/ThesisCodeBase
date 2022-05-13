import MDAnalysis as md
import numpy as np
# import PySimpleGUI as sg
import freud

def Calculate_NN(trajfile,grofile,Metal,Halide,Ref_Type,Start_Points,
                    End_Points,TimePerFrame,First_Shell_Switch,Second_Shell_Switch,
                    QlEditMin,QlEdit,pbc_on):
    
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
        
    # sg.one_line_progress_meter('', 0.1, tot_calcs, 'key',
    #                        'Calculating Nearest Neighbours',no_titlebar=False,
    #                        orientation="h")
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
            
            # sg.one_line_progress_meter('', sbidx+0.5, tot_calcs, 'key',
            #                        'Calculating Shells',no_titlebar=False,
            #                        orientation="h")
            
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
    
        # sg.one_line_progress_meter('', sbidx+0.5, tot_calcs, 'key',
        #                        'Calculating Nearest Neighbours',no_titlebar=False,
        #                        orientation="h")
        
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
            # sg.one_line_progress_meter('', sbidx, tot_calcs, 'key',
            #                        'Calculating Nearest Neighbours',no_titlebar=False,
            #                        orientation="h")
            
        NN_result_TimePoints.append(NN_result_TimeWindow)
        
    return NN_result_TimePoints# -*- coding: utf-8 -*-

