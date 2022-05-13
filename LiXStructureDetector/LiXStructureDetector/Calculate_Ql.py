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

import MDAnalysis as md
import numpy as np
# import PySimpleGUI as sg
import freud
import warnings


def Calculate_Ql(trajfile,grofile,Metal,Halide,Ref_Type,Start_Points,
                    End_Points,TimePerFrame,First_Shell_Switch,Second_Shell_Switch,
                    SelRadius,NearestN,QlEditMin,QlEdit,pbc_on,L_list):
    
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
        
    # sg.one_line_progress_meter('', 0.1, tot_calcs, 'key',
    #                        'Calculating Order Parameters',no_titlebar=False,
    #                        orientation="h")
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
            
            # sg.one_line_progress_meter('', sbidx+0.5, tot_calcs, 'key',
            #                        'Calculating Shells',no_titlebar=False,
            #                        orientation="h")
            
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
    
        # sg.one_line_progress_meter('', sbidx+0.5, tot_calcs, 'key',
        #                        'Calculating Order Parameters',no_titlebar=False,
        #                        orientation="h")
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
            # sg.one_line_progress_meter('', sbidx, tot_calcs, 'key',
            #                        'Calculating Order Parameters',no_titlebar=False,
            #                        orientation="h")
            
        Ql_result_TimePoints.append(np.concatenate(Ql_result_TimeWindow,axis=1))
        
    return Ql_result_TimePoints
