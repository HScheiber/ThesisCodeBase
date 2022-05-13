import MDAnalysis as md
import numpy as np
# import PySimpleGUI as sg
import freud

def Slice_Traj(trajfile,grofile,outfile,Start_Point,End_Point,TimePerFrame,
               Filters,pbc_on,Viewer):
        
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
        
        # sg.one_line_progress_meter('', id_ts, num_ts, 'key',
        #                        'Calculating Shells',no_titlebar=False,
        #                         orientation="h")
        
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
            # sg.one_line_progress_meter('', id_ts, num_ts, 'key',
            #                        'Calculating Shells',no_titlebar=False,
            #                        orientation="h")
        
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
        # sg.one_line_progress_meter('', filtdx, num_filters, 'key',
        #                        'Applying Global Filters',no_titlebar=False,
        #                        orientation="h")
            
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
                # sg.one_line_progress_meter('', filtdx, num_filters, 'key',
                #        'Applying Time-Based Filters',no_titlebar=False,
                #        orientation="h")
            
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