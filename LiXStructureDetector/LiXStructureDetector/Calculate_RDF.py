import MDAnalysis as md
import numpy as np
# import PySimpleGUI as sg
import freud

def Calculate_RDF(trajfile,grofile,Metal,Halide,Ref_Type,Start_Points,
                    End_Points,TimePerFrame,pbc_on):
    
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
        
    # sg.one_line_progress_meter('', 1, tot_calcs, 'key',
    #                        'Calculating RDF(s)',no_titlebar=False,
    #                        orientation="h")
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
            # sg.one_line_progress_meter('', sbidx, tot_calcs, 'key',
            #                        'Calculating RDF(s)',no_titlebar=False,
            #                        orientation="h")

        RDF_result_TimePoints.append(np.array([rdf.bin_centers,rdf.rdf]))
        CDF_result_TimePoints.append(np.array([rdf.bin_centers,rdf.n_r]))
        
    return [np.array(RDF_result_TimePoints),np.array(CDF_result_TimePoints)]