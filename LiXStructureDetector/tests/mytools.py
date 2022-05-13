


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

import LiXStructureDetector
WorkDir = r'C:\Users\Hayden\Documents\Patey_Lab\Melting_Point_Studies\NaCl\Test_R_TF_NPT\T_750.0000'
MLModelDir = 'C:/Users/Hayden/Documents/Patey_Lab/Machine_Learning_Structure_Finder/LiX_Structure_Selector_Optimization'
Salt = 'NaCl'

Test = LiXStructureDetector.Check_Structures(WorkDir, Salt, SystemName='Test_R_TF_NPT',
                        SaveTrajectory=True, SaveFeatures=False, 
                        SavePredictions=False, SavePredictionsImage=True, 
                        ML_TimeLength=5, ML_TimeStep=1, TimePerFrame=10, 
                        FileType='gro', Verbose=False)



Test = LiXStructureDetector.Calculate_Liquid_Fraction(WorkDir, MLModelDir, Salt, SystemName='Test_R_TF_NPT', T_Ref=750.0000,
                              RefStructure='Liquid', CheckFullTrajectory=True, 
                              SaveTrajectory=False, SaveFeatures=False, 
                              SavePredictions=False, SavePredictionsImage=True, 
                              Threshold=0.75, ML_TimeLength = 5, ML_TimeStep = 1, 
                              TimePerFrame=10, FileType='gro', Verbose=True,
                              SlopeThreshold=0.0125)


