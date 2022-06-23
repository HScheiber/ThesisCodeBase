import sys
sys.path.append(r"C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\LiXStructureDetector")
#sys.path.append(r"/home/user/ThesisCodeBase/LiXStructureDetector")
import LiXStructureDetector


# WorkDir = r'D:\Example_Nucleations\NaCl\Interface_at_MP'
# Salt='NaCl'
# SystemName = 'Interface_at_MP'
# RefChangeThreshold = 0.25

# NaCl Example
WorkDir = r'D:\Example_Nucleations\NaCl\ExampleNuc_L_JC_NPT'
Salt = 'NaCl'
SystemName = 'ExampleNuc_L_JC_NPT'
SaveTrajectory=False
SaveFeatures=True
SavePredictions=False 
SavePredictionsImage=True
ML_TimeLength=0
ML_TimeStep=0
TimePerFrame=1
FileType='gro'
Verbose=True
Temporal_Cutoff = 0
StartPoint = 6100
EndPoint = 6400
Version = 2
SaveDir = None
InMemory = False
Temporal_Cutoff = 0
Voronoi = True
Qlm_Average = False
Prob_Interfacial = 0.9
Spatial_Reassignment = False
Spatial_Interfacial = None
SaveTrajectoryAux = 2
LoadFeatures = True


# # NaCl interface Example
# WorkDir = r'D:\Example_Nucleations\NaCl\Interface_at_MP'
# Salt = 'NaCl'
# SystemName = 'Interface_at_MP'
# SaveTrajectory=True
# SaveFeatures=True
# SavePredictions=True 
# SavePredictionsImage=True, 
# ML_TimeLength=21
# ML_TimeStep=5
# TimePerFrame=5
# FileType='gro'
# Verbose=True
# Temporal_Cutoff = 0
# StartPoint = 0-Temporal_Cutoff
# EndPoint = 2000+Temporal_Cutoff
# Version = 2
# SaveDir = None
# InMemory = False


# # CsCl Example
# WorkDir = r'D:\Example_Nucleations\CsCl\ExampleNuc_L_TF_NPT'
# Salt = 'CsCl'
# SystemName = 'ExampleNuc_L_TF_NPT'
# SaveTrajectory=True
# SaveFeatures=True
# SavePredictions=False 
# SavePredictionsImage=True
# ML_TimeLength=10
# ML_TimeStep=1
# TimePerFrame=1
# FileType='gro'
# Verbose=True
# Temporal_Cutoff = 0
# StartPoint = 12200
# EndPoint = 12600
# Version = 2
# SaveDir = None
# InMemory = False
# Temporal_Cutoff = 0
# Voronoi = False
# Qlm_Average = True
# Prob_Interfacial = None
# Spatial_Reassignment = False
# Spatial_Interfacial = None


# # LiBr Example
# WorkDir = r'D:\Example_Nucleations\LiBr\ExampleNuc_L_TF_NPT'
# Salt = 'LiBr'
# SystemName = 'ExampleNuc_L_TF_NPT'
# SaveTrajectory=True
# SaveFeatures=False
# SavePredictions=True 
# SavePredictionsImage=True
# ML_TimeLength=0
# ML_TimeStep=0
# TimePerFrame=50
# FileType='gro'
# Verbose=True
# Temporal_Cutoff = 0 #
# StartPoint = 17400-Temporal_Cutoff
# EndPoint = 17800+Temporal_Cutoff
# Version = 2
# SaveDir = None
# InMemory = False
# Voronoi = False
# Qlm_Average = True
# Probability_Interfacial = None
# Spatial_Reassignment = True
# Spatial_Interfacial = 0.5

# # LiBr Interface Example
# WorkDir = r'D:\Example_Nucleations\LiBr\Interface_at_MP'
# Salt = 'LiBr'
# SystemName = 'Interface_at_MP'
# SaveTrajectory=True
# SaveFeatures=True
# SavePredictions=False 
# SavePredictionsImage=True
# ML_TimeLength=20
# ML_TimeStep=5
# TimePerFrame=5
# FileType='gro'
# Verbose=True
# Temporal_Cutoff = 0
# StartPoint = 0
# EndPoint = None
# Version = 2
# SaveDir = None
# InMemory = False
# Temporal_Cutoff = 0
# Voronoi = False
# Qlm_Average = True
# Prob_Interfacial = None
# Spatial_Reassignment = False
# Spatial_Interfacial = None


# ## LiI Example
# WorkDir = r'D:\Example_Nucleations\LiI\ExampleNuc_L_TF_NPT'
# Salt = 'LiI'
# SystemName = 'ExampleNuc_L_TF_NPT'
# SaveTrajectory=True
# SaveFeatures=True
# SavePredictions=False 
# SavePredictionsImage=True
# ML_TimeLength=10
# ML_TimeStep=1
# TimePerFrame=1
# FileType='gro'
# Verbose=True
# Temporal_Cutoff = 0
# StartPoint = 19400
# EndPoint = 20400
# Version = 2
# SaveDir = None
# InMemory = False
# Temporal_Cutoff = 0
# Voronoi = False
# Qlm_Average = True
# Prob_Interfacial = None
# Spatial_Reassignment = False
# Spatial_Interfacial = None


# WorkDir = r'C:\Users\Hayden\Documents\Patey_Lab\Testing'
# Salt = 'LiBr'
# SystemName = 'Prod2_W_TF_NPT'
# SaveTrajectory=True
# SaveFeatures=True
# SavePredictions=True 
# SavePredictionsImage=True
# ML_TimeLength=21
# ML_TimeStep=5
# TimePerFrame=5
# FileType='gro'
# Verbose=True
# Temporal_Cutoff = 0 #
# StartPoint = Temporal_Cutoff
# EndPoint = 1250-Temporal_Cutoff
# Version = 2
# SaveDir = None
# InMemory = False






# WorkDir = r'C:\Users\Hayden\Documents\Patey_Lab\Testing\T_1291.7288'
# Salt = 'NaCl'
# SystemName = 'Set60_Rep_1_R_JC_NPT'
# SaveTrajectory=False
# SaveFeatures=False
# SavePredictions=False 
# SavePredictionsImage=True
# ML_TimeLength=20
# ML_TimeStep=5
# TimePerFrame=5
# FileType='gro'
# Verbose=True
# Temporal_Cutoff = 0
# StartPoint = None
# EndPoint = None
# Version = 2
# SaveDir = None
# InMemory = False
# Temporal_Cutoff = 0
# Voronoi = False
# Qlm_Average = True
# Prob_Interfacial = None
# Spatial_Reassignment = False
# Spatial_Interfacial = None
# T = 1291.7288
# T_Ref = 1291.7288
# RefChangeThreshold = 0.25
# CheckFullTrajectory = True



# [system_froze,system_melted,time_to_phase_change,final_ref_frac] = LiXStructureDetector.Calculate_Liquid_Fraction(WorkDir, Salt, SystemName=SystemName, T=T,
#                               T_Ref=T_Ref, RefStructure='Liquid', CheckFullTrajectory=CheckFullTrajectory, 
#                               SaveTrajectory=SaveTrajectory, SaveFeatures=SaveFeatures, 
#                               SavePredictions=SavePredictions, SavePredictionsImage=SavePredictionsImage,
#                               InitialRefFrac=None, RefChangeThreshold=RefChangeThreshold, 
#                               SlopeThreshold=1e10, SlopeCheckBegin=0.1,
#                               ML_TimeLength=ML_TimeLength, ML_TimeStep=ML_TimeStep, TimePerFrame=TimePerFrame, 
#                               FileType=FileType, Verbose=Verbose, Version=Version,
#                               Temporal_Cutoff=Temporal_Cutoff, Voronoi=Voronoi, Qlm_Average=Qlm_Average,
#                               Prob_Interfacial=Prob_Interfacial, Spatial_Reassignment=Spatial_Reassignment,
#                               Spatial_Interfacial=Spatial_Interfacial)

LiXStructureDetector.Check_Structures(WorkDir, Salt, SystemName=SystemName,
                        SaveTrajectory=SaveTrajectory, SaveFeatures=SaveFeatures, 
                        SavePredictions=SavePredictions, SavePredictionsImage=SavePredictionsImage, 
                        ML_TimeLength=ML_TimeLength, ML_TimeStep=ML_TimeStep, TimePerFrame=TimePerFrame, 
                        FileType=FileType, Verbose=Verbose, StartPoint = StartPoint,
                        EndPoint=EndPoint, Version=Version, SaveDir=SaveDir,
                        InMemory=InMemory, Temporal_Cutoff=Temporal_Cutoff,
                        Voronoi=Voronoi, Qlm_Average=Qlm_Average,
                        Prob_Interfacial=Prob_Interfacial,
                        Spatial_Reassignment=Spatial_Reassignment,
                        Spatial_Interfacial=Spatial_Interfacial,
                        SaveTrajectoryAux=SaveTrajectoryAux,
                        LoadFeatures=LoadFeatures)