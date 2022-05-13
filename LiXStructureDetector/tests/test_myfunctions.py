import sys
sys.path.append(r"C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\LiXStructureDetector")
#sys.path.append(r"/home/user/ThesisCodeBase/LiXStructureDetector")
import LiXStructureDetector



#WorkDir = r'C:\Users\Hayden\Documents\Patey_Lab\Thesis_Projects\figures\Project_4-5\NaCl_TF_Example'
#Salt='NaCl'
#SystemName= r'ExpRef_L_TF_NPT'




# LiXStructureDetector.Check_Structures(WorkDir, Salt, SystemName=SystemName,
#                         SaveTrajectory=True, SaveFeatures=False, 
#                         SavePredictions=False, SavePredictionsImage=True, 
#                         ML_TimeLength=51, ML_TimeStep=5, TimePerFrame=1, 
#                         FileType='gro', Verbose=True, StartPoint = 15300,
#                         EndPoint = 15380, Version = 2, SaveDir = None,
#                         InMemory = False)


WorkDir = r'C:\Users\Hayden\Documents\Patey_Lab\Testing'
Salt='NaCl'
SystemName = 'Prod2_R_JC_NPT'
RefChangeThreshold = 0.25

[system_froze,system_melted,time_to_phase_change] = LiXStructureDetector.Calculate_Liquid_Fraction(WorkDir, Salt, SystemName=SystemName, T=None,
                              T_Ref=None, RefStructure='Rocksalt', CheckFullTrajectory=False, 
                              SaveTrajectory=False, SaveFeatures=False, 
                              SavePredictions=False, SavePredictionsImage=True,
                              InitialRefFrac=None, RefChangeThreshold=RefChangeThreshold, 
                              SlopeThreshold=1e10, SlopeCheckBegin=0.1,
                              ML_TimeLength=21, ML_TimeStep=5, TimePerFrame=10, 
                              FileType='gro', Verbose=True, Version = 2)

