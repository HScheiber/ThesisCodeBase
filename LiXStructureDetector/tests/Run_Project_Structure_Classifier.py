import LiXStructureDetector
import numpy as np
import os
import time

# Run_Project_Structure_Classifier
Project_dir = r'/home/haydensc/project/Molten_Salts_MD';
Scratch_dir = r'/home/haydensc/scratch/Molten_Salts_MD';
ML_Models = np.array( [[1,0,0],[1,5,1],[2,0,0],[2,0,1],[2,3,1],[2,5,1],[2,11,1],[2,11,5],[2,21,5],[2,51,5],[2,51,10]] );
Salts = ['LiF','LiCl','LiBr','LiI','NaCl','CsCl'];

for ML_idx, ML_times in enumerate(ML_Models):
    
    Version = ML_times[0]
    ML_TimeLength = ML_times[1]
    ML_TimeStep = ML_times[2]
    ML_Model = 'V' + str(Version) + 'NN-' + str(ML_TimeLength) + '-' + str(ML_TimeStep)
    
    for idx,Salt in enumerate(Salts):
        
        calcs = [f for f in os.listdir(os.path.join(Project_dir,Salt)) if os.path.isdir(os.path.join(Project_dir,Salt, f))]
        
        for jdx, JobName in enumerate(calcs):
            WorkDir = os.path.join(Project_dir,Salt,JobName)
            SaveDir = os.path.join(Scratch_dir,Salt,JobName)
            
            if not os.path.isdir(SaveDir):
                os.makedirs(SaveDir)
            
            # Check if job already completed
            outfile = os.path.join(SaveDir,ML_Model + '_' + Salt + '_' + JobName + '.png')
            
            if not os.path.isfile(outfile):
                print('Starting Job: ' + ML_Model + ' ' + Salt + ' ' + JobName + '.')
                tclock = time.time()
                LiXStructureDetector.Check_Structures(WorkDir, Salt, SystemName=JobName,
                    SaveTrajectory=True, SaveFeatures=False, 
                    SavePredictions=False, SavePredictionsImage=True, 
                    ML_TimeLength=ML_TimeLength, ML_TimeStep=ML_TimeStep, TimePerFrame=50, 
                    FileType='gro', Verbose=False, StartPoint = None,
                    EndPoint = None, Version = Version, SaveDir = SaveDir,
                    InMemory = False)
                print('Job Complete: ' + ML_Model + ' ' + Salt + ' ' + JobName 
                      + '. Time Elapsed: {:.1f} s'.format(time.time() - tclock))
            else:
                print('Job ' + ML_Model + ' ' + Salt + ' ' + JobName + ' already complete, skipping.')