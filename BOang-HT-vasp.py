
__author__ = 'Heesoo Park'
__copyright__ = 'Copyright ---------------------------'
__version__ = '0.2'
__maintainer__ = 'Heesoo Park'
__email__ = 'heesoo.p@gmail.com'
__date__ = 'Jul 11, 2021'


import PerovskiteBOang as PBOa
import math
import time
import os
if __name__ == "__main__":
    f=open('atoms.inp', 'r')
    lines= f.readlines() 
    perovskite = ""
    for i in range(len(lines)):
        if not lines[i].startswith("#"):
            perovskite = lines[i].strip()
            break
    x_dim = [(0.0, 1.0), [perovskite], (0.0, math.pi/8.0)]
   #x_dim = [(0.0, 1.0), ["6.30 6.30 6.30 Methylammonium Pb I I 2 2 2 "]]
    acq = "ei"    # 'ei'(Expected improvement ); 'pi' (Probability of Improvment); and 'lcb' (Lower confidence bound)
    predictor = "GaussianProcessRegressor"  # RandomForestRegressor; ExtraTreesRegressor; GradientBoostingRegressor; and GaussianProcessRegressor

    mc_config = {"x_dim":x_dim,
                 "acq":acq,
                 "predictor":predictor,
                 "maximize":False}

    BO_num_batches = 14
    batch_size=5
    vasp_incar = {'isif': 3,                      # Set ion relaxation mode
                  'MP_x':1 , 'MP_y':1, 'MP_z':1,  # Set Monkhorst-Pack grids
                  'ediffg':0.02,                  # Set Ion relaxation tolerance
                  'FixInorg':True                 # Fix the coordnated of the inorganic network 
                 }
    start_time=time.time()
    os.system("echo 'Starting',`date '+%H:%M:%S'` >> ~/HPC/profiling.csv")
    PBOa.RunBO(mc_config, BO_num_batches, vasp_incar, batch_size, math.pi/16.0)
    os.system("echo 'Completed',`date '+%H:%M:%S'` >> ~/HPC/profiling.csv")
    end_time=time.time()
    total_time_duration=end_time-start_time
    print("TOTAL TIME TAKEN is ",total_time_duration)

