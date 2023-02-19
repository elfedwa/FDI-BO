
__author__ = 'Heesoo Park, Abdul Wahab Ziaullah'
__copyright__ = 'Copyright ---------------------------'
__version__ = '0.3'
__maintainer__ = 'Heesoo Park, Abdul Wahab Ziaullah'
__email__ = 'heesoo.p@gmail.com, awahab@hbku.edu.qa'
__date__ = 'Jan 15, 2023'


import PerovskiteBOang as PBOa
import math
import time
import os
import shutil
import sys

def get_dir_name(method):
    root=os.getcwd()
    method = method
    for root, subdirs, files in os.walk(root, topdown=False):
        if len(subdirs)>0:
            print(subdirs)
    compare=0
    for names in subdirs:
        if method in names:
            if int(names.split('_')[1]) > compare:
                compare=int(names.split('_')[1])
    if compare==0:
        dir_name=method+'_1'
        os.makedirs(dir_name)
    else:
        addition=str((compare+1))
        dir_name=method+ '_'+ addition
        os.makedirs(dir_name)

    path1=root+'/Pickle_file/perovskites.pkl'
    path2=root+'/Pickle_file/perovskites.csv'

    path3=root+'/'+dir_name+'/'
    shutil.copy2(path1,path3)
    shutil.copy2(path2,path3)
    return dir_name
if __name__ == "__main__":
    n =len(sys.argv)
    print("Total Arguments passed",n,sys.argv[0])
    if n > 1:
        print("Methodology Used",sys.argv[1])
    if n >2:
        print("batch_size used",sys.argv[2])

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
     # RandomForestRegressor; ExtraTreesRegressor; GradientBoostingRegressor; and GaussianProcessRegressor

    if str(sys.argv[1]) == 'FDI-BO' or str(sys.argv[1]) == 'TOPK-BO':
        BO_num_batches = 8
        batch_size = 5
        predictor = "GaussianProcessRegressor"
        print("selected batch_size and BO_num_batches and predictor", str(sys.argv[1]), batch_size, BO_num_batches, predictor)
    elif str(sys.argv[1]) == 'S-BO' or str(sys.argv[1]) == 'random':
        BO_num_batches = 100
        batch_size = 1
        if sys.argv[1] == 'random':
            predictor = "random"
        else:
            predictor = "GaussianProcessRegressor"
        print("selected batch_size and BO_num_batches and predictor", str(sys.argv[1]), batch_size, BO_num_batches, predictor)

    mc_config = {"x_dim":x_dim,
                 "acq":acq,
                 "predictor":predictor,
                 "maximize":False}
    # methods FDI-BO, S-BO, TOPK-BO
    method= sys.argv[1]
    folder_name=get_dir_name(method)

    vasp_incar = {'isif': 3,                      # Set ion relaxation mode
                  'MP_x':1 , 'MP_y':1, 'MP_z':1,  # Set Monkhorst-Pack grids
                  'ediffg':0.02,                  # Set Ion relaxation tolerance
                  'FixInorg':True                 # Fix the coordnated of the inorganic network 
                 }
    start_time=time.time()
    root = os.getcwd()
    os.system("echo 'Starting',`date '+%H:%M:%S'` >>" +root+"/"+folder_name+"/profiling.csv")
    PBOa.RunBO(mc_config, BO_num_batches, vasp_incar, batch_size,folder_name, math.pi/16.0)
    os.system("echo 'Completed',`date '+%H:%M:%S'` >>" +root+"/"+folder_name+"/profiling.csv")

    end_time=time.time()
    total_time_duration=end_time-start_time
    print("TOTAL TIME TAKEN is ",total_time_duration)

