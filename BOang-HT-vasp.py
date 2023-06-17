
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
from utils import get_default_opttask_kwargs
import numpy as np

def get_dir_name(method,perovskite_):
    root=os.getcwd()
    method = method
    for root, subdirs, files in os.walk(root, topdown=False):
        {}
    compare=0
    for names in subdirs:
        if method in names:
            if int(names.split('_')[2]) > compare:
                compare=int(names.split('_')[2])
    if compare==0:
        dir_name=method+'_'+perovskite_+'_1'
        os.makedirs(dir_name)
    else:
        addition=str((compare+1))
        dir_name=method+'_'+perovskite_+'_'+ addition
        os.makedirs(dir_name)
    if perovskite_== 'FEA':
        path1=root+'/FEA/perovskites.pkl'
        path2=root+'/FEA/perovskites.csv'
        atoms=root+'/FEA/atoms.inp'
    else:
        path1=root+'/MA-EA/perovskites.pkl'
        path2=root+'/MA-EA/perovskites.csv'
        atoms=root+'/MA-EA/atoms.inp'

    path3=root+'/'+dir_name+'/'
    shutil.copy2(path1,path3)
    shutil.copy2(path2,path3)
    return dir_name,atoms
if __name__ == "__main__":
    n =len(sys.argv)
    if n > 1:
        print("Methodology Used",sys.argv[1])
    if n >2:
        print("batch_size used",sys.argv[2])
    config = get_default_opttask_kwargs()
    if config['mongodb_uri']== None:
        mongodb_uri=input("Please insert your MongoDB URI e.g. 'username:password@cluster0.abcdefg.mongodb.net'")
        config['mongodb_uri']=mongodb_uri
    else:
        mongodb_uri=config['mongodb_uri']
        print("Your MongoDB URI is",mongodb_uri)

    perovskite_=input("Type '1' for 'Methylammonium-Ethylammonium' Cation or press Enter for default: 'Methylammonium-FluorinatedEA' Cation ")
    if perovskite_ == '':
        perovskite_ = "FEA"
    else:
        perovskite_ ="MA-EA"
    initial_concentration_=input("Type initial concentration of cation 0-1 or press Enter for default = 0.5: ")
    if initial_concentration_ == '':
        initial_concentration_ = 0.5
    if str(sys.argv[1]) == 'FDI-BO' or str(sys.argv[1]) == 'TOPK-BO':
        iterations_=input("Please insert number of iterations or press Enter for default = 8: ")
        if iterations_=="":
            iterations_=8
        batch_=input("Please insert batch size or press Enter for default = 5: ")
        if batch_=="":
            batch_=8
        iterations_ = int(iterations_)
        batch_size = int(batch_)
        predictor = "GaussianProcessRegressor"
    elif str(sys.argv[1]) == 'S-BO' or str(sys.argv[1]) == 'random':
        iterations_=input("Please insert number of iterations or press Enter for default = 100: ")
        if iterations_=="":
            iterations_=100
        iterations_ = int(iterations_)
        batch_size = 1
        if sys.argv[1] == 'random':
            predictor = "random"
        else:
            predictor = "GaussianProcessRegressor"
    print("Selected method: {}, Batch_size: {}, No. of iterations: {}, Predictor: {}, Perovskite system: {}".format(str(sys.argv[1]),batch_size,iterations_,predictor,perovskite_))
   
    method= sys.argv[1]
    folder_name,atoms_=get_dir_name(method,perovskite_)
    f=open(atoms_, 'r')
    lines= f.readlines() 
    perovskite = ""
    for i in range(len(lines)):
        if not lines[i].startswith("#"):
            perovskite = lines[i].strip()
            break
    x_dim = [(0.0, 1.0), [perovskite], (0.0, math.pi/8.0)]
    acq = "ei"    # 'ei'(Expected improvement ); 'pi' (Probability of Improvment); and 'lcb' (Lower confidence bound)
     # RandomForestRegressor; ExtraTreesRegressor; GradientBoostingRegressor; and GaussianProcessRegressor
    perovskite=perovskite[:-1]
    perovskite=perovskite+str(initial_concentration_)
    print("Initial Structure",perovskite)

    mc_config = {"x_dim":x_dim,
                 "acq":acq,
                 "predictor":predictor,
                 "maximize":False}

    vasp_incar = {'isif': 3,                      # Set ion relaxation mode
                  'MP_x':1 , 'MP_y':1, 'MP_z':1,  # Set Monkhorst-Pack grids
                  'ediffg':0.02,                  # Set Ion relaxation tolerance
                  'FixInorg':True                 # Fix the coordnated of the inorganic network 
                 }
    start_time=time.time()
    root = os.getcwd()
    os.system("echo 'Starting',`date '+%H:%M:%S'` >>" +root+"/"+folder_name+"/profiling.csv")
    PBOa.RunBO(mc_config, iterations_, vasp_incar, batch_size,folder_name,mongodb_uri, math.pi/16.0)
    os.system("echo 'Completed',`date '+%H:%M:%S'` >>" +root+"/"+folder_name+"/profiling.csv")
    end_time=time.time()
    total_time_duration=end_time-start_time
    print("Total Computation time {} hrs".format(np.round(total_time_duration/3600,2)))

