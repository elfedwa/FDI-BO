# FDIBO
Faux data Injection Bayesian Optimization

`pip install -r requirements.txt`


This is the new code. Rocketsled is part of this code so no need to install it. Install Firework and Pymatgen as before. Configure HPC as below:

conda create -n rocketsled_env python==3.6

conda config --append channels conda-forge

conda activate rocketsled_env


#Installing right packages in following order


pip install pymatgen

pip install numpy

pip install pymongo

pip install dnspython

copy paste fireworks template writer files to current fireworks directory

copy paste POTCAR FOLDER -> /cray_home/user_name/vasp_potcar/MY_PSP

copy paste templates to 

gedit ~/.pmgrc.yaml

#Add these lines:

PMG_DEFAULT_FUNCTIONAL: PBE_54

PMG_VASP_PSP_DIR: /cray_home/user_name/vasp_potcar/MY_PSP/

#Run in login-node using this which will automatically submit job to compute node. the terminal can be closed and the job will remain active with disown -h

python BOang-HT-vasp.py (method) </dev/null> name.log 2>1 & disown -h "$!"

available methods are : FDI-BO
                       TOPK-BO
                       S-BO
                       random
