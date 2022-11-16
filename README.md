# FDIBO
Running Faux data Injection BO

#Setting up the enviornment & Copy the files to HPC's home

conda create -n rocketsled_env python==3.6.15

conda config --append channels conda-forge

conda activate rocketsled_env

scp FDIBO.zip username@raad2.qatar.tamu.edu:

#Installing right packages in following order

pip install rocketsled

pip install pymatgen

pip install dnspython

copy paste fireworks template writer

copy paste POTCAR FOLDER -> /cray_home/user_name/vasp_potcar/MY_PSP

gedit ~/.pmgrc.yaml

#Add these lines:

PMG_DEFAULT_FUNCTIONAL: PBE_54

PMG_VASP_PSP_DIR: /cray_home/user_name/vasp_potcar/MY_PSP/

#Run in login-node using this which will automatically submit job to compute node. the terminal can be closed and the job will remain active with disown -h

python BOang-HT-vasp.py </dev/null >file.log 2>&1 & disown -h "$!"

