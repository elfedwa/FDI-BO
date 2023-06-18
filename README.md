## Faux data Injection Bayesian Optimization
## Table of Contents


1. [Overview](#overview)
2. [Environment Setup](#env)
3. [Installation](#installation)
4. [Execution](#exe)
5. [Citation](#cite)


## Overview
<a name="overview"></a>
Our work on Bayesian optimization with faux-data injection for target material discovery is presented in [Faux-Data Injection Optimization for Accelerating Data-Driven Discovery of Materials](https://link.springer.com/article/10.1007/s40192-023-00301-x)

It makes use of faux-data injection in Bayesian optimization loop as illustrated in the figure below: 

![image info](img/bo_loop.png)

## Environment Setup
<a name="env"></a>

This code submits SLURM jobs to run propritary DFT software, [VASP](https://www.vasp.at/). A seperate `run.bash` must be provided based on your VASP's location, HPC's configuration and other parameters. A sample `run.bash` file is provided which can be replaced based on your setup. 

This code also makes use of MongoDB to store values from Bayesian optimization, a MongoDB URI can be provided in `defaults.yaml` as `mongodb_uri:`. You can add your URI e.g `username:password@cluster0.abcdefg.mongodb.net`. If not provided in `defaults.yaml`, the user will be prompted to add in the command-line. 

Create Conda Environment as follows:

`conda create -n fdibo_env python==3.6`

`conda config --append channels conda-forge`

`conda activate fdibo_env`

## Installation
<a name="installation"></a>
Once the virtual environment is active, install the required python packages: 

`pip install -r requirements.txt`

Modified version of Rocketsled library is part of this code so no need to install it. 

Copy/paste `POTCAR` folder to `~/vasp_potcar/MY_PSP`

Copy/paste `templates` folder to `..python-version/site-packages/fireworks/user_objects/firetasks/templates`

`gedit ~/.pmgrc.yaml`

#Add these lines:

`PMG_DEFAULT_FUNCTIONAL: PBE_54`

`PMG_VASP_PSP_DIR: ~/vasp_potcar/MY_PSP/`


## Execution 

The computational experiments can be executed by running the following command in the terminal, replacing (method) with either **FDI-BO**, 
                       **TOPK-BO**,
                       **S-BO** or
                       **random** 

`python BOang-HT-vasp.py (method)`

e.g: `python BOang-HT-vasp.py FDI-BO`

To run with no-hangup:

`python BOang-HT-vasp.py (method) </dev/null> name.log 2>1 & disown -h "$!"`

                       
## Citation 
<a name="cite"></a>

Please cite this work as:

<pre>@article{ziaullah2023faux,
  title={Faux-Data Injection Optimization for Accelerating Data-Driven Discovery of Materials},
  author={Ziaullah, Abdul Wahab and Chawla, Sanjay and El-Mellouhi, Fedwa},
  journal={Integrating Materials and Manufacturing Innovation},
  pages={1--14},
  year={2023},
  publisher={Springer}
}
</pre>

## DOI
https://doi.org/10.1007/s40192-023-00301-x
