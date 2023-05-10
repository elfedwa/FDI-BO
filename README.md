## Faux data Injection Bayesian Optimization
## Table of Contents


1. [Overview](#overview)
2. [Installation](#installation)
3. [Environment Setup](#env)
3. [Execution](#exe)
4. [Citation](#cite)


## Overview
<a name="overview"></a>
Bayesian optimization with Faux-data injection loop is illustrated in the figure below: 

![image info](img/bo_loop.png)

## Installation
<a name="installation"></a>
To install the required python packages run: 

`pip install -r requirements.txt`

Rocketsled is part of this code so no need to install it. 

## Environment Setup
<a name="env"></a>
Add your corresponding VASP path in `~./bashrc` using the variable `VASP_PATH` such as below:

`export VASP_PATH='lustre/software/vasp/vasp.5.4.4.pl2/bin/'`

create conda virtural environment as below:

`conda create -n fdibo_env python==3.6`

`conda config --append channels conda-forge`

`conda activate fdibo_env`
<a name="installation"></a>

copy/paste fireworks template writer files to current fireworks directory

copy paste POTCAR FOLDER -> `/cray_home/user_name/vasp_potcar/MY_PSP`

copy paste templates to `~/site-packages/fireworks/user_objects/firetasks/templates`

`gedit ~/.pmgrc.yaml`

#Add these lines:

`PMG_DEFAULT_FUNCTIONAL: PBE_54`

`PMG_VASP_PSP_DIR: /cray_home/user_name/vasp_potcar/MY_PSP/`

## Execution 
<a name="exe"></a>
The computational experiment is initialized in  `atoms.inp`. So to run computations on  Methylammonium lead halide change with supercell size of `2x2x2` and lattice parameters of `6.30x6.30x6.30` , and initial concentration of `0.275` we add the correseponding line:
<pre>
 6.30 6.30 6.30    Methylammonium   Pb I  I      2 2 2  Ethylammonium    0.275
</pre>
The computational experiments can be executed by running the following command-line, replacing (method) with either **FDI-BO**, 
                       **TOPK-BO**,
                       **S-BO** or
                       **random** 

`python BOang-HT-vasp.py (method) </dev/null> name.log 2>1 & disown -h "$!"`

                       
## Citation 
<a name="cite"></a>
Please cite this work as:

<pre>@article{faux2023
  title={Faux-Data Injection Optimization for Accelerating
   Data-Driven Discovery of Materials},
  author={Ziaullah, Abdul Wahab and Chawla, Sanjay and El-Mellouhi, Fedwa},
  journal={Integrating Materials and Manufacturing Innovation},
  publisher={Springer}
}</pre>
