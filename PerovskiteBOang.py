
__author__ = 'Heesoo Park'
__copyright__ = 'Copyright ---------------------------'
__version__ = '0.2'
__maintainer__ = 'Heesoo Park'
__email__ = 'heesoo.p@gmail.com'
__date__ = 'Jul 11, 2021'

from pymatgen.core.structure import Structure
from fireworks import Firework, LaunchPad, FWorker, ScriptTask, TemplateWriterTask
from fireworks.core.firework import Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.utilities.fw_utilities import explicit_serialize
from MixCations import SupercellPoscarTask, PoscarSelTask, PoscarSelInOrgTask, PerturbPoscarTask, MixedCationTask
from MixCations import PerovskiteAnalysis, DeltaEMixfromPickle, WaitTask

from rocketsled import MissionControl, OptTask
from fireworks.utilities.fw_utilities import explicit_serialize

import random
import sys
random.seed()
import time

from fireworks import FireTaskBase, Firework, FWAction, LaunchPad, Workflow
from fireworks.core.rocket_launcher import rapidfire
from fireworks.features.multi_launcher import launch_multiprocess
from fireworks.core.firework import FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from pymatgen.io.vasp import Poscar
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice
import pymatgen as mg
import os
import numpy as np


# list of B or C elements in ABC3 perovskite
Bs = ['V', 'Nb', 'Ta', 'W', 'Bi', 'V', 'Ge', 'Sn', 'Pb']
Cs = ['O', 'S', 'Se', 'Te', 'F', 'Cl', 'Br', 'I' ]
BCs = Bs + Cs


launchpad = LaunchPad(host="mongodb+srv://user:user@cluster0.4wrvs.mongodb.net/fireworks",uri_mode=True)
opt_label = "opt_complex_parallel"
db_info = {"launchpad": launchpad, "opt_label": opt_label}

# Make a MissionControl object
mc = MissionControl(**db_info)
#   wf_creator = wf_creator
mc.reset(hard=True)

# Reset the launchpad and optimization db for this example
launchpad.reset(password=None, require_password=False, max_reset_wo_password=2000)



organic_molecules = ['methane', 'methylammonium', 'ammonium', 'hydronium', 'hydroxylammonium',\
                     'sulfonium', 'phosphonium', 'diamine', 'formamide', 'ethylammonium',\
                     'formamidinium','methylsulfonium', 'methylphosphonium', 'guanidinium', 'dimethylammonium',\
                     'fluorinatedma', 'chlroinatedma', 'brominatedma', 'fluorinatedea']

orgcations = ['Li', 'Na', 'K', 'Rb', 'Cs',\
              'Methane', 'Methylammonium', 'Ammonium', 'Hydronium', 'Hydroxylammonium',\
              'Sulfonium', 'Phosphonium', 'Diamine', 'Formamide', 'Ethylammonium',\
              'Formamidinium','Methylsulfonium', 'Methylphosphonium', 'Guanidinium', 'Dimethylammonium',\
              'FluorinatedMA', 'ChlroinatedMA', 'BrominatedMA', 'FluorinatedEA']
shortname = ['Li', 'Na', 'K', 'Rb', 'Cs',\
             'MT', 'MA', 'AM', 'HY', 'HA',\
             'SF', 'PH', 'HZ', 'FO', 'EA',\
             'FA', 'MS', 'MP', 'GA', 'DM',\
             'FMA', 'ClMA', 'BrMA', 'FEA']
SN = dict(zip(orgcations, shortname)) 



@explicit_serialize
class DeltaEMixBO(FiretaskBase):
    """
    Calculate $\Delta E _\mathrm{mix}$,
    reading the dataframe, perovskites.pkl
    """
   
    _fw_name = "MixEenergyBO Task"

#   def _load_params(self, d):

#       self.concentration = d['concentration']


    def run_task(self, fw_spec):
        x = fw_spec["_x"]
        print("inside DeltaEMIXBO x",x)
        # x = x[:-1]
        # x = x[:-1] + [random.uniform(0, 1)]

        #       if self.get("use_global_spec"):
#           self._load_params(fw_spec)
#       else:
#           self._load_params(self)

        from pathlib import Path
        import pandas as pd

        try:
            my_file = Path("../perovskites.pkl")
            df = pd.read_pickle('../perovskites.pkl')
        except:
            print('Cannot find perovksite.pkl datavase to calculate the enthlapy of cation mix')
        
       #yA2 = self.concentration
        yA2 = x[0]
        print(df) 
        returntest = df['yA2'].mean()
        recentdf = df.iloc[[-1]]
        print("Recent df is...")
        print(recentdf) 

        # Energy of pure perovskites and set their minimum as the reference
        # Comp_A1 = df[(df['yA2']==0.0) & (df['Converged']==True)]
        # print("Comp_A1",Comp_A1)
        # MinE_A1 = Comp_A1['TotEng'].min()
        # print("MinE_A1",MinE_A1)
        #
        # Comp_A2 = df[(df['yA2']==1.0) & (df['Converged']==True)]
        # print("Comp_A2",Comp_A2)
        #
        # MinE_A2 = Comp_A2['TotEng'].min()
        # print("MinE_A2",MinE_A2)
        #
        # recent_yA2 = df.at[df.index[-1],'yA2']
        # print("recent_yA2",recent_yA2)
        #
        # enthalpy_mix = df.at[df.index[-1],'TotEng']
        # print("enthalpy_mix",enthalpy_mix)
        #
        # enthalpy_ref = (MinE_A1 * (1.0 - yA2)) + (MinE_A2 * yA2)
        # print("enthalpy_ref",enthalpy_ref)
        #
        # enthalpy_mix = (enthalpy_mix - enthalpy_ref) * 1000     # in meV/FU
        # print("Computed Dleta H_Mix (meV/F.U.): ", enthalpy_mix)
        #
        # #df_Emix = df[['Comp','A1','A2','yA2','TotEng']]
        # df_Emix = df.loc[:,('Comp','A1','A2','yA2','TotEng')]
        # print("df_Emix", df_Emix)
        # df_Emix['DHmix'] = df_Emix['TotEng'] - MinE_A1 * (1.0 - df_Emix['yA2']) - MinE_A2 * (df_Emix['yA2'])
        # print("df_Emix['DHmix']", df_Emix['DHmix'])
        Comp_A1 = df[(df['yA2'] < (0.0 + 1 / df['numAllCat'])) & (df['Converged'] == True)]

        print("INSIDE MIXCAT, Comp_A1", Comp_A1)
        MinE_A1 = Comp_A1['TotEng'].min()
        print("INSIDE MIXCAT, MinE_A1", MinE_A1)

        Comp_A2 = df[(df['yA2'] > (1.0 - 1 / df['numAllCat'])) & (df['Converged'] == True)]
        print("INSIDE MIXCAT, Comp_A2", Comp_A2)

        MinE_A2 = Comp_A2['TotEng'].min()
        print("INSIDE MIXCAT, MinE_A2", MinE_A2)

        recent_yA2 = df.at[df.index[-1], 'yA2']
        print("INSIDE MIXCAT, recent_yA2", recent_yA2)

        enthalpy_comp = df.at[df.index[-1], 'TotEng']
        print("INSIDE MIXCAT, enthalpy_comp", enthalpy_comp)

        enthalpy_ref = (MinE_A1 * (1.0 - yA2)) + (MinE_A2 * yA2)
        print("INSIDE MIXCAT, enthalpy_ref", enthalpy_ref)

        enthalpy_mix = (enthalpy_comp - enthalpy_ref) * 1000  # in eV/FU
        print("INSIDE MIXCAT, enthalpy_mix", enthalpy_mix)

        print("Computed Dleta H_Mix (meV/F.U.): ", enthalpy_mix)

        # df_Emix = df[['Comp','A1','A2','yA2','TotEng']]
        df_Emix = df.loc[:, ('Comp', 'A1', 'A2', 'yA2', 'TotEng')]
        print("INSIDE MIXCAT, df_Emix", df_Emix)

        df_Emix['DHmix'] = df_Emix['TotEng'] - MinE_A1 * (1.0 - df_Emix['yA2']) - MinE_A2 * (df_Emix['yA2'])
        print("INSIDE MIXCAT, df_Emix['DHmix']", df_Emix['DHmix'])

        df_Emix['DHmix_in_meV'] = df_Emix['DHmix'] * 1000      # in meV/FU
        df_Emix.to_csv('../DHmix.csv')

        recent_DHmix = df_Emix.at[df_Emix.index[-1],'DHmix_in_meV']
        y = recent_DHmix
        print("x: ", x)
        print()
        print("y: ", y)
       
        return FWAction(update_spec={"_y": y, "_x": x})


#######################################################
##
##   Set up the VASP's input and run parameters
##   in wf_creator(x)
##
####################################################### 
def cell_with_Cs(cation_site):
    cation = cation_site.lower()
    check_molecule = [x for x in organic_molecules if x == cation]
    if check_molecule != []:
        print("The compond consists of organic molecule, " + cation_site + ". The cell optimization will start with Cs atoms instead of the molecule.")
        return 'Cs'
    else:
        return cation_site



def wf_creator(x):

    print("=================================")
    print("New batch starts, given ", x )
    print("- - - - - - - - - - - - - - - - -")

    x1 = x[0]
    BO_input = x[1]
    vc = vasp_config  # vasp_config is a global variable, set in BO-HT-vasp.py
    molalign = x[2]
    # path_=x[3]

    #BO_input = sys.argv[1]
    rfile = BO_input.split()
    lattices = rfile[0:3]
    lattices = list(map(float, lattices))
    atoms = rfile[3:7] + rfile[10:]
    supercell_repeat = rfile[7:10]
    supercell_repeat = list(map(int, supercell_repeat))
    cation_site = atoms[0]
    ## Random/BO in atoms.inp: the concetration values are given by the workflow
    if (atoms[5] == "Random") or (atoms[5] == "BO"):
        concentration = x1
        organic_molecule_to_add = atoms[0:5:4] + [concentration]
    else:
        concentration_inp = float(atoms[5])
        print("In the atoms.inp, {} is given as the concentration of A2".format(concentration_inp))
        print("The BO method will not use this concentration. But the BO-generated concentration will be used.")
       #organic_molecule_to_add = atoms[0:5:4] + [concentration]
        concentration = x1
        organic_molecule_to_add = atoms[0:5:4] + [concentration]
        print("ORGANIC MOLECULE TO ADD",organic_molecule_to_add)
    # Cs will be replaced with the organic molecule later. And Cs's coordinates are the center of the molecule.
    atoms[0] = cell_with_Cs(atoms[0])
    print("atoms[0] = cell_with_Cs(atoms[0]) ",atoms[0],cell_with_Cs(atoms[0]) )
    atoms[4] = cell_with_Cs(atoms[4])
    print("atoms[4] = cell_with_Cs(atoms[4]) ",atoms[4],cell_with_Cs(atoms[4]) )

    print('Input lattices, element, and supercell repetition: {}, {}, and {}'\
           .format(lattices, atoms, supercell_repeat))
    print('Cations to be mixed: ', organic_molecule_to_add)
    vaspinp={'lattices':lattices, 'atoms':atoms, 'cation_site':cation_site, 
                      'organic_cation':organic_molecule_to_add, 'supercell_repeat':supercell_repeat}


    lattices = vaspinp['lattices']
    atoms = vaspinp['atoms']
    cation_site = vaspinp['cation_site']
    organic_cation = vaspinp['organic_cation']
    supercell_repeat = vaspinp['supercell_repeat']
    concentration = x[0]

    pi = 3.1415926535 # the ratio of a circle's circumference to its diameter
    specie_to_remove = ['Cs']

    sc_repeat = supercell_repeat
    compound = '{}{}{}2{}_{}x{}x{}_{}{}'.format(SN[organic_cation[0]], atoms[1], atoms[2], atoms[3], 
               sc_repeat[0], sc_repeat[1], sc_repeat[2], SN[organic_cation[1]], concentration)

    print("CONCENTRATION IN POSCAR",concentration)
    print("CONCENTRATION IN POSCAR",concentration)
    print("CONCENTRATION IN POSCAR",concentration)
    print("CONCENTRATION IN POSCAR",concentration)

    poscar_context = {'system': compound,
                    'lattice_A': lattices[0], 'lattice_B': lattices[1], 'lattice_C': lattices[2],
                    'atom_A': atoms[0], 'atom_B': atoms[1], 'atom_C': atoms[2], 'atom_Q': atoms[3]}
    incar_context = {'system': 'Relax box with fixed frac coords '+compound, 'ediffg': vc['ediffg'], 'isif': vc['isif'], 'npar': 1, 'ivdw': "20", \
                     'ibrion': 2, 'ispin': 1, 'lreal': "Auto", 'ialgo': 38, 'nsw': 800}
    kpoints_context = {'method': 'Monkhorst-Pack', 'num_X': vc['MP_x'], 'num_Y': vc['MP_y'], 'num_Z': vc['MP_z'] }
    slurm_job_name = "RunByBO"
    pbs_context = {'job_name': slurm_job_name, 'job_time': "172:00:00", 'nnode': 1, 'ncpu': 28, 'cmd': 'vasp_std' }
    firetask1 = TemplateWriterTask({'context': poscar_context, 'template_file': 'poscar_cell_shape', 'output_file': 'POSCAR'})
    supercelltask = SupercellPoscarTask({'A_repeat': sc_repeat[0], 
                    'B_repeat': sc_repeat[1], 'C_repeat': sc_repeat[2], 'compound': compound })
    
    postask = PoscarSelInOrgTask()    # Selective dynamics F F F for inorganic frame
    if vc['FixInorg'] == False:
        postask = PoscarSelTask()         # Selective dynamics T T T for all the atoms
    perturbtask = PerturbPoscarTask()
    firetask2 = TemplateWriterTask({'context': incar_context, 'template_file': 'incar_relax', 'output_file': 'INCAR'})
    firetask3 = TemplateWriterTask({'context': kpoints_context, 'template_file': 'kpoints_template', 'output_file': 'KPOINTS'})
    firetask5 = TemplateWriterTask({'context': pbs_context, 'template_file': 'slurm_script', 'output_file': 'vasp.sh'})

    # Rotate the molecules: 2*pi is the maximum rotation angle. Axis [x,y,x] sets the rotataion basis. For example, [1,0,0] rotates the moleceule along the x-axis.
    orgmoltask = MixedCationTask({'remove_atom': specie_to_remove, 'organic_mol': organic_cation,
                       'angle': molalign, 'axis': [1, 1, 1],
                       'center': [0., 0., 0.]})

    #runvaspcal = "sbatch ~/HPC/run.bash  > vasp_aw.log"
    #firetask6 = ScriptTask.from_str(runvaspcal)
    #start_time=time.time()
    firetask6=ScriptTask.from_str('~/HPC/./check.job')
    #end_time=time.time()
    #process_name=
    firetask7 = ScriptTask.from_str('pwd > where.txt')
    #time_delay_firetask = WaitTask()

    firetask_analysis = PerovskiteAnalysis({'concentration':concentration})
    # firetask_analysis = PerovskiteAnalysis({'concentration':concentration})

####################################################################################
    fw1 = Firework([firetask1, supercelltask,# perturbtask,
        postask, firetask2, firetask3, firetask5, orgmoltask,
        firetask6, firetask7,#time_delay_firetask,
        firetask_analysis], spec={"_x": x},
       #firetask_analysis], spec={"_pass_job_info": True},
        name="VASP compuation")
########################################################################################

###############################################################################################
    #dEmixtask = DeltaEMixBO({'concentration':concentration})
    dEmixtask = DeltaEMixBO()
    fw2 = Firework([dEmixtask],
                  name="Compute_dEmix", spec={"_x": x})
   
    optimization = Firework([OptTask(**db_info)], name="optimization")
 
    workflow = Workflow([fw1,fw2,optimization],
                         {fw1:fw2, fw2:optimization}  )
    return(workflow)
# #End of wf_creator



def RunBO(mc_config, numbatch, VC,batch_size,molalign=0.0):
    # Configure the optimization db with MissionControl
    x_dim = mc_config['x_dim']
    acq = mc_config['acq']
    predictor = mc_config['predictor']
    maximize = mc_config['maximize']
    batch_size=batch_size
    global vasp_config   # Sets up the VASP input
    vasp_config = VC

    mc.configure(
        wf_creator=wf_creator,
        dimensions=x_dim,
        acq=acq,
        predictor=predictor,
        maximize=maximize,
        batch_size=batch_size,
        enforce_sequential=True,

    )

    # Sets high-thoughtpu VASP
    concentration = x_dim[0]
    print("concentration",concentration)
    BO_input = x_dim[1][0]
    
    # This below line sets the input for BO instead of reading atoms.inp file    
    #BO_input = sys.argv[1]
    rfile = BO_input.split()
    lattices = rfile[0:3]
    lattices = list(map(float, lattices))
    atoms = rfile[3:7] + rfile[10:]
    supercell_repeat = rfile[7:10]
    supercell_repeat = list(map(int, supercell_repeat))
    cation_site = atoms[0]

    ## Random/BO in atoms.inp: the concetration values are given by the workflow
    if (atoms[5] == "Random") or (atoms[5] == "BO"):
        concentration = random.random()
        organic_molecule_to_add = atoms[0:5:4] + [concentration]
    else:
        concentration = float(atoms[5])
        organic_molecule_to_add = atoms[0:5:4] + [concentration]

    # Cs will be replaced with the organic molecule later. And Cs's coordinates are the center of the molecule.
    atoms[0] = cell_with_Cs(atoms[0]) 
    atoms[4] = cell_with_Cs(atoms[4])   
    print('Input lattices, element, and supercell repetition: {}, {}, and {}'\
           .format(lattices, atoms, supercell_repeat))
    print('Cations to be mixed: ', organic_molecule_to_add)
    wf_creator_input={'lattices':lattices, 'atoms':atoms, 'cation_site':cation_site, 
                      'organic_cation':organic_molecule_to_add, 'supercell_repeat':supercell_repeat}

    # Run 30 workflows + optimization
    print("printing range",range(batch_size))
    if batch_size>1:
        for bs in range(batch_size):
            print("bs",bs)
            if bs != 0:
                print("bs not equal to 0",bs)
                launchpad.add_wf(wf_creator([random.random(), BO_input, random.random()]))
            else:
                print("bs equal to 0", bs)
                launchpad.add_wf(wf_creator([concentration, BO_input, random.random()]))
    else:
        launchpad.add_wf(wf_creator([concentration, BO_input, molalign]))

    print("NUMBER OF LAUNCHES THAT WIll HAPPeN IS",numbatch,batch_size)
    #rapidfire(launchpad, nlaunches=(numbatch*batch_size*3), sleep_time=0)

    launch_multiprocess(launchpad,FWorker(),nlaunches=(numbatch*batch_size*3), sleep_time=0, loglvl="CRITICAL",num_jobs=batch_size)

    # Examine and plot the optimization
    plt = mc.plot(print_pareto=True)
#   plt.show()
    plt.savefig("bayesian_optimization.pdf")
