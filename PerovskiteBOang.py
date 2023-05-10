
__author__ = 'Heesoo Park , Abdul Wahab Ziaullah'
__copyright__ = 'Copyright ---------------------------'
__version__ = '2'
__maintainer__ = 'Heesoo Park'
__email__ = 'heesoo.p@gmail.com'
__date__ = 'May 30, 2023'

import matplotlib.pyplot as plt
from pymatgen.core.structure import Structure
from fireworks import Firework, LaunchPad, FWorker, ScriptTask, TemplateWriterTask
from fireworks.core.firework import Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.utilities.fw_utilities import explicit_serialize
from MixCations import SupercellPoscarTask, PoscarSelTask, PoscarSelInOrgTask, PerturbPoscarTask, MixedCationTask
from MixCations import PerovskiteAnalysis, DeltaEMixfromPickle

# from rocketsled import MissionControl, OptTask
from task import OptTask
from control import MissionControl
from fireworks.utilities.fw_utilities import explicit_serialize

import random
import sys
random.seed()


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

    def _load_params(self, d):

        self.dir_name = d['dir_name']
        self.mission_control = d['mission_control']
        print("parameter_loaded",self.dir_name)


    def run_task(self, fw_spec):

        if self.get("use_global_spec"):
            self._load_params(fw_spec)
        else:
            self._load_params(self)
        root = os.path.abspath(os.path.join(os.getcwd(), os.path.pardir))

        x = fw_spec["_x"]
        if len(x)==4:
            x = x[:-1]
        if len(x)==5:
            x = x[:-2]
        if len(x)==6:
            x = x[:-3]

        from pathlib import Path
        import pandas as pd

        try:
            my_file = Path(root + '/' + self.dir_name +"/perovskites.pkl")
            df = pd.read_pickle(root + '/' + self.dir_name +"/perovskites.pkl")
        except:
            print('Cannot find perovksite.pkl datavase to calculate the enthlapy of cation mix')
        
        yA2 = x[0]
        print(df) 
        recentdf = df.iloc[[-1]]
        print("Recent df is...")
        print(recentdf)


        Comp_A1 = df[(df['yA2Actual'] == 0.0 ) & (df['Converged']==True)]

        MinE_A1 = Comp_A1['TotEng'].min()
        Comp_A2 = df[(df['yA2Actual']== 1.0 ) & (df['Converged']==True)]
        MinE_A2 = Comp_A2['TotEng'].min()
        recent_yA2 = df.at[df.index[-1], 'yA2Actual']
        enthalpy_comp = df.at[df.index[-1], 'TotEng']
        enthalpy_ref = (MinE_A1 * (1.0 - recent_yA2)) + (MinE_A2 * recent_yA2)
        enthalpy_mix = (enthalpy_comp - enthalpy_ref) * 1000  # in eV
        df_Emix = df.loc[:, ('Comp', 'A1', 'A2','numA1','numA2', 'yA2', 'yA2Actual','angle','TotalAtoms','TotEng','TotEngperAtom')]
        df_Emix['DHmix'] = df_Emix['TotEng'] - MinE_A1 * (1.0 - df_Emix['yA2Actual']) - MinE_A2 * (df_Emix['yA2Actual'])
        df_Emix['DHmix_in_meV'] = df_Emix['DHmix'] * 1000      # in meV/
        df_Emix['DHmix_in_meV/atom'] =  (df_Emix['DHmix'] * 1000  )/ df_Emix['TotalAtoms']      # in meV/atom
        df_Emix.to_csv(root + '/' + self.dir_name + "/DHmix.csv")
        recent_DHmix = df_Emix.at[df_Emix.index[-1],'DHmix_in_meV/atom']
        plt.figure()
        plt.scatter(range(len(df_Emix['DHmix_in_meV/atom'])),df_Emix['DHmix_in_meV/atom'])
        plt.savefig(root + '/' + self.dir_name +'/liveplot.png')
        plt.close()
        y = recent_DHmix
        print("x: ", x)
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
    dir_name=global_dir_name
    db_info=global_db_info
    mc=mission_control_to_pass
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
        concentration = x1
        organic_molecule_to_add = atoms[0:5:4] + [concentration]
    # Cs will be replaced with the organic molecule later. And Cs's coordinates are the center of the molecule.
    atoms[0] = cell_with_Cs(atoms[0])
    atoms[4] = cell_with_Cs(atoms[4])
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
    root = os.getcwd()
    if "launcher" in root:
        root = os.path.abspath(os.path.join(os.getcwd(), os.path.pardir))
    firetask6=ScriptTask.from_str(root+'./check.job '+ global_dir_name)
    firetask7 = ScriptTask.from_str('pwd > where.txt')
    firetask_analysis = PerovskiteAnalysis({'concentration':concentration,'angle':molalign,'dir_name':dir_name})

####################################################################################

    fw1 = Firework([firetask1, supercelltask,# perturbtask,
        postask, firetask2, firetask3, firetask5, orgmoltask,
       firetask6, firetask7,
        firetask_analysis], spec={"_x": x},
        name="VASP compuation")
########################################################################################

###############################################################################################
    dEmixtask = DeltaEMixBO({'dir_name':dir_name,'mission_control':mc})
    fw2 = Firework([dEmixtask],
                  name="Compute_dEmix", spec={"_x": x})
   
    optimization = Firework([OptTask(**db_info)], name="optimization")
 
    workflow = Workflow([fw1,fw2,optimization],
                         {fw1:fw2, fw2:optimization}  )
    return(workflow)
# #End of wf_creator



def RunBO(mc_config, numbatch, VC,batch_size,dir_name,molalign=0.0):
    # Configure the optimization db with MissionControl
    launchpad = LaunchPad(host="mongodb+srv://user:user@cluster0.ahbrsyi.mongodb.net/" +dir_name, uri_mode=True)
    # launchpad = LaunchPad(name="rsled")

    opt_label = "opt_complex"
    global global_db_info

    global_db_info = {"launchpad": launchpad, "opt_label": opt_label}
    db_info=global_db_info

    # Make a MissionControl object
    mc = MissionControl(**db_info)
    #   wf_creator = wf_creator
    mc.reset(hard=True)
    global mission_control_to_pass
    mission_control_to_pass = mc

    # Reset the launchpad and optimization db for this example
    launchpad.reset(password=None, require_password=False, max_reset_wo_password=2000)

    x_dim = mc_config['x_dim']
    acq = mc_config['acq']
    predictor = mc_config['predictor']
    maximize = mc_config['maximize']
    global global_dir_name
    batch_size=batch_size
    global_dir_name = dir_name
    global vasp_config   # Sets up the VASP input
    vasp_config = VC

    mc.configure(
        wf_creator=wf_creator,
        dimensions=x_dim,
        acq=acq,
        predictor=predictor,
        maximize=maximize,
        batch_size=batch_size,
        dir_name=global_dir_name


    )

    # Sets high-thoughtpu VASP
    concentration = x_dim[0]
    print("concentration",concentration)
    BO_input = x_dim[1][0]
    print("BO_input",x_dim[1][0])
    # This below line sets the input for BO instead of reading atoms.inp file    
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
            launchpad.add_wf(wf_creator([concentration, BO_input, random.random(),dir_name,db_info,mc]))
    else:
        launchpad.add_wf(wf_creator([concentration, BO_input, molalign,dir_name,db_info,mc]))

    print("NUMBER OF LAUNCHES THAT WIll HAPPeN IS",numbatch,batch_size)
    print("TESTING TRUTH CONDITION",str(dir_name.split('_')[0]) , 'FDI-BO', 'FDI-BO'==str(dir_name.split('_')[0]))
    if str(dir_name.split('_')[0]) == 'FDI-BO':
        print("launching multiprocess",dir_name.split('_')[0],batch_size,numbatch)

        launch_multiprocess(launchpad,FWorker(),nlaunches=(numbatch*batch_size*3), sleep_time=0, loglvl="CRITICAL",num_jobs=batch_size)
    elif str(dir_name.split('_')[0]) == 'TOPK-BO':
        print("launching multiprocess",dir_name.split('_')[0],batch_size,numbatch)

        launch_multiprocess(launchpad,FWorker(),nlaunches=(numbatch*batch_size*3), sleep_time=0, loglvl="CRITICAL",num_jobs=batch_size)
    else:
        print("launching single",dir_name.split('_')[0],batch_size,numbatch)

        rapidfire(launchpad, nlaunches=(numbatch * batch_size * 3), sleep_time=0)

    # Examine and plot the optimization
    plt = mc.plot(print_pareto=True)
#   plt.show()
    import os
    root = os.getcwd()
    path_sv = root + '/' + dir_name + "/bayesian_optimization.pdf"
    plt.savefig(path_sv)
