#!/usr/bin/python

__author__ = 'Heesoo Park'
__copyright__ = 'Copyright ---------------------------'
__version__ = '0.2'
__maintainer__ = 'Heesoo Park'
__email__ = 'heesoo.p@gmail.com'
__date__ = 'Jul 11, 2021'

import random
random.seed()
import sys
import os
import re
import numpy as np
from numpy import array as npa
from numpy.linalg import norm
from scipy.spatial import ConvexHull
import math

from fireworks import Firework, LaunchPad, FWorker, ScriptTask, TemplateWriterTask,explicit_serialize
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.core.firework import Workflow
from fireworks.core.firework import FiretaskBase, FWAction

from pymatgen.io.vasp import Poscar
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice
from pymatgen.core.composition import Composition
from pymatgen.io.vasp.outputs import Vasprun, Procar, Poscar
from pymatgen.analysis.structure_analyzer import VoronoiConnectivity
import pymatgen as mg

import pandas as pd

# list of B or C elements in ABC3 perovskite
Bs = ['V', 'Nb', 'Ta', 'W', 'Bi', 'V', 'Ge', 'Sn', 'Pb']
Cs = ['O', 'S', 'Se', 'Te', 'F', 'Cl', 'Br', 'I' ]
BCs = Bs + Cs

@explicit_serialize
class PoscarSelTask(FiretaskBase):
   """
   Set coordinations to be relaxed
   by change 'F  F  F' to 'T  T  T' in POSCAR
   """

   _fw_name = "Poscar Sel Task"

   def run_task(self, fw_spec):
       poscar = Poscar.from_file('POSCAR')
       for i in range(len(poscar.structure)):
           poscar.structure.replace(i, poscar.structure[i].specie.name, list(poscar.structure[i].frac_coords), False,
              {'selective_dynamics':[True, True, True], 'velocities':[0.0,0.0,0.0]} )
       for i in range(len(poscar.structure)):
           #Set selective dynamics by F F F for the first Metal ion
           if poscar.structure[i].specie.name in Bs : 
              poscar.structure.replace(i, poscar.structure[i].specie.name, list(poscar.structure[i].frac_coords), False,
              {'selective_dynamics':[False, False, False], 'velocities':[0.0,0.0,0.0]} )
              break


       print("Writing POSCAR")
       print(poscar) 
       poscar.write_file('POSCAR')

@explicit_serialize
class PoscarSelInOrgTask(FiretaskBase):
    """
    Set coordinations to be fixed for only Metals
    by change 'T  T  T' to 'F  F  F' in POSCAR
    """

    _fw_name = "Poscar SelINOrg Task"

    def run_task(self, fw_spec):
        poscar = Poscar.from_file('POSCAR')
        for i in range(len(poscar.structure)):
            #Set selective dynamics as F F F for metals
            if poscar.structure[i].specie.name in BCs: 
               poscar.structure.replace(i, poscar.structure[i].specie.name, list(poscar.structure[i].frac_coords), False,
               {'selective_dynamics':[False, False, False], 'velocities':[0.0,0.0,0.0]} )


        #print("Writing POSCAR")
        #print(poscar) 
        poscar.write_file('POSCAR')


@explicit_serialize
class PoscarSelDynamTask(FiretaskBase):
    """
    To fix elements for only selected elements by the name of atoms
    To set selective dynamics by giving the list the value, for example [False,False,False], or [True,True,False]
    """

    _fw_name = "Poscar SelectDynamics Task"

    def run_task(self, fw_spec):
  
        if self.get("use_global_spec"):
            self._load_params(fw_spec)
        else:
            self._load_params(self)
        poscar = Poscar.from_file('POSCAR')
        element_species = self.elements_name_to_sel
        fix_coords = np.array(self.frac_coords_to_fix)
        dynamics_val = self.seldynamics
        for i in range(len(poscar.structure)):
            #Set selective dynamics as F F F for metals
            if poscar.structure[i].specie.name in element_species: 
               poscar.structure.replace(i, poscar.structure[i].specie.name, list(poscar.structure[i].frac_coords), False,
               {'selective_dynamics': dynamics_val, 'velocities': [0.0,0.0,0.0]} )
        
        for i in range(len(poscar.structure)):
            #Set the selective dynamics as F F F if the atom is at frac_coords_to_fix.
            if poscar.structure[i].frac_coords in fix_coords:
               poscar.structure.replace(i, poscar.structure[i].specie.name, list(poscar.structure[i].frac_coords), False,
               {'selective_dynamics': dynamics_val, 'velocities': [0.0,0.0,0.0]} )

        print("Writing POSCAR")
        print(poscar) 
        poscar.write_file('POSCAR')

    def _load_params(self, d):
         self.elements_name_to_sel = d['element_names']
         self.seldynamics = d['selective_dynamics']
         self.frac_coords_to_fix = d['frac_coords_to_fix']

@explicit_serialize
class PerturbPoscarTask(FiretaskBase):
   """
    Perturb thecoordinates in POSCAR 
   """

   _fw_name = "Perturb Poscar Task"

   def run_task(self, fw_spec):
       poscar = Poscar.from_file('POSCAR')
       dist_to_perturb = 0.40  # Angstrom
       poscar.structure.perturb(dist_to_perturb)
       #print("Writing POSCAR")
       #print(poscar) 
       Poscar.write_file(poscar, 'POSCAR')

@explicit_serialize
class AddMoleculeTask(FiretaskBase):
    """
      Replace cation with organic molecule in POSCAR 
    """

    _fw_name = "Add Molecule  Task"

    #Add extention to atom name
    def potcar_dec(self, atom):
        if atom == 'Cs':
            return 'Cs_sv'
        elif atom == 'Nb':
            return 'Nb_pv'
        elif atom == 'Ta':
            return 'Ta_pv'
        else:
            return atom
    
    # From string, molecule, this builds a molecule.
    def organic_molecule(self, molecule):
        mol = molecule.lower()
        if mol == 'methane':
            #methane
            coords = [[0.000000, 0.000000, 0.000000],
                      [0.000000, 0.000000, 1.089000],
                      [1.026719, 0.000000, -0.363000],
                      [-0.513360, -0.889165, -0.363000],
                      [-0.513360, 0.889165, -0.363000]]
            methane = Molecule(["C", "H", "H", "H", "H"], coords)
            return methane
        elif mol == 'ammonium':
            #ammonium
            coords = [[-1.43417, 0.20591, -0.09007],
                      [-0.57544, 0.43296,  0.42642],
                      [-1.45470, -0.80070, -0.29555],
                      [-1.45450, 0.73536, -0.97066],
                      [-2.25247, 0.45627, 0.47884]]
            ammonium = Molecule(["N", "H", "H", "H", "H"], coords)
            return ammonium
        elif mol == 'hydroxylammonium':
            #hydroxylammonium
            coords = [[-1.48351, 0.90614, 0.02342],
                      [-2.01241,-0.39907,  0.03535],
                      [-0.59520, 0.85120, -0.38195],
                      [-2.93656,-0.29874,  0.47905],
                      [-1.44498,-1.04675,  0.60457],
                      [-2.14562,-0.77638, -0.91613]]
            hydroxylammonium = Molecule(["O", "N", "H", "H", "H", "H"], coords)
            return hydroxylammonium
        elif mol == 'methylammonium':
            #methylammonium
            coords = [[-0.8088,  0.0000,  0.0000 ],
                      [ 0.7078, -0.0001,  0.0001 ],
                      [-1.1485,  0.7527, -0.7104 ],
                      [-1.1510,  0.2389,  1.0065 ],
                      [-1.1503, -0.9911, -0.2968 ],
                      [ 1.0836,  0.9163,  0.2707 ],
                      [ 1.0830, -0.2233, -0.9292 ],
                      [ 1.0845, -0.6927,  0.6575 ]]
            methylammonium = Molecule(["C", "N", "H", "H", "H", "H", "H", "H"], coords)
            return methylammonium
        elif mol == 'sulfonium':
            #sulfonium
            coords = [[ -2.27135,  0.74704,  0.15185 ],
                      [ -0.95883,  0.56506,  0.45656 ],
                      [ -2.16158,  0.29807, -1.12672 ],
                      [ -2.13857,  2.06143, -0.16998 ]]
            sulfonium = Molecule(["S", "H", "H", "H"], coords)
            return sulfonium
        elif mol == 'phosphonium':
            #phosphonium
            coords = [[ -2.01586,  0.88447, -0.03039 ],
                      [ -0.61753,  0.96852, -0.01997 ],
                      [ -2.54336,  1.96518, -0.74898 ],
                      [ -2.49553,  0.91988,  1.28528 ],
                      [ -2.40705, -0.31575, -0.63777 ]]
            phosphonium = Molecule(["P", "H", "H", "H", "H"], coords)
            return phosphonium
        elif mol == 'hydronium':
            #hydronium
            coords = [[ -2.05164,  0.88815, -0.06855 ],
                      [ -1.07600,  0.86290, -0.17855 ],
                      [ -2.39495,  0.18383,  0.52364 ],
                      [ -2.39565,  1.78186,  0.14930 ]]
            hydronium = Molecule(["O", "H", "H", "H"], coords)
            return hydronium
        elif mol == 'diamine':
            #diamine
            coords = [[ -2.50990, -0.06161, -0.00786 ],
                      [ -1.07120, -0.24443, -0.00665 ],
                      [ -2.89971, -0.59064,  0.78315 ],
                      [ -2.88613, -0.49671, -0.86027 ],
                      [ -0.70567,  0.16908,  0.85134 ],
                      [ -0.69238,  0.26565, -0.80477 ],
                      [ -2.85770,  0.90999,  0.04505 ]]
            diamine = Molecule(["N", "N", "H", "H", "H", "H", "H"], coords)
            return diamine
        elif mol == 'formamide':
            #formamide
            coords = [[ -2.04189, -0.29242, -0.01554 ],
                      [ -0.97107,  0.85689,  0.17850 ],
                      [ -2.71301, -0.06234, -0.76208 ],
                      [ -2.57715, -0.46004,  0.84823 ],
                      [  0.16094,  0.58456,  0.02872 ],
                      [ -1.45555,  1.80566,  0.44120 ],
                      [ -1.54222, -1.15556, -0.27044 ]]
            formamide = Molecule(["N", "C", "H", "H", "O", "H", "H"], coords)
            return formamide
        elif mol == 'ethylammonium':
            #ethylammonium
            coords = [[ -3.89859,  0.72737, -0.01905 ],
                      [ -2.58479, -0.03370,  0.00150 ],
                      [ -3.99157,  1.40991,  0.83173 ],
                      [ -4.71686,  0.00558,  0.04897 ],
                      [ -4.03358,  1.28872, -0.94926 ],
                      [ -1.40003,  0.93258, -0.11587 ],
                      [ -2.48447, -0.72630, -0.83674 ],
                      [ -2.42817, -0.58068,  0.93343 ],
                      [ -1.44565,  1.46667, -0.99148 ],
                      [ -0.50028,  0.43936, -0.09701 ],
                      [ -1.40077,  1.61111,  0.65444 ]]
            ethylammonium = Molecule(["C", "C", "H", "H", "H", "N", "H", "H", "H", "H", "H"], coords)
            return ethylammonium
        elif mol == 'dimethylammonium':
            #dimethylammonium
            coords = [[-2.8172,  0.0026, -0.0445],
                      [-1.4773,  0.6810,  0.1082],
                      [-2.8315, -0.5187, -1.0016],
                      [-3.5999,  0.7612, -0.0125],
                      [-2.9401, -0.7059,  0.7750],
                      [-0.2956, -0.2582,  0.0706],
                      [-1.4637,  1.2019,  0.9919],
                      [-0.2941, -0.7709, -0.8915],
                      [-0.4023, -0.9763,  0.8833],
                      [ 0.6185,  0.3235,  0.1931],
                      [-1.3706,  1.3842, -0.6312]]
            dimethylammonium = Molecule(["C", "N", "H", "H", "H", "C", "H", "H", "H", "H", "H"], coords)
            return dimethylammonium
        elif mol == 'formamidinium':
            #formamidinium
            coords = [[  0.0000,  0.4250,  0.0000 ],
                      [ -1.1672, -0.1778,  0.0000 ],
                      [ -2.0189,  0.3697,  0.0001 ],
                      [ -1.2697, -1.1868, -0.0005 ],
                      [  1.1672, -0.1778, -0.0000 ],
                      [ -0.0002,  1.5120, -0.0001 ],
                      [  2.0189,  0.3695, -0.0001 ],
                      [  1.2694, -1.1869,  0.0005 ]]
            formamidinium = Molecule(["C", "N", "H", "H", "N", "H", "H", "H"], coords)
            return formamidinium
        elif mol == 'guanidinium':
            #guanidinium
            coords = [[ -1.62470,  0.27633, -0.19017 ],
                      [ -0.83759,  1.20972,  0.35326 ],
                      [ -1.62997, -0.96093,  0.31510 ],
                      [ -2.40441,  0.57370, -1.23390 ],
                      [ -0.80572,  2.15400, -0.00148 ],
                      [ -0.24218,  0.99668,  1.13987 ],
                      [ -2.21213, -1.68843, -0.07238 ],
                      [ -1.04961, -1.20992,  1.10223 ],
                      [ -3.00667, -0.11933, -1.65329 ],
                      [ -2.42371,  1.49886, -1.63685 ]]
            guanidinium = Molecule(["C", "N", "N", "N", "H", "H", "H", "H", "H", "H"], coords)
            return guanidinium
        elif mol == 'fluorinatedma':
            #fluorinatedma
            coords = [[ -6.8786,  1.2714,  0.6604 ],
                      [ -6.4633,  0.1620,  0.0245 ],
                      [ -4.9355,  0.2083,  0.0743 ],
                      [ -6.7797, -0.7460,  0.5435 ],
                      [ -6.7553,  0.1653, -1.0284 ],
                      [ -4.5068, -0.6013, -0.3927 ],
                      [ -4.6139,  0.2293,  1.0510 ],
                      [ -4.5934,  1.0650, -0.3806 ]]
            fluorinatedma = Molecule(["F", "C", "N", "H", "H", "H", "H", "H"], coords)
            return fluorinatedma
        elif mol == 'fluorinatedea':
            #fluorinatedea
            coords = [[ -5.8185,  0.0772, -0.1407 ],
                      [ -3.5169, -0.1178, -0.0534 ],
                      [ -4.8463, -0.8707,  0.0494 ],
                      [ -2.3285, -1.0572,  0.1023 ],
                      [ -3.4322,  0.6366,  0.7308 ],
                      [ -3.4196,  0.3620, -1.0290 ],
                      [ -2.3337, -1.5303,  1.0143 ],
                      [ -1.4415, -0.5451,  0.0293 ],
                      [ -4.9708, -1.3342,  1.0360 ],
                      [ -4.9284, -1.6451, -0.7249 ],
                      [ -2.3211, -1.7853, -0.6225 ]]
            fluorinatedea = Molecule(["F", "C", "C", "N", "H", "H", "H", "H", "H", "H", "H"], coords)
            return fluorinatedea
        elif mol == 'chlorinatedma':
            #chlorinatedma
            coords = [[ -7.0785,  1.5595,  0.7888 ],
                      [ -6.4323,  0.1250, -0.0010 ],
                      [ -4.9117,  0.1731,  0.0612 ],
                      [ -6.7530, -0.7790,  0.5167 ],
                      [ -6.7161,  0.1015, -1.0531 ],
                      [ -4.4896, -0.6501, -0.3896 ],
                      [ -4.5892,  0.2078,  1.0372 ],
                      [ -4.5560,  1.0162, -0.4083 ]]
            chlorinatedma = Molecule(["Cl", "C", "N", "H", "H", "H", "H", "H"], coords)
            return chlorinatedma
        elif mol == 'brominatedma':
            #brominatedma
            coords = [[ -7.1487,  1.6671,  0.8472 ],
                      [ -6.4192,  0.1043, -0.0157 ],
                      [ -4.9034,  0.1578,  0.0550 ],
                      [ -6.7453, -0.7951,  0.5036 ],
                      [ -6.7039,  0.0867, -1.0664 ],
                      [ -4.4744, -0.6741, -0.3738 ],
                      [ -4.5870,  0.2166,  1.0317 ],
                      [ -4.5446,  0.9907, -0.4298 ]]
            brominatedma = Molecule(["Br", "C", "N", "H", "H", "H", "H", "H"], coords)
            return brominatedma
        else:
            print("Cannot find molecule, the name of "+molecule+".")
    
    def run_task(self, fw_spec):
        if self.get("use_global_spec"):
            self._load_params(fw_spec)
        else:
            self._load_params(self)
        #Specie of cation to remove
        specie_to_remove = self.remove_atom
        organic_molecule_to_add = self.organic_mol
        mol_rot_angle = self.angle
        mol_rot_axis = self.axis
        mol_rot_center = self.cen
    
        ##-------------------------------------------------------------------------------------
        ##From this line, it starts to replace and store to build organic-inorganic perovskite.
        ##-------------------------------------------------------------------------------------
    
        #Read POSCAR file
        p = Poscar.from_file('POSCAR')
        offset = p.structure.sites[0].coords.tolist()
    
        #Store the coordinates of cations before remove
        cation_coords = []
        for i in range(len(p.structure.sites)):
            if p.structure.sites[i].specie.name == specie_to_remove[0]:
                cation_coords += [p.structure.sites[i].coords.tolist() ]
    
        #Remove cations
        p.structure.remove_species(specie_to_remove)
    
        ## Molecule to add instead of the removed cations
        org_molecule = self.organic_molecule(organic_molecule_to_add)
        org_molecule = org_molecule.get_centered_molecule()
    
        #Add molecule
        for i in range(len(cation_coords)):
             #rotate molecules
            org_molecule.rotate_sites(list(range(len(org_molecule))), mol_rot_angle*random.random(), mol_rot_axis, mol_rot_center)
             #translate molecules
            org_molecule.translate_sites(list(range(len(org_molecule))), cation_coords[i])
             # add molecule
            for j in range(len(org_molecule.sites)):
                atom_to_add = org_molecule.species[j].name
                coords_to_add = org_molecule.cart_coords[j]
                p.structure.append(atom_to_add, coords_to_add, coords_are_cartesian=True, validate_proximity=False, 
                properties={'selective_dynamics':[True, True, True], 'velocities':[0.0, 0.0, 0.0]})
            inverse = [-1, -1, -1]
            cation_coords_res = [val*cation_coords[i][k] for k, val in enumerate(inverse)]
            org_molecule.translate_sites(list(range(len(org_molecule))), cation_coords_res)
    
        print("Writing POSCAR")
        print(poscar) 
        #Write POSCAR
        p.structure.sort(key=lambda s: s.species_string)
        p.write_file('POSCAR')
    
        #Write POTCAR
        print("Writing POTCAR")
        atoms = [p.structure[i].species_string for i in range(len(p.structure))]
        atoms = list(dict.fromkeys(atoms))
        for i in range(len(atoms)):
            atoms[i] = self.potcar_dec(atoms[i])
        atoms = " ".join(atoms)
        os.system("pmg potcar -s "+atoms)


    def _load_params(self, d):

        self.remove_atom = d['remove_atom']
        self.organic_mol = d['organic_mol']
        self.angle = d['angle']
        self.axis = d['axis']
        self.cen = d['center']


@explicit_serialize
class SupercellPoscarTask(FiretaskBase):
    """
     Perturb thecoordinates in POSCAR 
    """

    _fw_name = "Supercell Poscar Task"

    def run_task(self, fw_spec):
        if self.get("use_global_spec"):
            self._load_params(fw_spec)
        else:
            self._load_params(self)
        poscar=Structure.from_file("./POSCAR")
        poscar.make_supercell([self.A_repeat, self.B_repeat, self.C_repeat])
        poscar = Poscar(poscar)
        poscar.comment= self.comp
        #print("Writing POSCAR")
        #print(poscar) 
        poscar.write_file("POSCAR")

    def _load_params(self, d):

        self.A_repeat = d['A_repeat']
        self.B_repeat = d['B_repeat']
        self.C_repeat = d['C_repeat']
        self.comp = d['compound']

@explicit_serialize
class MakeCubicCellTask(FiretaskBase):
    """
    Make the cell as a cubic
    """

    _fw_name = "Make Cell as a Cubic"

    def run_task(self, fw_spec):
        poscar=Structure.from_file("POSCAR.old")
        print(poscar) 
        vol = poscar.volume
        print(vol)
        latticeA = vol ** (1. / 3)
        print("New lattice constants of the cubic is ", latticeA )
        newlattice = [[latticeA, 0., 0.], [0., latticeA, 0.], [0., 0., latticeA]]
        newlattice = Lattice(newlattice)
        poscar.modify_lattice(newlattice)
        poscar = Poscar(poscar)

        print("Writing POSCAR")
        print(poscar) 
        poscar.write_file("POSCAR")

@explicit_serialize
class MixedCationTask(FiretaskBase):
    """
      Replace cation with organic moleculei by mixing the catios at the given fraction in POSCAR 
    """

    _fw_name = "Add Mixed-cation  Task"


    #Add extention to atom name
    def potcar_dec(self, atom):
        if atom == 'Cs':
            return 'Cs_sv'
        elif atom == 'Fr':
            return 'Fr_sv'
        elif atom == 'Rb':
            return 'Rb_pv'
        elif atom == 'K':
            return 'K_pv'
        elif atom == 'Na':
            return 'Na_pv'

        elif atom == 'Mg':
            return 'Mg_pv'
        elif atom == 'Ca':
            return 'Ca_pv'
        elif atom == 'Sr':
            return 'Sr_sv'
        elif atom == 'Ba':
            return 'Ba_sv'

        elif atom == 'Y':
            return 'Y_sv'
        elif atom == 'Lu':
            return 'Lu'

        elif atom == 'Zr':
            return 'Zr_sv'
        elif atom == 'Hf':
            return 'Hf_pv'

        elif atom == 'V':
            return 'V_pv'
        elif atom == 'Nb':
            return 'Nb_pv'
        elif atom == 'Ta':
            return 'Ta_pv'

        elif atom == 'Mo':
            return 'Mo_pv'
        elif atom == 'W':
            return 'W_sv'

        else:
            return atom
    
    # From string, molecule, this builds a molecule.
    def organic_molecule(self, molecule):
        mol = molecule.lower()
        if mol == 'methane':
            #methane
            coords = [[0.000000, 0.000000, 0.000000],
                      [0.000000, 0.000000, 1.089000],
                      [1.026719, 0.000000, -0.363000],
                      [-0.513360, -0.889165, -0.363000],
                      [-0.513360, 0.889165, -0.363000]]
            methane = Molecule(["C", "H", "H", "H", "H"], coords)
            return methane
        elif mol == 'ammonium':
            #ammonium
            coords = [[-1.43417, 0.20591, -0.09007],
                      [-0.57544, 0.43296,  0.42642],
                      [-1.45470, -0.80070, -0.29555],
                      [-1.45450, 0.73536, -0.97066],
                      [-2.25247, 0.45627, 0.47884]]
            ammonium = Molecule(["N", "H", "H", "H", "H"], coords)
            return ammonium
        elif mol == 'hydroxylammonium':
            #hydroxylammonium
            coords = [[-1.48351, 0.90614, 0.02342],
                      [-2.01241,-0.39907,  0.03535],
                      [-0.59520, 0.85120, -0.38195],
                      [-2.93656,-0.29874,  0.47905],
                      [-1.44498,-1.04675,  0.60457],
                      [-2.14562,-0.77638, -0.91613]]
            hydroxylammonium = Molecule(["O", "N", "H", "H", "H", "H"], coords)
            return hydroxylammonium
        elif mol == 'methylammonium':
            #methylammonium
            coords = [[-0.8088,  0.0000,  0.0000 ],
                      [ 0.7078, -0.0001,  0.0001 ],
                      [-1.1485,  0.7527, -0.7104 ],
                      [-1.1510,  0.2389,  1.0065 ],
                      [-1.1503, -0.9911, -0.2968 ],
                      [ 1.0836,  0.9163,  0.2707 ],
                      [ 1.0830, -0.2233, -0.9292 ],
                      [ 1.0845, -0.6927,  0.6575 ]]
            methylammonium = Molecule(["C", "N", "H", "H", "H", "H", "H", "H"], coords)
            return methylammonium
        elif mol == 'sulfonium':
            #sulfonium
            coords = [[ -2.27135,  0.74704,  0.15185 ],
                      [ -0.95883,  0.56506,  0.45656 ],
                      [ -2.16158,  0.29807, -1.12672 ],
                      [ -2.13857,  2.06143, -0.16998 ]]
            sulfonium = Molecule(["S", "H", "H", "H"], coords)
            return sulfonium
        elif mol == 'phosphonium':
            #phosphonium
            coords = [[ -2.01586,  0.88447, -0.03039 ],
                      [ -0.61753,  0.96852, -0.01997 ],
                      [ -2.54336,  1.96518, -0.74898 ],
                      [ -2.49553,  0.91988,  1.28528 ],
                      [ -2.40705, -0.31575, -0.63777 ]]
            phosphonium = Molecule(["P", "H", "H", "H", "H"], coords)
            return phosphonium
        elif mol == 'hydronium':
            #hydronium
            coords = [[ -2.05164,  0.88815, -0.06855 ],
                      [ -1.07600,  0.86290, -0.17855 ],
                      [ -2.39495,  0.18383,  0.52364 ],
                      [ -2.39565,  1.78186,  0.14930 ]]
            hydronium = Molecule(["O", "H", "H", "H"], coords)
            return hydronium
        elif mol == 'diamine':
            #diamine
            coords = [[ -2.50990, -0.06161, -0.00786 ],
                      [ -1.07120, -0.24443, -0.00665 ],
                      [ -2.89971, -0.59064,  0.78315 ],
                      [ -2.88613, -0.49671, -0.86027 ],
                      [ -0.70567,  0.16908,  0.85134 ],
                      [ -0.69238,  0.26565, -0.80477 ],
                      [ -2.85770,  0.90999,  0.04505 ]]
            diamine = Molecule(["N", "N", "H", "H", "H", "H", "H"], coords)
            return diamine
        elif mol == 'formamide':
            #formamide
            coords = [[ -2.04189, -0.29242, -0.01554 ],
                      [ -0.97107,  0.85689,  0.17850 ],
                      [ -2.71301, -0.06234, -0.76208 ],
                      [ -2.57715, -0.46004,  0.84823 ],
                      [  0.16094,  0.58456,  0.02872 ],
                      [ -1.45555,  1.80566,  0.44120 ],
                      [ -1.54222, -1.15556, -0.27044 ]]
            formamide = Molecule(["N", "C", "H", "H", "O", "H", "H"], coords)
            return formamide
        elif mol == 'ethylammonium':
            #ethylammonium
            coords = [[ -1.3034, -0.2833, -0.0001 ],
                      [ -0.0708,  0.6034,  0.0000 ],
                      [ -1.3641, -0.9028, -0.9006 ],
                      [ -2.1894,  0.3570,  0.0174 ],
                      [ -1.3468, -0.9266,  0.8848 ],
                      [  1.2051, -0.2470, -0.0000 ],
                      [ -0.0043,  1.2345,  0.8886 ],
                      [ -0.0046,  1.2346, -0.8884 ],
                      [  1.2459, -0.8501,  0.8298 ],
                      [  2.0518,  0.3328, -0.0070 ],
                      [  1.2391, -0.8597, -0.8230 ]]
            ethylammonium = Molecule(["C", "C", "H", "H", "H", "N", "H", "H", "H", "H", "H"], coords)
            return ethylammonium
        elif mol == 'dimethylammonium':
            #dimethylammonium
            coords = [[-2.8172,  0.0026, -0.0445],
                      [-1.4773,  0.6810,  0.1082],
                      [-2.8315, -0.5187, -1.0016],
                      [-3.5999,  0.7612, -0.0125],
                      [-2.9401, -0.7059,  0.7750],
                      [-0.2956, -0.2582,  0.0706],
                      [-1.4637,  1.2019,  0.9919],
                      [-0.2941, -0.7709, -0.8915],
                      [-0.4023, -0.9763,  0.8833],
                      [ 0.6185,  0.3235,  0.1931],
                      [-1.3706,  1.3842, -0.6312]]
            dimethylammonium = Molecule(["C", "N", "H", "H", "H", "C", "H", "H", "H", "H", "H"], coords)
            return dimethylammonium
        elif mol == 'formamidinium':
            #formamidinium
            coords = [[  0.0000,  0.4250,  0.0000 ],
                      [ -1.1672, -0.1778,  0.0000 ],
                      [ -2.0189,  0.3697,  0.0001 ],
                      [ -1.2697, -1.1868, -0.0005 ],
                      [  1.1672, -0.1778, -0.0000 ],
                      [ -0.0002,  1.5120, -0.0001 ],
                      [  2.0189,  0.3695, -0.0001 ],
                      [  1.2694, -1.1869,  0.0005 ]]
            formamidinium = Molecule(["C", "N", "H", "H", "N", "H", "H", "H"], coords)
            return formamidinium
        elif mol == 'guanidinium':
            #guanidinium
            coords = [[ -1.62470,  0.27633, -0.19017 ],
                      [ -0.83759,  1.20972,  0.35326 ],
                      [ -1.62997, -0.96093,  0.31510 ],
                      [ -2.40441,  0.57370, -1.23390 ],
                      [ -0.80572,  2.15400, -0.00148 ],
                      [ -0.24218,  0.99668,  1.13987 ],
                      [ -2.21213, -1.68843, -0.07238 ],
                      [ -1.04961, -1.20992,  1.10223 ],
                      [ -3.00667, -0.11933, -1.65329 ],
                      [ -2.42371,  1.49886, -1.63685 ]]
            guanidinium = Molecule(["C", "N", "N", "N", "H", "H", "H", "H", "H", "H"], coords)
            return guanidinium
        elif mol == 'fluorinatedma':
            #fluorinatedma
            coords = [[ -6.8786,  1.2714,  0.6604 ],
                      [ -6.4633,  0.1620,  0.0245 ],
                      [ -4.9355,  0.2083,  0.0743 ],
                      [ -6.7797, -0.7460,  0.5435 ],
                      [ -6.7553,  0.1653, -1.0284 ],
                      [ -4.5068, -0.6013, -0.3927 ],
                      [ -4.6139,  0.2293,  1.0510 ],
                      [ -4.5934,  1.0650, -0.3806 ]]
            fluorinatedma = Molecule(["F", "C", "N", "H", "H", "H", "H", "H"], coords)
            return fluorinatedma
        elif mol == 'fluorinatedea':
            #fluorinatedea
            coords = [[ -1.7488, -0.1497,  0.0082 ],
                      [  0.5334, -0.5164, -0.0107 ],
                      [ -0.5681,  0.5472, -0.0093 ],
                      [  1.9196,  0.1138,  0.0077 ],
                      [  0.4771, -1.1376, -0.9063 ],
                      [  0.4586, -1.1505,  0.8747 ],
                      [  2.0785,  0.7089, -0.8146 ],
                      [  2.6519, -0.6059,  0.0096 ],
                      [ -0.5230,  1.1742, -0.9091 ],
                      [ -0.5040,  1.1884,  0.8799 ],
                      [  2.0620,  0.6963,  0.8419 ]]
            fluorinatedea = Molecule(["F", "C", "C", "N", "H", "H", "H", "H", "H", "H", "H"], coords)
            return fluorinatedea
        elif mol == 'chlorinatedma':
            #chlorinatedma
            coords = [[ -7.0785,  1.5595,  0.7888 ],
                      [ -6.4323,  0.1250, -0.0010 ],
                      [ -4.9117,  0.1731,  0.0612 ],
                      [ -6.7530, -0.7790,  0.5167 ],
                      [ -6.7161,  0.1015, -1.0531 ],
                      [ -4.4896, -0.6501, -0.3896 ],
                      [ -4.5892,  0.2078,  1.0372 ],
                      [ -4.5560,  1.0162, -0.4083 ]]
            chlorinatedma = Molecule(["Cl", "C", "N", "H", "H", "H", "H", "H"], coords)
            return chlorinatedma
        elif mol == 'brominatedma':
            #brominatedma
            coords = [[ -7.1487,  1.6671,  0.8472 ],
                      [ -6.4192,  0.1043, -0.0157 ],
                      [ -4.9034,  0.1578,  0.0550 ],
                      [ -6.7453, -0.7951,  0.5036 ],
                      [ -6.7039,  0.0867, -1.0664 ],
                      [ -4.4744, -0.6741, -0.3738 ],
                      [ -4.5870,  0.2166,  1.0317 ],
                      [ -4.5446,  0.9907, -0.4298 ]]
            brominatedma = Molecule(["Br", "C", "N", "H", "H", "H", "H", "H"], coords)
            return brominatedma
        else:
            print("Cannot find molecule, the name of "+molecule+".")
            coords = [[ 0.0, 0.0, 0.0]]
            cation = Molecule([molecule], coords)
            return cation
    
    def run_task(self, fw_spec):
        from random import shuffle
        if self.get("use_global_spec"):
            self._load_params(fw_spec)
        else:
            self._load_params(self)
        #Specie of cation to remove
        specie_to_remove = self.remove_atom
        print("specie_to_remove", specie_to_remove)
        organic_molecule_to_add = self.organic_mol
        print("organic_molecule_to_add", organic_molecule_to_add)

        mol_rot_angle = self.angle
        mol_rot_axis = self.axis
        mol_rot_center = self.cen
    
        ##-------------------------------------------------------------------------------------
        ##From this line, it starts to replace and store to build organic-inorganic perovskite.
        ##-------------------------------------------------------------------------------------
    
        #Read POSCAR file
        p = Poscar.from_file('POSCAR')
        offset = p.structure.sites[0].coords.tolist()
    
        #Store the coordinates of cations before remove
        cation_coords = []
        # print("p.structure.sites",p.structure.sites)
        # print("specie_to_remove", specie_to_remove)

        for i in range(len(p.structure.sites)):
            # print("p.structure.sites i ",p.structure.sites)
            if p.structure.sites[i].specie.name == specie_to_remove[0]:
                # print("p.structure.sites[i].specie.name ",p.structure.sites[i].specie.name )
                cation_coords += [p.structure.sites[i].coords.tolist() ]
        # import sys
        # sys.exit()

        #Remove cations
        p.structure.remove_species(specie_to_remove)
        #Shuffle the list before choose the sites for each cation in mixed-pervoskites
        shuffle(cation_coords) 
    
        ## Molecule to add instead of the removed cations
        org_moleculeA = self.organic_molecule(organic_molecule_to_add[0])
        print("org_moleculeA",org_moleculeA)
        org_moleculeB = self.organic_molecule(organic_molecule_to_add[1])
        print("org_moleculeA",org_moleculeB)

        org_moleculeA = org_moleculeA.get_centered_molecule()
        org_moleculeB = org_moleculeB.get_centered_molecule()
        suborg_frac = float(organic_molecule_to_add[2])
        print("suborg_frac",suborg_frac)

    
        #Add molecule
        num_cation_coords = len(cation_coords)
        print("num_cation_coords",num_cation_coords)

        num_subcations = round(num_cation_coords * suborg_frac)
        print("num_subcations = round(num_cation_coords * suborg_frac)",num_cation_coords,suborg_frac,round(num_cation_coords * suborg_frac))

        num_maincations = num_cation_coords - num_subcations
        print("num_maincations",num_maincations)

        for i in range(num_maincations):
             # Set up the alignment 
            mol_rot_randaxis = np.array([1, 0, 0])
            org_moleculeA.rotate_sites(list(range(len(org_moleculeA))), 2 * math.pi * random.random(), mol_rot_randaxis, mol_rot_center)
            mol_rot_randaxis = np.array([0, 1, 0])
            org_moleculeA.rotate_sites(list(range(len(org_moleculeA))), (math.pi/8 + mol_rot_angle * (random.random() - 0.5)), mol_rot_randaxis, mol_rot_center)
            mol_rot_randaxis = np.array([0, 0, 1])
            org_moleculeA.rotate_sites(list(range(len(org_moleculeA))), (math.pi/8 + mol_rot_angle * (random.random() - 0.5)), mol_rot_randaxis, mol_rot_center)
             # Moleocules randomly romtated
            if mol_rot_axis[2] == 1:
                mol_rot_randaxis = np.array([0, 0, 1])
                org_moleculeA.rotate_sites(list(range(len(org_moleculeA))), math.pi/2 * random.randint(0,3), mol_rot_randaxis, mol_rot_center)
            if mol_rot_axis[1] == 1:
                mol_rot_randaxis = np.array([0, 1, 0])
                org_moleculeA.rotate_sites(list(range(len(org_moleculeA))), math.pi/2 * random.randint(0,3), mol_rot_randaxis, mol_rot_center)
            if mol_rot_axis[0] == 1:
                mol_rot_randaxis = np.array([1, 0, 0])
                org_moleculeA.rotate_sites(list(range(len(org_moleculeA))), math.pi/2 * random.randint(0,3), mol_rot_randaxis, mol_rot_center)
             #translate molecules
            org_moleculeA.translate_sites(list(range(len(org_moleculeA))), cation_coords[i])
             # add molecule
            for j in range(len(org_moleculeA.sites)):
                atom_to_add = org_moleculeA.species[j].name
                coords_to_add = org_moleculeA.cart_coords[j]
                p.structure.append(atom_to_add, coords_to_add, coords_are_cartesian=True, validate_proximity=False, 
                properties={'selective_dynamics':[True, True, True], 'velocities':[0.0, 0.0, 0.0]})
            inverse = [-1, -1, -1]
            cation_coords_res = [val*cation_coords[i][k] for k, val in enumerate(inverse)]
            org_moleculeA.translate_sites(list(range(len(org_moleculeA))), cation_coords_res)
        for i in range(num_maincations, num_maincations + num_subcations):
             # Set up the alignment 
            mol_rot_randaxis = np.array([1, 0, 0])
            org_moleculeB.rotate_sites(list(range(len(org_moleculeB))), 2 * math.pi * random.random(), mol_rot_randaxis, mol_rot_center)
            mol_rot_randaxis = np.array([0, 1, 0])
            org_moleculeB.rotate_sites(list(range(len(org_moleculeB))), (math.pi/8 + mol_rot_angle * (random.random() - 0.5)), mol_rot_randaxis, mol_rot_center)
            mol_rot_randaxis = np.array([0, 0, 1])
            org_moleculeB.rotate_sites(list(range(len(org_moleculeB))), (math.pi/8 + mol_rot_angle * (random.random() - 0.5)), mol_rot_randaxis, mol_rot_center)
             #rotate molecules
            if mol_rot_axis[2] == 1:
                mol_rot_randaxis = np.array([0, 0, 1])
                org_moleculeB.rotate_sites(list(range(len(org_moleculeB))), math.pi/2 * random.randint(0,3), mol_rot_randaxis, mol_rot_center)
            if mol_rot_axis[1] == 1:
                mol_rot_randaxis = np.array([0, 1, 0])
                org_moleculeB.rotate_sites(list(range(len(org_moleculeB))), math.pi/2 * random.randint(0,3), mol_rot_randaxis, mol_rot_center)
            if mol_rot_axis[0] == 1:
                mol_rot_randaxis = np.array([1, 0, 0])
                org_moleculeB.rotate_sites(list(range(len(org_moleculeB))), math.pi/2 * random.randint(0,3), mol_rot_randaxis, mol_rot_center)
             #translate molecules
            org_moleculeB.translate_sites(list(range(len(org_moleculeB))), cation_coords[i])
             # add molecule
            for j in range(len(org_moleculeB.sites)):
                atom_to_add = org_moleculeB.species[j].name
                coords_to_add = org_moleculeB.cart_coords[j]
                p.structure.append(atom_to_add, coords_to_add, coords_are_cartesian=True, validate_proximity=False, 
                properties={'selective_dynamics':[True, True, True], 'velocities':[0.0, 0.0, 0.0]})
            inverse = [-1, -1, -1]
            cation_coords_res = [val*cation_coords[i][k] for k, val in enumerate(inverse)]
            org_moleculeB.translate_sites(list(range(len(org_moleculeB))), cation_coords_res)
    
        p.structure.sort(key=lambda s: s.species_string)
        #Write POTCAR
        print("Writing POTCAR")
        atoms = [p.structure[i].species_string for i in range(len(p.structure))]
        atoms = list(dict.fromkeys(atoms))
        for i in range(len(atoms)):
            atoms[i] = self.potcar_dec(atoms[i])
        atoms = " ".join(atoms)
        os.system("pmg potcar -s "+atoms)

        #Write POSCAR
        print("Writing POSCAR")
        print(p) 
        p.write_file('POSCAR')
    


    def _load_params(self, d):

        self.remove_atom = d['remove_atom']
        self.organic_mol = d['organic_mol']
        self.angle = d['angle']
        self.axis = d['axis']
        self.cen = d['center']
 


 
#Shannon's radii
radii={"Cs":1.88, "Sr":1.44, "Rb":1.75, "Ba":1.42, "Na":1.18, "Li": 0.92,
       "H3O":1.0, "NH4":1.5, "CH3NH3":2.2,
       "V":0.54, "Nb":0.64, "Ta":0.64, "Bi":0.76, "Ti":0.605, 
       "P":0.38, "Sb":0.6, "Ge":0.73, "Sn": 0.69, "Pb":1.19,
       "O":1.35, "S":1.84, "Se":1.98, "Te":2.21,
       "F":1.33, "Cl":1.81, "Br":1.96, "I":2.20}

    # Goldshmidt's tolerance factor with Shannon radii
def tolerance_factor(A, B, C):
    tf = ( radii[A] + radii[C] ) / ( math.sqrt(2) * (radii[B] + radii[C]) )
    return tf


class OrganicMolecule(object):

    def __init__(self, molecule):
        self.name = molecule.lower()
        self.energy, self.num = self.num_and_energy()
        self.energy = self.num * self.energy

    def num_and_energy(self):
# Formation energy / atom is defined.
        if self.name == 'ammonium':
            return (-0.01347 , 5)
        elif self.name == 'hydroxylammonium':
            return (-0.06424 , 6)
        elif self.name == 'methylammonium':
            return (-0.06653 , 8)
        elif self.name == 'hydronium':
            return (-0.23758 , 4)
        elif self.name == 'phosphonium':
            return (0.30154 , 5)
        elif self.name == 'sulfonium':
            return (0.21744 , 4)
        elif self.name == 'diamine':
            return (0.08418 , 7)
        elif self.name == 'formamide':
            return (-0.28895 , 7)
        elif self.name == 'ethylammonium':
            return (-0.11705 , 11)
        elif self.name == 'methylenediamine':
            return (0.54132 , 8)
        elif self.name == 'formamidinium':
            return (0.54132 , 8)
        elif self.name == 'methylsulfonium':
            return (0.00789 , 7)
        elif self.name == 'methylphosphonium':
            return (0.03763 , 8)
        elif self.name == 'guanidinium':
            return (-0.17819 , 10)
        elif self.name == 'dimethylammonium':
            return (-0.0820 , 11)
        else:
            return (0.0, 0.0)

class TheCompCell(object):

    def __init__(self,theComp):
        theComp = theComp.replace("DM", "Dma")
        theComp = theComp.replace("FMA", "Fma")
        theComp = theComp.replace("ClMA", "Clma")
        theComp = theComp.replace("BrMA", "Brma")
        theComp = theComp.replace("FEA", "Fea")
        theComp = theComp.replace("HY", "Hy")
        theComp = theComp.replace("AM", "Am")
        theComp = theComp.replace("SF", "Sf")
        theComp = theComp.replace("PH", "Ph")
        theComp = theComp.replace("HA", "Ha")
        theComp = theComp.replace("MA", "Ma")
        theComp = theComp.replace("HZ", "Hz")
        theComp = theComp.replace("FA", "Fa")
        theComp = theComp.replace("FO", "Fo")
        theComp = theComp.replace("EA", "Ea")
        theComp = theComp.replace("GA", "Gua")
        comp_cell = theComp.split("_")
        comp = comp_cell[0]
        print(comp)
        self.comp = comp
        self.cell = comp_cell[1]
        self.numcell = [int(i) for i in self.cell.split("x")]
        mainAcat = re.findall(r'([A-Z][^A-Z]*)', comp_cell[0])
        self.mainA = mainAcat[0]
        subAcat = re.match(r'([a-zA-Z]+)([0-9.]+)', comp_cell[2])
        subAcat = subAcat.groups()
        self.subA = subAcat[0]
        self.subAfrac = float(subAcat[1])

def ImagejCart(a, b, latt, cutoff):
    x = np.array([0., 0., 0. ])
    y = np.array([0., 0., 0. ])
    ab = np.subtract(a,b) 
    for i in range(3):
        ac = latt[i]
        cosine_angle = np.dot(ab, ac) / (np.linalg.norm(ab) * np.linalg.norm(ac))
        if abs(cosine_angle) > 0.8:
            x = b+latt[i]
            y = b-latt[i]
    if np.linalg.norm(x-a) <= np.linalg.norm(y-a):
#       print("distance", np.linalg.norm(x-a))
        return [x, np.linalg.norm(x-a)]
    elif np.linalg.norm(y-a) < np.linalg.norm(x-a):
#       print("distance", np.linalg.norm(y-a))
        return [y, np.linalg.norm(y-a)]

def AvgCoordNum(b_atom, c_atom, contcar, v_cutoff):
    vcf = VoronoiConnectivity(contcar, cutoff=v_cutoff)
    b_sites_idx = []
    bc_bonds = []
    con_sites = contcar.sites
    for i in range(len(con_sites)):
        if con_sites[i].specie.symbol == b_atom:
            b_sites_idx.append(i)

    for M_idx in b_sites_idx:
        v_coordinated_site = vcf.get_connections()
###     print(v_coordinated_site)
        cs_carts = []
        for b_i in range(len(v_coordinated_site)):
            element1 = con_sites[v_coordinated_site[b_i][0]].specie.symbol
            element2 = con_sites[v_coordinated_site[b_i][1]].specie.symbol
            element1_idx = v_coordinated_site[b_i][0]
            element2_idx = v_coordinated_site[b_i][1]
            bond_pair = [element1, element2]
            if (element1 in b_atom) and (element2 in c_atom) and (element1_idx == M_idx):
                i_cart = con_sites[v_coordinated_site[b_i][0]].coords
                j_cart = con_sites[v_coordinated_site[b_i][1]].coords
                dist = np.linalg.norm(i_cart - j_cart)
                if dist < v_cutoff:
                    bc_bonds.append([element1,element2,dist])
                    cs_carts.append([j_cart, element1, element2])
            elif (element2 in b_atom) and (element1 in c_atom) and element2_idx == M_idx:
                i_cart = con_sites[v_coordinated_site[b_i][1]].coords
                j_cart = con_sites[v_coordinated_site[b_i][0]].coords
                j_cart = ImagejCart(i_cart, j_cart, vcf.s.lattice.matrix, v_cutoff)
                dist = np.linalg.norm(i_cart - j_cart[0])
                if j_cart[1] < v_cutoff:
                    bc_bonds.append([element1,element2,dist])
                    cs_carts.append([j_cart[0], element1, element2])

        for i in range(len(cs_carts)):
            a = np.array(cs_carts[i][0])
            b = np.array(con_sites[M_idx].coords)
#           print(a, b, np.linalg.norm(b - a))
    num_c_sites_to_b = len(bc_bonds)
    num_b_sites = len(b_sites_idx)
    avg_num_coordinated_bc = num_c_sites_to_b / num_b_sites
    return avg_num_coordinated_bc

def OctDeform(b_atom, c_atom, contcar, v_cutoff, lapa):
   
    def UniqueVec(list):
        rot_axis = [list[0]]
        for i in range(1,len(list)):       
            unique_vec = True
            for j in range(len(rot_axis)):
                ab =  list[i]
                ac =  rot_axis[j]
                cosine_angle = np.dot(ab, ac) / (np.linalg.norm(ab) * np.linalg.norm(ac))
                if cosine_angle > 0.9:
                    unique_vec = False 
            if unique_vec == True:
                rot_axis.append(list[i])
        return rot_axis

    def OctElong(list,d0):   
       #print(list)
       #print("the average is: ")
        sum_de = 0
        for i in range(len(list)):
           sum_de += (list[i]/d0)**2
        lambda_oct = sum_de / 6
        return lambda_oct

    def OctAngVar(list):
        oct_axses = []
        angle_i = []
        for i in range(len(list)):
            for j in range(i+1,len(list)):
                        a = list[i][3]
                        b = list[i][4]
                        c = list[j][4]
                        ab = b - a
                        ac = c - a
                        bc = c - b
                        cosine_angle = np.dot(ab, ac) / (np.linalg.norm(ab) * np.linalg.norm(ac))
                        angle = np.arccos(cosine_angle)
                        angle_i.append([bc,np.degrees(angle)])
        angle_i.sort(key=lambda x:x[1])
        vec_a = (angle_i[-1][0])
        vec_b = (angle_i[-2][0])
        vec_c = (angle_i[-3][0])
        oct_axses.append(vec_a)
        oct_axses.append(vec_b)
        oct_axses.append(vec_c)
      
#       oct_axses = UniqueVec(oct_axses)
        lapa_A, lapa_B, lapa_C = (lapa.matrix[i] for i in range(3))
        lapa_ABC = [lapa_A, lapa_B, lapa_C]
        oct_tilt = OctTilt(oct_axses, lapa_ABC)
        
        angles_for_var = []
        for i in range(12):
            angles_for_var.append(angle_i[i][1])
        
        sigma_sum = 0 
        for i in range(len(angles_for_var)):
            sigma_sum += ( angles_for_var[i] - 90 ) ** 2
        sigma_oct_square = sigma_sum / 11

        phi_sum = 0
        if oct_tilt != None:
            for i in range(3):
                phi_sum += oct_tilt[i] ** 2
            phi_deviation = math.sqrt(phi_sum) / 3
        else:
            phi_deviation = None

        return sigma_oct_square, phi_deviation

    vcf = VoronoiConnectivity(contcar, cutoff=v_cutoff)
    b_sites_idx = []
    con_sites = contcar.sites
    lambda_oct_i = []
    sigma_oct_i = []
    phi_oct_i = []
    long_bonds = 0
    for i in range(len(con_sites)):
        if con_sites[i].specie.symbol == b_atom:
            b_sites_idx.append(i)
    num_c_sites_to_b = 0.0
    num_b_sites = len(b_sites_idx)
    for M_idx in b_sites_idx:
        v_coordinated_site = vcf.get_connections()
        num_bc_i = 0
        bc_bonds = []
        for b_i in range(len(v_coordinated_site)):
            element1 = con_sites[v_coordinated_site[b_i][0]].specie.symbol
            element2 = con_sites[v_coordinated_site[b_i][1]].specie.symbol
            element1_idx = v_coordinated_site[b_i][0]
            element2_idx = v_coordinated_site[b_i][1]
            bond_pair = [element1, element2]
            if (element1 in b_atom) and (element2 in c_atom) and (element1_idx == M_idx):
                i_cart = con_sites[v_coordinated_site[b_i][0]].coords
                j_cart = con_sites[v_coordinated_site[b_i][1]].coords
                dist = np.linalg.norm(i_cart-j_cart)
                if dist < v_cutoff:
                    bc_bonds.append([element1,element2,dist, i_cart, j_cart])
                    num_c_sites_to_b += 1
                    num_bc_i += 1
            elif (element2 in b_atom) and (element1 in c_atom) and element2_idx == M_idx:
                i_cart = con_sites[v_coordinated_site[b_i][1]].coords
                j_cart = con_sites[v_coordinated_site[b_i][0]].coords
                j_cart = ImagejCart(i_cart, j_cart, vcf.s.lattice.matrix, v_cutoff)
                dist = np.linalg.norm(i_cart - j_cart[0])
                if j_cart[1] < v_cutoff:
                    bc_bonds.append([element1,element2,dist, i_cart, j_cart[0]])
                    num_c_sites_to_b += 1
                    num_bc_i += 1

        if num_bc_i == 6:
#           print("B--C bonds are: ")
#           print(bc_bonds[0:6])
            points_oct = [ bc_bonds[i][4] for i in range(6) ]
            print(points_oct)
#           points = np.array(points_oct)  # your points
            points = points_oct
            volume = ConvexHull(points).volume
            # d0 is the center-to-vertex distance of a regular polyhedron of the same volume 
            d0 = (volume * 3 / 4) ** (1.0/3) 
            bc_bonds_M_idx = bc_bonds[0:6]
            bc_lengths = [row[2] for row in bc_bonds_M_idx]
            bc_lengths.sort()
            print("The ratio of bonds beweet the shortest and longest: ")
            print((bc_lengths[0]/bc_lengths[5]))
            for bc_length in bc_lengths:
                if bc_length > (radii[b_atom] + radii[c_atom[0]]):
                   print("WARNING: the distance is longer than the anion radii") 
                   long_bonds = long_bonds + 1
            # build list of the multiple octahedrals in a supercelli
            lambda_oct_i.append(OctElong(bc_lengths,d0))
            sigma_oct_i.append(OctAngVar(bc_bonds)[0])
            phi_oct_i.append(OctAngVar(bc_bonds)[1])
        else:
            print(("Warning: found less or more anions than 6 around the transition metal: ", b_atom, M_idx, num_bc_i))
            
    lambda_oct_filter = []
    sigma_oct_filter = []
    phi_oct_filter = []
    for i in range(len(phi_oct_i)):
        if phi_oct_i[i] != None:
            lambda_oct_filter.append(lambda_oct_i[i])
            sigma_oct_filter.append(sigma_oct_i[i])
            phi_oct_filter.append(phi_oct_i[i])

    lambda_oct_avg = np.mean(lambda_oct_filter)
    sigma_oct_avg = np.mean(sigma_oct_filter)
    phi_oct_avg = np.mean(phi_oct_filter)
 
    if len(phi_oct_filter) != len(phi_oct_i):
        print(("Wanring, only ", len(phi_oct_filter), " octahedra are found."))

    return lambda_oct_avg, sigma_oct_avg, phi_oct_avg, long_bonds, len(phi_oct_filter)

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def OctTilt(octahedron, lattice):
    oct = [[],[],[]]
    lat = [[],[],[]]
    angles = []
    for i in range(3):
        octahedron[i] = [abs(x) for x in octahedron[i]]
        lattice[i] = [abs(x) for x in lattice[i]]
        oct_max_index = octahedron[i].index(max(octahedron[i]))
        lat_max_index = lattice[i].index(max(lattice[i]))
        oct[oct_max_index] = octahedron[i]
        lat[lat_max_index] = lattice[i]
    find_oct_axses = [] in oct
    if find_oct_axses != True:
        for i in range(3):
            angles.append(angle_between(oct[i], lat[i])*180/math.pi)
        return angles
    else:
        print(("Failed to find the axses of the octahedron", oct, " from the given octahedron, ", octahedron))
        return None

def rgbline(ax, k, e, red, green, blue, alpha=1.):
    # creation of segments based on
    # http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
    ######################################
    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
 
    nseg = len(k) - 1
    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linewidth=3)
    ax.add_collection(lc)
 

def get_convex_hull(material_composition):
    from pymatgen.analysis.phase_diagram import PhaseDiagram
    from pymatgen.analysis.phase_diagram import PDPlotter

    import collections
    from pandas import DataFrame

    entries_list=[]
    for i in comp.elements:
      entries_list.append(i.name)

    a = MPRester("KvRs3wlFimrk7Hyq")
    entries = a.get_entries_in_chemsys(entries_list)

    pd = PhaseDiagram(entries)
    print(pd)

    plotter = PDPlotter(pd, show_unstable=False)

    ## Analyze decomposition
    hull = pd.get_hull_energy(comp)
    print("Composition: ", comp)
    print("Convex Hull (Not normalized by atoms): ",hull)
    depos = pd.get_decomposition(comp)
    decomp_data = collections.defaultdict(list)
    for e in depos:
        decomp_data["Entry"].append(e.entry_id)
        decomp_data["ComputedEnergy"].append(e.energy)
        decomp_data["Uncorrected_CE"].append(e.uncorrected_energy)
        decomp_data["Correction"].append(e.correction)
        decomp_data["Composition"].append(str(e.composition))
        print(e.composition)
    decomp_df = DataFrame(decomp_data, columns=["Entry", "ComputedEnergy", "Uncorrected_CE", "Correction", "Composition"])
    print(decomp_df)

    data = collections.defaultdict(list)
    for e in entries:
        decomp, ehull = pd.get_decomp_and_e_above_hull(e)
        if ehull == 0:
            data["Materials ID"].append(e.entry_id)
            data["Composition"].append(e.composition.reduced_formula)
            data["CompEnergy"].append(e.energy)
            data["Ehull"].append(ehull)


    ##print data
    df = DataFrame(data, columns=["Materials ID", "Composition", "Ehull", "CompEnergy", "Decomposition"])
    ##print(df.to_string())
    print(df)
    return hull


@explicit_serialize
class PerovskiteAnalysis(FiretaskBase):
    """
    Analysis the computed perovskite strcuture.
    This calculates the octahedral deformation,
    and store the dataframe as a pandas piclke file.
    """
   
    _fw_name = "Analysis Task"

    def _load_params(self, d):

        self.concentration = d['concentration']
        self.angle= d['angle']
        self.dir_name=d['dir_name']



    def run_task(self, fw_spec):
        from pathlib import Path

        if self.get("use_global_spec"):
            self._load_params(fw_spec)
        else:
            self._load_params(self)

        # infile = open(self.path_+'POSCAR', 'r')
        # infile = open('/home/heesoo/high-throughput_wahab/MA-EA-AbdulWAhab/dummy/POSCAR', 'r')
        infile = open('POSCAR', 'r')

        firstLine = infile.readline()
        print(firstLine)
        theComp = firstLine

        # contcar=mg.core.structure.Structure.from_file(self.path_+"CONTCAR")
        # contcar=mg.core.structure.Structure.from_file("/home/heesoo/high-throughput_wahab/MA-EA-AbdulWAhab/dummy/CONTCAR")
        contcar=mg.core.structure.Structure.from_file("CONTCAR")

        composition1 = str(contcar.composition)
        composition2 = composition1.replace(' ','-')

        #Components for A, B and C sites
        composition3 = Composition(composition1) ## This can be used by  anonymized_formula
        comp = composition3.reduced_composition
        comp_amt = comp.get_el_amt_dict()
     

        vasprun = Vasprun("vasprun.xml")
        # vasprun = Vasprun(self.path_ + "vasprun.xml")
        # vasprun = Vasprun("/home/heesoo/high-throughput_wahab/MA-EA-AbdulWAhab/dummy/vasprun.xml")
        entry = vasprun.get_computed_entry()

        # set up matplotlib plot
     
        try:
            theComp
        except NameError:
            print('theComp is not defiend. entry.composition is stored')
            theComp = entry.composition

        # contcar = mg.core.structure.Structure.from_file("/home/heesoo/high-throughput_wahab/MA-EA-AbdulWAhab/dummy/CONTCAR")
        # contcar = mg.core.structure.Structure.from_file(self.path_+"CONTCAR")

        contcar = mg.core.structure.Structure.from_file("CONTCAR")

        lapa = contcar.lattice 
        
        #systemname was gernerated accodring to atoms.inp
        systemname=TheCompCell(theComp)
        print(systemname.comp)
        mainA = systemname.mainA
        print(systemname.subA)
        subA = systemname.subA
        print(systemname.subAfrac)
        print(systemname.numcell)
        theCell = systemname.cell
        num_of_cell = list(map(int, theCell.split('x')))
        num_of_cell = num_of_cell[0] * num_of_cell[1] * num_of_cell[2]
        num_of_sites = num_of_cell * 5

        num_of_mainA = round(num_of_cell * (1 - systemname.subAfrac))
        num_of_subA  = round(num_of_cell * systemname.subAfrac)
        # REVISION BY ABDULWAHAB: Removed num_of_cell to remove formula unit conversion in CorrectedEnergy and UncorrectedEnergy

        CorrectedEnergy = entry.energy  # (corrected total energy )  from DOS
        UncorrectedEnergy = entry.uncorrected_energy # (uncooredted total energy ) from DOS
        Correction = entry.correction
        allcomps = re.findall('[A-Z][a-z]*', systemname.comp)
        print(allcomps)
        B_site_atom = allcomps[1]
        C_site_atom = [allcomps[2]]
        print(("NUM_OF_CELL is ",num_of_cell))
        if num_of_cell == 1:
            if allcomps[2] == allcomps[3]:
               print((allcomps[0]+allcomps[1]+allcomps[2]+"3"))
               theComp = allcomps[0]+"-"+allcomps[1]+allcomps[2]+"3"
               theComp2 = allcomps[0]+allcomps[1]+allcomps[2]+"3"
            else:
                C_site_atom.append(allcomps[3])
        else:
            if allcomps[2] == allcomps[3]:
               print((allcomps[0]+str(num_of_cell)+allcomps[1]+str(num_of_cell)+allcomps[2]+str(3*num_of_cell)))
               theComp = allcomps[0]+str(num_of_mainA)+subA+str(num_of_subA)+'-'+allcomps[1]+str(num_of_cell)+allcomps[2]+str(3*num_of_cell)
               theComp2 = allcomps[0]+str(num_of_mainA)+subA+str(num_of_subA)+allcomps[1]+str(num_of_cell)+allcomps[2]+str(3*num_of_cell)
            else:
                C_site_atom.append(allcomps[3])

        print(("Total energy is :", entry.energy))
        print(("Total energy (uncorrected) is :", entry.uncorrected_energy))
     
        voronoi_cutoff_list = np.arange(3.3, 5.0, 0.1)
        voronoi_cutoff = 0
        for i in voronoi_cutoff_list: 
            print(("Voronoi cutoff is ",i))
            voronoi_cutoff = i
            bc_coord_num = AvgCoordNum(B_site_atom, C_site_atom, contcar, voronoi_cutoff)    
            if bc_coord_num >= 6 :
                break
        print(("Average coordination number around the B sites: ", bc_coord_num))
        #OctDeform 
        lambda_oct, sigma_oct , phi_oct, long_bonds, num_oct = OctDeform(B_site_atom, C_site_atom, contcar, voronoi_cutoff, lapa)    
        print(("Number of octahedra in the sueprcell: ", num_oct))
        print(("Mean Octahedral Elongation and Octahedra angle variance, and Tilted angle: ", lambda_oct, "and ", sigma_oct, "and ", phi_oct))

        line1 = "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        line2 = 'Composition   |super    |CorrectedE/FU|UncorE/FU|Correc- |Num of   | MOE   | OAV  |Tilt | y_A2          '
        line3 = '              | -cell   |     (eV)    |   (eV)  |tion(eV)|B-C Coord|       |      |dev  |(concentraion) '
        line4 = "--------------------------------------------------------------------------------------------------------------------------------------------------"
        line5 = '{0:<15} {1:<7} {2:^13.5f} {3:^9.5f} {4:^7.3f} {5:^7.3f} {6:^9.4f} {7:^7.3f} {8:^7.3f} {9:^7.4}'\
                .format(theComp, theCell, CorrectedEnergy,\
                UncorrectedEnergy, Correction, bc_coord_num,\
                lambda_oct, sigma_oct, phi_oct, self.concentration)
        line7 = "D band gap (Direct band gap); Form E (Formation energy); Uncor E (Uncorrected energy); MOE (Mean Octahedral Elongation); OAV (Octahedra angle variance)\n"
        print(line1)
        print(line2)
        print(line3)
        print(line4)
        print(line5)
        print(line7)
        print(line1)

        sum_file = open('summary.txt', 'w')
        sum_file.write(line1+'\n')
        sum_file.write(line2+'\n')
        sum_file.write(line3+'\n')
        sum_file.write(line4+'\n')
        sum_file.write(line5+'  '+'\n')
        sum_file.close()
        print("TotEng",[UncorrectedEnergy])
        print("TotEng",float(UncorrectedEnergy))
        print("TotEng",8*((num_of_mainA)*8+(num_of_subA)*11+4))
        total_atoms=8*((num_of_mainA/(num_of_mainA+num_of_subA))*8+(num_of_subA/(num_of_mainA+num_of_subA))*11+4)
        print("TOTAL ATOMS",total_atoms,UncorrectedEnergy,(UncorrectedEnergy)/total_atoms)
        # REVISION BY ABDULWAHAB: Total Atoms calculated by supercell * ((number of A1 * 8) + (number of A2 * 11) + 4 ) -> 2*2*2(A1 * 8 + A2 * 11 + 4)
        print("Converged",[vasprun.converged])
        working_directory=os.getcwd()
        summary = {"Comp":[theComp], "SuperCell":[theCell],\
                   "A1":[mainA], "A2":[subA], "numA1":[num_of_mainA], "numA2":[num_of_subA],\
                    "numAllCat": [num_of_mainA+num_of_subA], "yA2":[float(self.concentration)],"yA2Actual":[num_of_subA/(num_of_mainA+num_of_subA)],"TotalAtoms":[total_atoms],\
                   "TotEng":[UncorrectedEnergy],"Converged":[vasprun.converged],\
                   "lambda":[lambda_oct], "sigma":[sigma_oct], "tilt":[phi_oct],"angle":[self.angle],"working_directory":[working_directory],"material_id":[str(theComp)+"_"+str(theCell)+"_"+str(mainA)+"_"+str(num_of_subA/(num_of_mainA+num_of_subA))+"_"+str(total_atoms)+"_"+str(lambda_oct)+"_"+str(sigma_oct)+"_"+str(phi_oct)+"_"+str(self.angle)],"TotEngperAtom": [(UncorrectedEnergy)/total_atoms] }
        print(summary)

        print("Saving as csv...")

        df = pd.DataFrame(summary)
        df.index.name='index'
        print("path to save pickle",self.dir_name+'/perovskites.pkl')

        import os
        root = os.path.abspath(os.path.join(os.getcwd(), os.path.pardir))

        my_file = Path(root + '/' + self.dir_name +"/perovskites.pkl")
        if my_file.is_file():
            df_from_file = pd.read_pickle(root + '/' + self.dir_name +"/perovskites.pkl")
            df_from_file = df_from_file.append(df, ignore_index=True)
            pd.to_pickle(df_from_file, root + '/' + self.dir_name +"/perovskites.pkl")
            df_from_file.to_csv(root + '/' + self.dir_name +"/perovskites.csv")
        else:
            pd.to_pickle(df, root + '/' + self.dir_name +"/perovskites.pkl")
            df.to_csv(root + '/' + self.dir_name +"/perovskites.csv")

        return FWAction(stored_data=summary, mod_spec=[])

@explicit_serialize
class DeltaEMixfromPickle(FiretaskBase):
    """
    Calculate $\Delta E _\mathrm{mix}$,
    reading the dataframe, perovskites.pkl
    """
   
    _fw_name = "MixEenergy Task"

    def _load_params(self, d):

        self.concentration = d['concentration']


    def run_task(self, fw_spec):
        if self.get("use_global_spec"):
            self._load_params(fw_spec)
        else:
            self._load_params(self)

        from pathlib import Path
        import pandas as pd

        try:
            root = os.path.abspath(os.path.join(os.getcwd(), os.path.pardir))

            my_file = Path(root + '/' + self.dir_name +"/perovskites.pkl")
            df = pd.read_pickle(root + '/' + self.dir_name +"/perovskites.pkl")
        except:
            print('Cannot find perovksite.pkl database to calculate the enthlapy of cation mix')
        
        # yA2 = self.concentration
        yA2 = num_of_subA/(num_of_mainA+num_of_subA)
        print(df) 
        returntest = df['yA2'].mean()
        recentdf = df.iloc[[-1]]
        print("Recent df is...")
        print(recentdf) 

        # Energy of pure perovskites and set their minimum as the reference
        # Comp_A1 = df[(df['yA2'] < (0.0 + 1/df['numAllCat']) ) & (df['Converged']==True)]
        Comp_A1 = df[(df['yA2'] == 0.0 ) & (df['Converged']==True)]
        """problem above """
        print("INSIDE MIXCAT, Comp_A1",Comp_A1)
        MinE_A1 = Comp_A1['TotEng'].min()
        print("INSIDE MIXCAT, MinE_A1",MinE_A1)

        # Comp_A2 = df[(df['yA2'] > (1.0 - 1/df['numAllCat']) ) & (df['Converged']==True)]
        Comp_A2 = df[(df['yA2']== 1.0 ) & (df['Converged']==True)]

        print("INSIDE MIXCAT, Comp_A2",Comp_A2)

        MinE_A2 = Comp_A2['TotEng'].min()
        print("INSIDE MIXCAT, MinE_A2",MinE_A2)

        recent_yA2 = df.at[df.index[-1],'yA2']
        print("INSIDE MIXCAT, recent_yA2", recent_yA2)

        enthalpy_comp = df.at[df.index[-1],'TotEng']
        print("INSIDE MIXCAT, enthalpy_comp", enthalpy_comp)

        enthalpy_ref = (MinE_A1 * (1.0 - yA2)) + (MinE_A2 * yA2)
        print("INSIDE MIXCAT, enthalpy_ref", enthalpy_ref)

        enthalpy_mix = (enthalpy_comp - enthalpy_ref)  * 1000    # in eV/FU
        print("INSIDE MIXCAT, enthalpy_mix", enthalpy_mix)

        print("Computed Dleta H_Mix (meV/F.U.): ", enthalpy_mix)
        
        #df_Emix = df[['Comp','A1','A2','yA2','TotEng']]
        df_Emix = df.loc[:,('Comp','A1','A2','yA2','TotEng')]
        print("INSIDE MIXCAT, df_Emix", df_Emix)

        df_Emix['DHmix'] = df_Emix['TotEng'] - MinE_A1 * (1.0 - df_Emix['yA2']) - MinE_A2 * (df_Emix['yA2'])
        print("INSIDE MIXCAT, df_Emix['DHmix']", df_Emix['DHmix'])

        df_Emix['DHmix_in_meV'] = df_Emix['DHmix'] * 1000      # in meV/FU
        df_Emix.index.name='index'
        df_Emix.to_csv('../DHmix.csv')

        return FWAction(stored_data={"DHmix":[df_Emix.at[df_Emix.index[-1],'DHmix_in_meV']]}, mod_spec=[])

