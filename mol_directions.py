#!/usr/bin/env python
# -*- coding=utf-8 -*-

from pymatgen.core.structure import Structure
#from pymatgen.analysis.structure_matcher import StructureMatcher
import os
import fnmatch
#from pymatgen.analysis.graphs import StructureGraph
from custom_mol import StructureGraph
 # used for deciding which atoms are bonded
from pymatgen.analysis.local_env import JmolNN
import numpy as np
import math as m
from astropy.coordinates import cartesian_to_spherical
import pandas as pd 


df_vol = pd.read_csv('deform_volume.csv')


contcar_count = len(fnmatch.filter(os.listdir("."), 'CONTCAR*.vasp'))
cc = contcar_count
connum = '{:04d}'.format(cc-1)
lastfile = "CONTCAR"+connum+".vasp"
print("Last file is: ",lastfile)
structures = []


def get_mol_structure(s1):
    # Translate the atoms to find 8 molecules
    s1.translate_sites(range(len(s1)), (0.25,0.25,0.25))
    # Remove inorganic network frame
    s1.remove_species(['Pb', 'I'])
    
    # Mol direction C -> N    
    sg = StructureGraph.with_local_env_strategy(s1, JmolNN())
    my_molecules = sg.get_subgraphs_as_molecules()
    return my_molecules            

def get_longdist_pair(mol, N, C):
    Nindx = mol.indices_from_symbol(N)
    Cindx = mol.indices_from_symbol(C)
    longdist = 0.0
    longdistpair = []
    for n in Nindx:
        for c in Cindx:
            CN_dist = mol.get_distance(n,c)
            if CN_dist > longdist:
                longdist = CN_dist
                longdistpair = [n,c]
    vec = mol.sites[longdistpair[1]].coords - mol.sites[longdistpair[0]].coords
    vec_len = np.linalg.norm(vec)
    vec = vec / vec_len
    return (n,c,vec)


def op_cos_theta(ref_vec, vec):
    """ Dot product of two vectors in spherical coordinates
    return cos(omega)"""
    sin_th1 = np.sin(ref_vec[0])
    cos_th1 = np.cos(ref_vec[0])
    sin_th2 = np.sin(vec[0])
    cos_th2 = np.cos(vec[0])
    cos_ph1ph2 = np.cos(ref_vec[1]-vec[1])
    cos_omega = sin_th1 * sin_th2 * cos_ph1ph2 + cos_th1 * cos_th2
    return cos_omega


def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = m.sqrt(XsqPlusYsq + z**2)            # r
    elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # theta
    az = m.atan2(y,x)                        # phi
    #if az < 0:
    #    elev = (-1.0) * elev
    #    az = az + m.pi
    return r, elev, az


def get_directions(my_molecules):
    cations=[]
    sphcoords_theta=[]
    sphcoords_phi=[]
    # The reference vector of alignment order parameter [theta,phi]
    OP_ref_vec = [np.pi/4, np.pi/4]
    # OP_ref_vec = [0.0, 0.0]
    OPcos = []
    for i in my_molecules:
        rdfml = i.composition.reduced_formula
        if rdfml == 'H6CN':
            nc = get_longdist_pair(i,'N','C')
            x, y, z =  nc[2]
            cations.append('MA')
            r, theta, phi = cart2sph(x,y,z)
            r, theta, phi = cartesian_to_spherical(x,y,z)
            op = op_cos_theta(OP_ref_vec, [theta,phi])
            sphcoords_theta.append(theta.radian)
            sphcoords_phi.append(phi.radian - m.pi)
            OPcos.append(op**2)
            
        elif rdfml == 'H8C2N':
            nc = get_longdist_pair(i,'N','C')
            x, y, z =  nc[2]
            cations.append('EA')
            r, theta, phi = cart2sph(x,y,z)
            r, theta, phi = cartesian_to_spherical(x,y,z)
            op = op_cos_theta(OP_ref_vec, [theta,phi])
            sphcoords_theta.append(theta.radian)
            sphcoords_phi.append(phi.radian - m.pi)
            OPcos.append(op**2)
            
        elif rdfml == 'H7C2NF':
            nc = get_longdist_pair(i,'N','F')
            x, y, z =  nc[2]
            cations.append('FEA')
            r, theta, phi = cart2sph(x,y,z)
            r, theta, phi = cartesian_to_spherical(x,y,z)
            op = op_cos_theta(OP_ref_vec, [theta,phi])
            sphcoords_theta.append(theta.radian)
            sphcoords_phi.append(phi.radian - m.pi)
            OPcos.append(op**2)
            

    MA_n = float(len([1 for i in cations if i == 'MA']))
    EA_n = float(len([1 for i in cations if i == 'EA']))
    FEA_n = float(len([1 for i in cations if i == 'FEA']))
    print("MA, EA, FEA: ", MA_n, EA_n, FEA_n)

    return [cations,MA_n,EA_n,FEA_n], sphcoords_theta, sphcoords_phi, OPcos


##-----------
# Start main
##-----------
tigercolors = ["#FA8D2F", "#84A660", "#1A99E9"]

structures = []
volumes = []
for i in range(0,cc):
    print("Read file: CONTCAR" + str(i+1) + ".vasp" )
    connum = '{:04d}'.format(i+1)
    contcar_struct = Structure.from_file("CONTCAR" + connum  + ".vasp")
    structures.append(contcar_struct)
    contcar_volume = contcar_struct.volume
    volumes.append(round(contcar_volume,3))

columns = ['CONTCAR', 'Cation', 'Theta', 'Phi', 'OPcos', 'Volume', 'dEform',
           'lambda', 'sigma', 'tilt', 'nMA', 'nEA', 'nFEA']
df = pd.DataFrame(columns=columns)
df = df[0:0]

for i in range(0,cc):
    my_molecules = get_mol_structure(structures[i])
    cations, sphcoords_theta, sphcoords_phi, opcos = get_directions(my_molecules)
    connum = '{:04d}'.format(i+1)
    filename = "CONTCAR" + connum + ".vasp"     
    filename = [filename] * len(cations[0])
    # match dEform from a external file by matching the volume
    volume = [volumes[i]] * len(cations[0])
    dEform = df_vol[df_vol['CellVol'] ==  volumes[i] ]['dEform'].values
    dEform = [dEform[0]] * len(cations[0])
    elongation = df_vol[df_vol['CellVol'] ==  volumes[i] ]['lambda'].values
    elongation = [elongation[0]] * len(cations[0])
    anglevariance = df_vol[df_vol['CellVol'] ==  volumes[i] ]['sigma'].values
    anglevariance = [anglevariance[0]] * len(cations[0])
    tilt = df_vol[df_vol['CellVol'] ==  volumes[i] ]['tilt'].values
    tilt = [tilt[0]] * len(cations[0])
    index = [str(i+1)] * len(cations[0])

    dict_to_add = {columns[0]: filename, columns[1]: cations[0],
            columns[2]: sphcoords_theta, columns[3]: sphcoords_phi,
            columns[4]: opcos, columns[5]: volume, 
            columns[6]: dEform, columns[7]: elongation, 
            columns[8]: anglevariance, columns[9]: tilt,
            columns[10]: [cations[1]]*8,
            columns[11]: [cations[2]]*8,
            columns[12]: [cations[3]]*8
            }
    print(dict_to_add)
    df_to_add = pd.DataFrame(dict_to_add, index=index)
    
    # MA/EA <- cations[3]==0, MA/FEA <- cations[2]==0
    # if (cations[2]) == 0 and (cations[3] == 8):
    # if (cations[2]) == 0:
    if True:                  # all compounds
        # print("add: ", cations)
        df = df.append(df_to_add)
    


print(df)
print(df.describe())
#df = df[df['Cation']=='MA']
    

print(df.describe())
xticks = [(-1)*m.pi, (-0.5)*m.pi, 0, 0.5*m.pi, m.pi]
yticks = [(-0.5)*m.pi, 0, 0.5*m.pi]
ax = df.plot.hexbin(x="Phi", y="Theta", gridsize=20, sharex=False,
#               cmap="viridis",
               xticks=xticks, yticks=yticks,
               xlim=[-m.pi,m.pi],
               ylim=[-0.5*m.pi,0.5*m.pi]);
ax.set_xticklabels(['-$\pi$', '-1/2$\pi$', '0', '1/2$\pi$', '$\pi$'],
                   rotation=0)
ax.set_yticklabels(['-1/2$\pi$', '0', '1/2$\pi$'],
                   rotation=0)