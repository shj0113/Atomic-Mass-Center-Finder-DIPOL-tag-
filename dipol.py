#!/usr/bin/env python3

from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Poscar

import numpy as np
from scipy import ndimage
import os
import re

path = os.getcwd()
structure_name = "POSCAR"
accuracy = 3 

structure = Structure.from_file(path + "/" + structure_name)
structure_dict = structure.as_dict()

a = structure_dict['lattice']['a']
b = structure_dict['lattice']['b']
c = structure_dict['lattice']['c']

num_of_a_grid = round(a)*accuracy
num_of_b_grid = round(b)*accuracy
num_of_c_grid = round(c)*accuracy

a_grid = np.linspace(0,a,num=num_of_a_grid)
b_grid = np.linspace(0,b,num=num_of_b_grid)
c_grid = np.linspace(0,c,num=num_of_c_grid)

cell_array = np.zeros((len(c_grid), len(b_grid), len(a_grid)))

def find_close_point(target, points):
    
    # points 중 target과 가장 가까운 grid points의 index(=target_index) 찾기
    distance_list = []
    error_list = []
    for i in points:
        distance = np.sum((target - i)**2, axis=0)
        distance_list.append(distance)
    target_index = np.argmin(distance_list)
#     print('오차', distance_list[target_index])
    error_list.append(distance_list[target_index])
    
    # grid points의 index가 3d array에 어느 index에 들어가야 하는지 찾기
    grid_index = []
    num = 0
    for i,a_value in enumerate(a_grid):
        for j,b_value in enumerate(b_grid):
            for k,c_value in enumerate(c_grid):
                num+=1
                if num == target_index:
                    grid_index.append(k)
                    grid_index.append(j)
                    grid_index.append(i)
                    
    return grid_index, error_list

def grep_mass(element):

    with open(path+'/POTCAR', 'r') as f:
        lines = f.readlines()
        pattern = r'POMASS\s*=\s*([\d.]+)'

        for i, line  in enumerate(lines):
            if re.search('= PAW_PBE '+element, line ):
#                 print(i, line)
#                 print(lines[i+4])

                match = re.search(pattern, lines[i+4])
                if match:
                    pomass_value = float(match.group(1))
#                     print("POMASS value:", pomass_value)
    return pomass_value

grid_points = []
for i,a_value in enumerate(a_grid):
    for j,b_value in enumerate(b_grid):
        for k,c_value in enumerate(c_grid):
            grid_points.append(np.array([c_value,b_value,a_value]))
#             print(c_value,b_value,a_value)

atom_position_abc = []
for i in structure_dict['sites']:
    atom_position_abc.append(i['abc'])

## periodic boundary 처리하기

for num, sublist in enumerate(atom_position_abc):
    if str(sublist[0]).startswith('0.0') or str(sublist[1]).startswith('0.0') or str(sublist[2]).startswith('0.0'):
        atom_position_abc[num][0] = 0
        atom_position_abc[num][1] = 0
        atom_position_abc[num][2] = 0
    if str(sublist[0]).startswith('0.9') or str(sublist[1]).startswith('0.9') or str(sublist[2]).startswith('0.9'):
        atom_position_abc[num][0] = 0
        atom_position_abc[num][1] = 0
        atom_position_abc[num][2] = 0
        
for num, sublist in enumerate(atom_position_abc):
    atom_position_abc[num][0] = sublist[0]*structure_dict['lattice']['a']
    atom_position_abc[num][1] = sublist[1]*structure_dict['lattice']['b']
    atom_position_abc[num][2] = sublist[2]*structure_dict['lattice']['c']
    
count_pd_atoms = atom_position_abc.count([0, 0, 0])
print('Periodic boundary에 있는 원자 수: ', count_pd_atoms)

atom_position_a = np.array(atom_position_abc).T[0]
atom_position_b = np.array(atom_position_abc).T[1]
atom_position_c = np.array(atom_position_abc).T[2]

atom_position_cba = []
for c, b, a in zip(atom_position_c, atom_position_b, atom_position_a):
    atom_position_cba.append([c,b,a])

total_error = []
for num, i in enumerate(structure_dict['sites']):
    if not (atom_position_cba[num][0] == 0 and atom_position_cba[num][1] == 0 and atom_position_cba[num][2] == 0):
        grid_index, error_list = find_close_point(atom_position_cba[num], grid_points)
        total_error.append(error_list)
        cell_array[grid_index[0],grid_index[1],grid_index[2]] = grep_mass(i['species'][0]['element'])

center = list(ndimage.center_of_mass(cell_array))

Dipol_cartesian = []
for i,a_value in enumerate(a_grid):
    for j,b_value in enumerate(b_grid):
        for k,c_value in enumerate(c_grid):
            if k == round(center[0]) and j == round(center[1]) and i == round(center[2]):
                Dipol_cartesian.append(a_value)
                Dipol_cartesian.append(b_value)
                Dipol_cartesian.append(c_value)

Dipol_x = round(Dipol_cartesian[0]/structure_dict['lattice']['a'], 2)
Dipol_y = round(Dipol_cartesian[1]/structure_dict['lattice']['b'], 2)
Dipol_z = round(Dipol_cartesian[2]/structure_dict['lattice']['c'], 2)

print('DIPOL = ', Dipol_x, Dipol_y, Dipol_z)
