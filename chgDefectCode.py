# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 13:36:08 2023

@author: evanp
"""

#!/usr/bin/env python3
import os

import ase.io
import ase.units

from spinney.structures.defectivesystem import DefectiveSystem
from spinney.defects.diagrams import extract_formation_energies_from_file, Diagram

group_direc = 'D:\Summer Research 2023\Spinney/'
path_ef = os.path.join(group_direc, 'data_EJP2')
path_defects = os.path.join(path_ef, 'data_defects')
path_pristine = os.path.join(path_ef, 'pristine', 'OUTCAR')
path_rb = os.path.join(path_ef, 'Rb', 'OUTCAR')
path_sb = os.path.join(path_ef, 'Sb', 'OUTCAR')
path_ef_Rb3Sb2I9txt = os.path.join(path_ef, 'formation_energies_Rb3Sb2I9.txt')
path_translvl = os.path.join(path_ef, 'transition_levels.txt')
path_translvl_ext = os.path.join(path_ef, 'transition_levels_extended.txt')
path_diagram = os.path.join(path_ef, 'Defect Graph Point 3')

# initialize the defective system, where calculations have been done with VASP
defective_system = DefectiveSystem(path_ef, 'vasp')
ase_pristine = ase.io.read(path_pristine, format='vasp-out')
ase_rb = ase.io.read(path_rb, format='vasp-out')
ase_sb = ase.io.read(path_sb, format='vasp-out')

dmu_rb = -2.487 #-2.487 -2.629 -3.530 -3.437
dmu_sb = 0     #0 0 -2.704 -2.850

# calculate chemical potentials
mu_prist = 14*ase_pristine.get_total_energy()/ase_pristine.get_global_number_of_atoms()
mu_rb = ase_rb.get_total_energy()/ase_rb.get_global_number_of_atoms() # Ga-rich
mu_sb = ase_sb.get_total_energy()/ase_sb.get_global_number_of_atoms() # Ga-rich
mu_i = mu_prist - 3*mu_rb - 2*mu_sb # N-poor
mu_rb = mu_rb + dmu_rb
mu_sb = mu_sb + dmu_sb
print(mu_prist, mu_rb, mu_sb, mu_i)

# feed the data to the instance
defective_system.chemical_potentials = {'Rb3':mu_rb, 'Sb2':mu_sb, 'I9':mu_i}

print(defective_system.chemical_potentials)

# eigenvalue of the valence band maximum
vbm = 0.973503  # 2.7515
# calculated dielectric tensor
e_rx = 12.46
e_rz = 16.47
e_r = [[e_rx, 0, 0], [0, e_rx, 0], [0, 0, e_rz]]

# feed the data
defective_system.vbm = vbm
defective_system.dielectric_tensor = e_r
defective_system.correction_scheme = 'ko'  # Kumagai and Oba
# defective_system.correction_scheme = 'fnv'  # Freysoldt, Neugebauer and Van de Walle

# calculate defect formation energies for each point defect
defective_system.calculate_energies(verbose=False)
df = defective_system.data # data frame with calculated formation energies
print(df)

defective_system.write_formation_energies(path_ef_Rb3Sb2I9txt)

defective_system.gap_range = (0, 1.960909)  # E_F from 0 to PBE band gap value
defective_system.extended_gap_range = (-0.85, - 0.85 + 3.51)  # Fermi level will range in the HSE range
# calculate charge transition levels
defective_system.diagram.transition_levels

#print defect transisition levels

#data_file = 'formation_energies_GaN.txt'
defect_dictionary = extract_formation_energies_from_file(path_ef_Rb3Sb2I9txt )
dgm = Diagram(defect_dictionary, (0, 1.960909)) # E_F from 0 to PBE to HSE band gap value
dgm.write_transition_levels(path_translvl_ext)

# use some prettier labels in the plot
defective_system.diagram.labels = {
    r'$I_{Rb}$':'I1',
    r'$Sb_{Rb}$':'I2',
    r'$Rb_I$':'I3',
    r'$Sb_I$':'I4',
    r'$I_{Sb}$':'I5',
    r'$Rb_{Sb}$':'I6',
    r'$Vac_{I}$':'I7',
    r'$Vac_{Rb}$':'I8',
    r'$Vac_{Sb}$':'I9'
}

# personalize colors to use in the plot
ef_colors = {
    r'$I_{Rb}$':'green',
    r'$Sb_{Rb}$':'red',
    r'$Rb_I$':'blue',
    r'$Sb_I$':'purple',
    r'$I_{Sb}$':'orange',
    r'$Rb_{Sb}$':'pink',
    r'$Vac_{I}$':'cyan',
    r'$Vac_{Rb}$':'brown',
    r'$Vac_{Sb}$':'black'
}

defective_system.diagram.plot(save_flag=True, save_title=path_diagram,
                              title=r'$Rb_3Sb_2I_9$  Charge  Defects', legend=True,
                              colors_dict=ef_colors, x_label=r'$E_F$ (eV)')

#-37.09789452374999 -3.4764014125 -4.05931461 -26.011061066249994
#-37.09789452374999 -3.6184014125 -4.05931461 -26.011061066249994
#-37.09789452374999 -4.5194014125 -6.76331461 -26.011061066249994