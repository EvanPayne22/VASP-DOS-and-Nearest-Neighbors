# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:59:13 2023

@author: evanp
"""

import numpy as np
import matplotlib.pyplot as plt

# Load DOSCAR file
doscar_file = "D:/Summer Research 2023/Rb3Sb2I9/Supercell/DefectCalc/ItoSb/DOS/DOSCAR"
poscar = "D:/Summer Research 2023/Rb3Sb2I9/Supercell/DefectCalc/ItoSb/POSCAR"

f = open(doscar_file)
DOSCAR = f.readlines()

f = open(poscar)
POSCAR = f.readlines()

VBM = 0.98
CBM = 2.98

numAtoms  = POSCAR[6].split()
for i in range (len(numAtoms)):
    numAtoms[i] = int(numAtoms[i])

totalAtoms=[]
for i in range (len(numAtoms)+1):
    if (i>0):
        totalAtoms.append(totalAtoms[i-1] + numAtoms[i-1])
    else:
        totalAtoms.append(0)

for i in range(len(DOSCAR) - 6):
    DOSCAR[i+6] = DOSCAR[i+6].split()

energy=[]
for i in range (301):    
    energy.append(float(DOSCAR[i+6][0]) - VBM)

dos=[]
for i in range (301):    
    dos.append(float(DOSCAR[i+6][1]))

def sOrbital(atom):
    dosS = []
    for n in (atom):
        for i in range (301):
            if(n == atom[0]):
                z = 0
            if(z == 0):
                dosS.append(float(DOSCAR[i+6+302*n][1]))
            elif(z != 0):
                dosS[i] = dosS[i] + float(DOSCAR[i+6+302*n][1])
            z = z + 1
    return (dosS)

def pOrbital(atom):         
    dosP = []
    for n in (atom):
        for i in range (301):
            if(n == atom[0]):
                z = 0
            if (z == 0):
                dosP.append((float(DOSCAR[i+6+302*n][2])+float(DOSCAR[i+6+302*n][3])+float(DOSCAR[i+6+302*n][4]))/3)
            elif(z != 0):
                dosP[i] = dosP[i] + (float(DOSCAR[i+6+302*n][2])+float(DOSCAR[i+6+302*n][3])+float(DOSCAR[i+6+302*n][4]))/3
            z = z + 1
    return (dosP)

def dOrbital(atom):
    dosD = []
    for n in (atom):
        for i in range (301):
            if(n == atom[0]):
                z = 0
            if (z == 0):
                dosD.append((float(DOSCAR[i+6+302*n][5])+float(DOSCAR[i+6+302*n][6])+float(DOSCAR[i+6+302*n][7])
                                 +float(DOSCAR[i+6+302*n][8])+float(DOSCAR[i+6+302*n][9]))/5)
            elif(z != 0):
                dosD[i] = dosD[i] + (float(DOSCAR[i+6+302*n][5])+float(DOSCAR[i+6+302*n][6])
                                     +float(DOSCAR[i+6+302*n][7]) + float(DOSCAR[i+6+302*n][8])+float(DOSCAR[i+6+302*n][9]))/5
            z = z + 1
        return (dosD)


#Specifies element (1-3 depending on location in POSCAR)
atomType = 1
#Number of atom in element group
atomNumber = 14

atom1 = totalAtoms[atomType] + atomNumber


atom = [36,96]


plt.figure(figsize=(10, 6))    
# plt.plot(energy, dos, label='Total DOS', color = 'teal')
plt.plot(energy, sOrbital(atom), label='s Orbital', color = 'firebrick')
plt.plot(energy, pOrbital(atom), label='p Orbital', color = 'purple')
plt.plot(energy, dOrbital(atom), label='d Orbital', color = 'deeppink')

plt.axvline(0, color='black', linestyle='dashed')
plt.axvline(CBM - VBM, color='black', linestyle='dashed')

print(CBM - VBM)

size = 15

# Customize the plot
plt.xlabel('Energy (eV)', fontsize = size)
plt.ylabel('Density of States (# of States/eV)', fontsize = size)
plt.title('Density of States for the Defect Atom in $Rb_3Sb_2I_9$ ($Rb_{I}$)', fontsize = size)
plt.xticks(fontsize = size)
plt.yticks(fontsize = size)
plt.legend(fontsize = size)

# Display the plot
plt.show()
    
