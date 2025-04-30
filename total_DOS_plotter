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

VBM = 0.95
CBM = 2.93

numAtoms  = POSCAR[6].split()
for i in range (len(numAtoms)):
    numAtoms[i] = int(numAtoms[i])

totalAtoms=[]
for i in range (len(numAtoms)+1):
    if (i>0):
        totalAtoms.append(totalAtoms[i-1] + numAtoms[i-1])
    else:
        totalAtoms.append(0)
print(totalAtoms)


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
    for n in range (totalAtoms[atom], totalAtoms[atom+1]):
        for i in range (301):
            if(n==totalAtoms[atom]):
                dosS.append(float(DOSCAR[i+6+301+1][1]))
            else:
                dosS[i] = dosS[i] + float(DOSCAR[i+6+302*n][1])
    return (dosS)

def pOrbital(atom):         
    dosP = []
    for n in range (totalAtoms[atom], totalAtoms[atom+1]):
        for i in range (301):
            if(n==totalAtoms[atom]):
                dosP.append((float(DOSCAR[i+6+301+1][2])+float(DOSCAR[i+6+301+1][3])+float(DOSCAR[i+6+301+1][4]))/3)
            else:
                dosP[i] = dosP[i] + (float(DOSCAR[i+6+302*n][2])+float(DOSCAR[i+6+302*n][3])+float(DOSCAR[i+6+302*n][4]))/3
    return (dosP)

def dOrbital(atom):
    dosD = []
    for n in range (totalAtoms[atom], totalAtoms[atom+1]):
        for i in range (301):
            if(n==totalAtoms[atom]):
                dosD.append((float(DOSCAR[i+6+301+1][5])+float(DOSCAR[i+6+301+1][6])+float(DOSCAR[i+6+301+1][7])+ 
                             float(DOSCAR[i+6+301+1][8]) + float(DOSCAR[i+6+301+1][9]))/5)
            else:
                dosD[i] = dosD[i] + (float(DOSCAR[i+6+302*n][5])+float(DOSCAR[i+6+302*n][6])+float(DOSCAR[i+6+302*n][7])
                                     +float(DOSCAR[i+6+302*n][8])+float(DOSCAR[i+6+302*n][9]))/5
        return (dosD)




def element1():    
    element1 = []
    for i in range(301):
        element1.append(sOrbital(0)[i] + pOrbital(atom)[i] + dOrbital(0)[i])
    return (element1)
    
def element2():     
    element2 = []
    for i in range(301):    
        element2.append(sOrbital(1)[i] + pOrbital(atom+1)[i] + dOrbital(1)[i])
    return (element2)
    
def element3():     
    element3 = []    
    for i in range(301):    
        element3.append(sOrbital(2)[i] + pOrbital(atom+2)[i] + dOrbital(2)[i])
    return (element3)

atom = 0

plt.figure(figsize=(10, 6))    
plt.plot(energy, dos, label='Total DOS', color = 'teal')
# plt.plot(energy, element1(), label=' Rb DOS', color = 'firebrick')
# plt.plot(energy, element2(), label='Sb DOS', color = 'purple')
# plt.plot(energy, element3(), label='I DOS', color = 'deeppink')
plt.plot(energy, sOrbital(atom), label='s Orbital', color = 'firebrick')
plt.plot(energy, pOrbital(atom), label='p Orbital', color = 'purple')
plt.plot(energy, dOrbital(atom), label='d Orbital', color = 'deeppink')

plt.axvline(0, color='black', linestyle='dashed')
plt.axvline(CBM - VBM, color='black', linestyle='dashed')

print(CBM - VBM)

# Customize the plot
plt.xlabel('Energy (eV)')
plt.ylabel('Density of States (# of States/eV)')
plt.title('Total Density of States in Rb3Sb2I9 with Defect (I to Sb)')
plt.legend()

# Display the plot
plt.show()
    
