# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 22:59:16 2023

@author: evanp
"""

#Libraries
import numpy as np
import matplotlib.pyplot as plt
import math

#Functions
#Function to convert POSCAR into Angstoms
def convertToAngstrom(latticeX, latticeY, latticeZ, array1, array2, array3, arrayX, arrayY, arrayZ, latticeConstant, atoms):
    for a in range (0, atoms):
        x = latticeConstant*(latticeX[0]*array1[a] + latticeX[1]*array2[a] + latticeX[2]*array3[a])
        y = latticeConstant*(latticeY[0]*array1[a] + latticeY[1]*array2[a] + latticeY[2]*array3[a])
        z = latticeConstant*(latticeZ[0]*array1[a] + latticeZ[1]*array2[a] + latticeZ[2]*array3[a])
        arrayX.append(x)
        arrayY.append(y)
        arrayZ.append(z)
  
#Function to obtain factor in order to include periodicty in distance calculations        
def secondAtom(latticeX, latticeY, latticeZ, array1, array2, array3, period, latticeConstant, atoms):
    for a in range (0, atoms):
        o=0
        for b in (-1, 0, 1):
            for c in (-1, 0, 1):
                for d in (-1, 0, 1):
                    x = latticeConstant*(latticeX[0]*(array1[a]+b) + latticeX[1]*(array2[a]+c) + latticeX[2]*(array3[a]+d))
                    y = latticeConstant*(latticeY[0]*(array1[a]+b) + latticeY[1]*(array2[a]+c) + latticeY[2]*(array3[a]+d))
                    z = latticeConstant*(latticeZ[0]*(array1[a]+b) + latticeZ[1]*(array2[a]+c) + latticeZ[2]*(array3[a]+d))
                    period[o][a]=[x,y,z]
                    o = o + 1


#Makes Histograms for Same Element comparisons             
def distanceCalculator(array1, array2, array3, period, totalAtoms, i):
    array4 = []      
    for a in range (totalAtoms[i], totalAtoms[i+1]):
        for c in range (a+1, totalAtoms[i+1]):
            distance = []
            for b in range(27):
                distance.append(math.sqrt((array1[a]-period[b][c][0])**2+(array2[a]-period[b][c][1])**2+(array3[a]-period[b][c][2])**2))        
            w = min(distance)
            array4.append(w)
    
    return(array4)

#Makes Histograms for Same Element comparisons             
def difElements(array1, array2, array3, period, totalAtoms, i, j):
    array4 = []   
    for a in range (totalAtoms[i], totalAtoms[i+1]):
        for c in range (totalAtoms[j], totalAtoms[j+1]):
            distance = []
            for b in range(27):
                distance.append(math.sqrt((array1[a]-period[b][c][0])**2+(array2[a]-period[b][c][1])**2+(array3[a]-period[b][c][2])**2))        
            w = min(distance)
            if (w !=0):
                array4.append(w)
    
    return(array4)

#Beginning of Code
# poscar = "D:/Summer Research 2023/VESTA Files/MD/e01_solid-cd-Si/POSCAR"
# poscar = "D:/Summer Research 2023/VESTA Files/bulk/e04_cd-Si/POSCAR.vasp"
# poscar = "D:/Summer Research 2023/Rb3Sb2I9/Supercell/DefectCalc/Bulk/POSCAR"
poscar = "D:/Summer Research 2023/BiI3/ExperimentalPoscar/POSCAR"

data = np.loadtxt(poscar, skiprows=8, usecols=[0,1,2])
latticeInfo = np.loadtxt(poscar, skiprows=2, max_rows=3, usecols=[0,1,2])
latticeConstant = np.loadtxt(poscar, skiprows=1, max_rows=1)
atoms = np.loadtxt(poscar, skiprows=6, max_rows=1)


f = open(poscar)
POSCAR = f.readlines()

numAtoms  = POSCAR[6].split()
for i in range (len(numAtoms)):
    numAtoms[i] = int(numAtoms[i])

numberAtoms = 0
for i in range (len(numAtoms)):
    numberAtoms = numberAtoms + numAtoms[i]

totalAtoms=[]
for i in range (len(numAtoms)+1):
    if (i>0):
        totalAtoms.append(totalAtoms[i-1] + numAtoms[i-1])
    else:
        totalAtoms.append(0)
        
print(totalAtoms)

#Lattice Multiplication
latticeX = latticeInfo[:, 0]
latticeY = latticeInfo[:, 1]
latticeZ = latticeInfo[:, 2]

#Positions of Atoms (before converting)
position1 = data[:, 0]
position2 = data[:, 1]
position3 = data[:, 2]

#Positions of Atoms in Angstroms
arrayX = []
arrayY = []
arrayZ = []

#Contains [atom1, atom2, distance between the 2 atoms]
finaldata = []

convertToAngstrom(latticeX, latticeY, latticeZ, position1, position2, position3, arrayX, arrayY, arrayZ, latticeConstant, numberAtoms)

#contains information needed to include periodicty of cell
period = [[0 for x in range (numberAtoms)] for y in range(27)]

secondAtom(latticeX, latticeY, latticeZ, position1, position2, position3, period, latticeConstant, numberAtoms)

size = 20

plt.figure(figsize=(10, 6))  
plt.title("Nearest Neighbors for I to I", fontsize = size)
plt.xlabel("Relative Distance ($\AA$)", fontsize = size)
plt.ylabel("# of Atoms at Distance", fontsize = size)
plt.xticks(fontsize = size)
plt.yticks(fontsize = size)

finaldata = difElements(arrayX, arrayY, arrayZ, period, totalAtoms, 1, 1)

# plt.hist(distanceCalculator(arrayX, arrayY, arrayZ, period, totalAtoms, 1), bins = 50, color = 'teal')
plt.hist(finaldata, bins = 50, color = 'teal')
