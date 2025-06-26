# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 17:19:22 2025

@author: epayne
"""

import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Arguments for DOS visualization",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-doscar", nargs='?', default = "./DOSCAR", help="sets the location of the DOSCAR file")
parser.add_argument("-poscar", nargs='?', default = "./POSCAR", help="sets the location of the POSCAR file to determine elements")
parser.add_argument("-saveloc", nargs='?', default = "./doscar_plot.png", help="sets the save location of plot")
parser.add_argument("-colors", nargs='+', 
                    default=["red", "lightgreen", "dodgerblue", "plum",
                             "orange", "mediumorchid", "gold", "mediumaquamarine",
                             "deepskyblue", "salmon", "turquoise", "darkkhaki"]
                    , help="color array for charge neutrality plot")
parser.add_argument("-fontsize", nargs='?', type=int, default = 12, help="sets fontsize")
parser.add_argument("-NEDOS", nargs='?', type=int, default = 301, help="the NEDOS value set in VASP")
parser.add_argument("-plotheight", nargs='?', type=int, default = 4, help="the NEDOS value set in VASP")
parser.add_argument("-plotwidth", nargs='?', type=int, default = 8, help="the NEDOS value set in VASP")
parser.add_argument("-total", dest="total", action='store_true', help="show total DOS")
parser.add_argument("-nototal", dest="total", action='store_false', help="hide total DOS")
parser.set_defaults(total=True)
parser.add_argument("-orbitals", action='store_true', help="show orbitals (default: False)")
parser.add_argument("-elements", action='store_true', help="show elements (default: False)")
parser.add_argument("-subplots", action='store_true', help="splits elements and orbitals into subplots (default: False)")
args = parser.parse_args()
config = vars(args)


data = config["doscar"]
colors = config["colors"]
fontSize = config["fontsize"]
colorCounter = 0

f = open(data)
doscar = f.readlines()

if(config["elements"]):
    atomCounter = 0 #Counter to count atoms in DOSCAR
    elementCounter = 0 #Counter to seperate elements in DOSCAR
    totalElementDensity = [0] * config["NEDOS"] #Array for elemental contribution
    data2 = config["poscar"]
    f = open(data2)
    poscar = f.readlines()
    
    elementNames = poscar[5].split() #Array with the names of all the Elements in POSCAR
    numOfElements = [int(x) for x in poscar[6].split()] #Array with the number of all the Elements in POSCAR

# Determines the VBM to center band gap
VBM = float(doscar[5].split()[3])

energy = [] 
totalDensity = []
sDensity = [0] * config["NEDOS"] #Array for s-orbital contribution
pDensity = [0] * config["NEDOS"] #Array for p-orbital contribution
dDensity = [0] * config["NEDOS"] #Array for d-orbital contribution

# Settings for plot display
if(config["subplots"]):
    plt.subplots(nrows=(len(numOfElements) + 1), ncols=1, figsize=(config["plotwidth"]*len(numOfElements),config["plotwidth"]*len(numOfElements))) #Sets the number of subplots
else:
    plt.figure(figsize=(8, 4))

#Gets the Values for the Total DOS
for i in range (6, config["NEDOS"] + 6):
    energy.append(float(doscar[i].split()[0]) - VBM)
    totalDensity.append(float(doscar[i].split()[1]))

ymax = max(totalDensity)*1.1

# Plot and fill each curve
if(config["total"]):
    if(config["subplots"]):
        plt.subplot((len(numOfElements) + 1), 1, 1)
        plt.ylim(0, ymax)
        plt.xlim(energy[0], energy[-1])
        ax = plt.gca()
        plt.xlabel("Energy (eV)", fontsize=fontSize)
        plt.ylabel("Density of States", fontsize=fontSize, labelpad=10)
        plt.tick_params(labelleft=False, left=False)
        
    plt.plot(energy, totalDensity, color=colors[colorCounter], label="Total DOS")
    plt.fill_between(energy, totalDensity, color=colors[colorCounter], alpha=0.3)
    colorCounter = (colorCounter + 1)%len(colors)
    
    if(config["subplots"]):
        plt.legend()

#Gets the contribution of each orbital
for i in range (config["NEDOS"] + 6, len(doscar)):
    doscar[i] = doscar[i].split()
    j = (i - (config["NEDOS"] + 6))%(config["NEDOS"] + 1)
    if(j != 0):
        sDensity[j - 1] = sDensity[j - 1] + float(doscar[i][1])
        pDensity[j - 1] = pDensity[j - 1] + float(doscar[i][2]) + float(doscar[i][3]) + float(doscar[i][4])
        dDensity[j - 1] = dDensity[j - 1] + float(doscar[i][5]) + float(doscar[i][6]) + float(doscar[i][7]) + float(doscar[i][8]) + float(doscar[i][9])
        if(config["elements"]):
            totalElementDensity[j - 1] = sDensity[j - 1] + pDensity[j - 1] + dDensity[j - 1]
    elif(j == 0 and config["elements"]):
        if(atomCounter%numOfElements[elementCounter] == 0 and atomCounter != 0):
            #Plots the orbitals of each element
            if(config["orbitals"] and config["subplots"] == False):
                # Plot and fill each curve
                label1 = elementNames[elementCounter] + "(s)"
                plt.plot(energy, sDensity, color=colors[colorCounter], label=label1)
                plt.fill_between(energy, sDensity, color=colors[colorCounter], alpha=0.3)
                colorCounter = (colorCounter + 1)%len(colors)
                
                label2 = elementNames[elementCounter] + "(p)"
                plt.plot(energy, pDensity, color=colors[colorCounter], label=label2)
                plt.fill_between(energy, pDensity, color=colors[colorCounter], alpha=0.3)
                colorCounter = (colorCounter + 1)%len(colors)
                
                label3 = elementNames[elementCounter] + "(d)"
                plt.plot(energy, dDensity, color=colors[colorCounter], label=label3)
                plt.fill_between(energy, dDensity, color=colors[colorCounter], alpha=0.3)
                colorCounter = (colorCounter + 1)%len(colors)
            #Plots the elements on seperate subplots w/ or w/o orbitals
            elif(config["subplots"]):
                plt.subplot(len(elementNames) + 1, 1, elementCounter + 2)
                plt.ylim(0, ymax)
                plt.xlim(energy[0], energy[-1])
                ax = plt.gca()
                plt.xlabel("Energy (eV)", fontsize=fontSize)
                plt.ylabel("Density of States", fontsize=fontSize, labelpad=10)
                plt.tick_params(labelleft=False, left=False)
                plt.plot(energy, totalElementDensity, color=colors[colorCounter], label=elementNames[elementCounter])
                plt.fill_between(energy, totalElementDensity, color=colors[colorCounter], alpha=0.3)
                colorCounter = (colorCounter + 1)%len(colors)
                
                #Plots the individual orbitals for each element
                if(config["orbitals"]):
                    # Plot and fill each curve
                    label1 = elementNames[elementCounter] + "(s)"
                    plt.plot(energy, sDensity, color=colors[colorCounter], label=label1)
                    plt.fill_between(energy, sDensity, color=colors[colorCounter], alpha=0.3)
                    colorCounter = (colorCounter + 1)%len(colors)
                    
                    label2 = elementNames[elementCounter] + "(p)"
                    plt.plot(energy, pDensity, color=colors[colorCounter], label=label2)
                    plt.fill_between(energy, pDensity, color=colors[colorCounter], alpha=0.3)
                    colorCounter = (colorCounter + 1)%len(colors)
                    
                    label3 = elementNames[elementCounter] + "(d)"
                    plt.plot(energy, dDensity, color=colors[colorCounter], label=label3)
                    plt.fill_between(energy, dDensity, color=colors[colorCounter], alpha=0.3)
                    colorCounter = (colorCounter + 1)%len(colors)
                
                plt.legend()
            
            #Plots the total dos for each element
            else:
                plt.plot(energy, totalElementDensity, color=colors[colorCounter], label=elementNames[elementCounter])
                plt.fill_between(energy, totalElementDensity, color=colors[colorCounter], alpha=0.3)
                colorCounter = (colorCounter + 1)%len(colors)
            elementCounter = elementCounter + 1
            atomCounter = 0
            totalElementDensity = [0] * config["NEDOS"]
            sDensity = [0] * config["NEDOS"]
            pDensity = [0] * config["NEDOS"]
            dDensity = [0] * config["NEDOS"]
        atomCounter = atomCounter + 1

if(config["elements"]):
    #Plots the orbitals of each element
    if(config["orbitals"] and config["subplots"] == False):
        # Plot and fill each curve
        label1 = elementNames[elementCounter] + "(s)"
        plt.plot(energy, sDensity, color=colors[colorCounter], label=label1)
        plt.fill_between(energy, sDensity, color=colors[colorCounter], alpha=0.3)
        colorCounter = (colorCounter + 1)%len(colors)
        
        label2 = elementNames[elementCounter] + "(p)"
        plt.plot(energy, pDensity, color=colors[colorCounter], label=label2)
        plt.fill_between(energy, pDensity, color=colors[colorCounter], alpha=0.3)
        colorCounter = (colorCounter + 1)%len(colors)
        
        label3 = elementNames[elementCounter] + "(d)"
        plt.plot(energy, dDensity, color=colors[colorCounter], label=label3)
        plt.fill_between(energy, dDensity, color=colors[colorCounter], alpha=0.3)
        colorCounter = (colorCounter + 1)%len(colors)
    #Plots the elements on seperate subplots w/ or w/o orbitals
    elif(config["subplots"]):
        plt.subplot(len(elementNames) + 1, 1, elementCounter + 2)
        plt.ylim(0, ymax)
        plt.xlim(energy[0], energy[-1])
        ax = plt.gca()
        plt.xlabel("Energy (eV)", fontsize=fontSize)
        plt.ylabel("Density of States", fontsize=fontSize, labelpad=10)
        plt.tick_params(labelleft=False, left=False)
        plt.plot(energy, totalElementDensity, color=colors[colorCounter], label=elementNames[elementCounter])
        plt.fill_between(energy, totalElementDensity, color=colors[colorCounter], alpha=0.3)
        colorCounter = (colorCounter + 1)%len(colors)
        
        #Plots the individual orbitals for each element
        if(config["orbitals"]):
            # Plot and fill each curve
            label1 = elementNames[elementCounter] + "(s)"
            plt.plot(energy, sDensity, color=colors[colorCounter], label=label1)
            plt.fill_between(energy, sDensity, color=colors[colorCounter], alpha=0.3)
            colorCounter = (colorCounter + 1)%len(colors)
            
            label2 = elementNames[elementCounter] + "(p)"
            plt.plot(energy, pDensity, color=colors[colorCounter], label=label2)
            plt.fill_between(energy, pDensity, color=colors[colorCounter], alpha=0.3)
            colorCounter = (colorCounter + 1)%len(colors)
            
            label3 = elementNames[elementCounter] + "(d)"
            plt.plot(energy, dDensity, color=colors[colorCounter], label=label3)
            plt.fill_between(energy, dDensity, color=colors[colorCounter], alpha=0.3)
            colorCounter = (colorCounter + 1)%len(colors)
    #Plots the total dos for each element
    else:
        plt.plot(energy, totalElementDensity, color=colors[colorCounter], label=elementNames[elementCounter])
        plt.fill_between(energy, totalElementDensity, color=colors[colorCounter], alpha=0.3)
        colorCounter = (colorCounter + 1)%len(colors)
    elementCounter = elementCounter + 1

# Plot and fill each curve
if(config["orbitals"] and config["elements"] != True):
    plt.plot(energy, sDensity, color=colors[colorCounter], label="s")
    plt.fill_between(energy, sDensity, color=colors[colorCounter], alpha=0.3)
    colorCounter = (colorCounter + 1)%len(colors)
    
    plt.plot(energy, pDensity, color=colors[colorCounter], label="p")
    plt.fill_between(energy, pDensity, color=colors[colorCounter], alpha=0.3)
    colorCounter = (colorCounter + 1)%len(colors)
    
    plt.plot(energy, dDensity, color=colors[colorCounter], label="d")
    plt.fill_between(energy, dDensity, color=colors[colorCounter], alpha=0.3)
    colorCounter = (colorCounter + 1)%len(colors)

plt.ylim(0)
plt.xlim(energy[0], energy[-1])
ax = plt.gca()
plt.xlabel("Energy (eV)", fontsize=fontSize)
plt.ylabel("Density of States", fontsize=fontSize, labelpad=10)
plt.tick_params(labelleft=False, left=False)
plt.legend()
plt.savefig(config["saveloc"])
plt.show()    
