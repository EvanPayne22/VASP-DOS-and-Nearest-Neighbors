# VASP DOS and Nearest Neighbors Histogram
This repository is supposed to be helpful for displaying the Density of States of a system from the VASP DOSCAR and POSCAR. There is additional script that will determine the distances between the nearest neighbors from a given POSCAR.

## doscar_plotter.py
### Description
This script produces plots of the density of states of a system. There are tags to display the orbitals, the elements, and just the total density of states. There are additional visualization features to display the orbitals of each element on seperate subplots. If you are using the element tag a POSCAR will be needed as well as the DOSCAR.

### Sample Plots
#### Orbitals
<img src="images/orbitals.png" alt="Displays the DOSCAR with orbital contributions" width="600"/>

#### Elements
<img src="images/orbitals.png" alt="Displays the DOSCAR with element contributions" width="600"/>

#### Subplots
<img src="images/subplots.png" alt="Displays the DOSCAR in subplots" width="600"/>

## nearest_neighbor_calculator.py
### Description 
A code used to generate nearest neighbor graphs between the same or different elements in a supercell. Useful for understanding the bond lengths inside of the system. 

