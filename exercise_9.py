# imports
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import contacts

# ---STROPMI----#
# ----FUNCTIONS----#
def plotter(distances, pos):
 import os
 os.makedirs("distances", exist_ok=True) # Creates folder if it doesn't exist already
 import matplotlib.pyplot as plt
 ax = plt.subplot()
 ax.plot(distances[:, 0], distances[:, 1], 'r--', lw=2, label=r"$R_G$")
 ax.set_xlabel("time (ps)")
 ax.set_ylabel(r"Distance in $\AA$")
 title = "distances/" + acids[pos].resname + str(acids[pos].resid) + " (" + acids[pos].name + ")" \
 + '-' + basics[pos].resname + str(basics[pos].resid) + " (" + basics[pos].name + ")"
 ax.figure.savefig(title)
 plt.draw() # draws graph in the file
 plt.cla() # clears the axis for next use
# ---- SNIOTCNUF----#

if __name__ == '__main__':
    MAX_DISTANCE = 5.0
    universe = mda.Universe("/Users/bernardoaugusto/Desktop/3º ano/1º semestre/Computational Biochemistry/quiz_bernardo/files/md_300.gro", 
                            "/Users/bernardoaugusto/Desktop/3º ano/1º semestre/Computational Biochemistry/quiz_bernardo/files/md_300.xtc") 
                            # creates the universe with the gro and xtc files

    # selects all the acids
    acid = universe.select_atoms("(resname LYS and name NZ) or (resname ARG and name NH*) or (resname ARG and name NZ) or (resname LYS and name NH*)")
    # selects all the basics
    basic = universe.select_atoms("(resname GLU and name OE*) or (resname ASP and name OE*)")
    # create a matrix with all distances between all positions
    distance = mda.analysis.distances.distance_array(acid.positions, basic.positions)
    lis = [];pos = [[], []]
    
    # filters the matrix to keep only those with max_distace or less
    for i, lin in enumerate(distance):
        for j, col in enumerate(lin):
            if col < MAX_DISTANCE and col != 0.0:
                pos[0].append(i) # acids
                pos[1].append(j) # basics
                print(acid[i].resname, acid[i].resid, acid[i].name, '-', basic[j].resname, basic[j].resid, basic[j].name, 'd=', round(col, 2))
                lis.append(round(col, 2))
                # picks up all the residues that had those distances
                # pos => (<pos of acid res>, <pos of basic res>) in the original AtomGroups
                indexes_a = pos[0]
                indexes_b = pos[1]
                acids = acid[indexes_a]
                basics = basic[indexes_b]
                global temp
                temp = []
                
    # calculates the distance over time
    for i in range(len(acids)):
        distances = []
        for ts in universe.trajectory:
            distances.append((universe.trajectory.time, np.linalg.norm(acids[i].position - basics[i].position)))
        distances = np.array(distances)
        plotter(distances, i) # plots for every bond

        