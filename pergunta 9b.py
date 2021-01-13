import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import contacts
import math

from numpy.core import numeric
from numpy.core.defchararray import zfill


universe = mda.Universe("/Users/bernardoaugusto/Desktop/3º ano/1º semestre/Computational Biochemistry/quiz_bernardo/files/md_300.gro", 
                        "/Users/bernardoaugusto/Desktop/3º ano/1º semestre/Computational Biochemistry/quiz_bernardo/files/md_300.xtc")

res_cys_sg = universe.select_atoms("(resname CYS and name SG)")
res_cys_hg = universe.select_atoms("(resname CYS and name HG)")
res_glu = universe.select_atoms("(resname GLU and name OE*)")

sg_pos = []
hg_pos = []
glu_pos = []
for x in range(len(res_cys_sg)):
    sg_pos.append(res_cys_sg[x].position)

for x in range(len(res_cys_hg)):
    hg_pos.append(res_cys_hg[x].position)

for z in range(len(res_glu)):
    glu_pos.append(res_glu[z].position)

vetorDH_list = []
vetorAH_list = []
for x in range(len(sg_pos)):
    for x in range(len(hg_pos)):
        for z in range(len(glu_pos)):
            vetorDH = sg_pos[x] - hg_pos[x]
            vetorAH = hg_pos[x] - glu_pos[z]
            cos = np.dot(vetorDH, vetorAH) / (np.linalg.norm(vetorDH) * np.linalg.norm(vetorAH))
            alfa = math.acos(cos)
            vetorAH_list.append(vetorAH)
            vetorDH_list.append(vetorDH)

            print(cos, alfa)