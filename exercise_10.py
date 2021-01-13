# imports
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from numpy.linalg.linalg import norm
universe = mda.Universe('/Users/bernardoaugusto/Desktop/3º ano/1º semestre/Computational Biochemistry/quiz_bernardo/files/md_300.gro', 
                        '/Users/bernardoaugusto/Desktop/3º ano/1º semestre/Computational Biochemistry/quiz_bernardo/files/md_300.xtc')
import numpy
from numpy import linalg
import math 


# Select the the phe resids
res_phe = universe.select_atoms("(resname PHE and name CD*) or (resname PHE and name CE*) or\
(resname PHE and name CZ) or (resname PHE and name CA)")



vetor_list = []
centroide = []


# Generate the centroids and the vetors necessary for the angles
for i in range(0,len(res_phe),6):
    piCoords1=[res_phe[i].position,res_phe[i + 1].position,
    res_phe[i + 2].position,res_phe[i + 3].position,
    res_phe[i + 4].position,res_phe[i + 5].position]
    media = (res_phe[i].position + res_phe[i + 1].position + 
    res_phe[i + 2].position + res_phe[i + 3].position + 
    res_phe[i + 4].position + res_phe[i + 5].position) / 6
    centroide.append(media)

    vector1=piCoords1[0]-piCoords1[5]
    vector2=piCoords1[1]-piCoords1[4]
    vector3=piCoords1[2]-piCoords1[3]
    prod1=numpy.cross(vector1,vector2)
    prod2=numpy.cross(vector1,vector3)
    prod3=numpy.cross(vector2,vector3)
    normalVector1=(prod1+prod2+prod3)/3
    vetor_list.append(normalVector1)


# Calculates the distance and the angle and then classifies the angle (if exists)
for i in range(len(centroide)):
    for j in range(len(centroide)):
        if j != i:
            dist = numpy.linalg.norm(centroide[i]-centroide[j])
            if dist < 5.0:
                cos = numpy.dot(vetor_list[i], vetor_list[j]) / (numpy.linalg.norm(vetor_list[i]) * numpy.linalg.norm(vetor_list[j]))
                alfa = math.acos(cos)
                if 0 <= alfa <= (math.pi/6): 
                    print("The angle between is parellel")
                elif  (math.pi/6) < alfa <= (math.pi/3):
                    print("The angle is oblique")
                elif (math.pi/3) < alfa <= (math.pi/2):
                    print("The angle is perpendicular")
        else:
            pass
                    