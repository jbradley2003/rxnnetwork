import numpy as np
import ac2mol
from rdkit.Chem import Draw 
import conversion as con


# ---------------- CONVERSION MATRICES TESTING ----------------------------------------------------
# ------------------------------------------------------------------------------------------
r = np.array([[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,0,1,0,0], [0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]])

test_matrix = np.array([[0,1,0],[0, 0, 1],[0, 0, 0]])

test2_matrix = np.empty((5,5))
result = con.one_cycle_conversion_matrices(test_matrix)



# --------------------  For all of our Tests, we will start with FORMALDEHYDE, CH2O, easy structure

# atoms list in the matrix: C, O, H, H
atoms = ['c', 'o', 'h', 'h']

# atoms list in integer form
atoms_int = ac2mol.int_atom(atoms)

#                            c  o  h  h
formaldehyde_AC = np.array([[0, 1, 1, 1],
                            [1, 0, 0, 0],
                            [1, 0, 0, 0],
                            [1, 0, 0, 0]])

# ---------------------- TESTING OUT THE BOND ORDER MATRIX CONVERSION ------------------------

# AC2BO function: params: AC: adjacency matrix
#                         atoms: int list of atoms
#                         charge: the charge of the structure
formaldehyde_BO = ac2mol.AC2BO(formaldehyde_AC, atoms_int, 0)


print(formaldehyde_BO[0])
print()
print(formaldehyde_BO[1])






# ---------------- TESTING OUT AC TO 3D MOLECULAR -----------------------------------------
# we are gonna call the AC2mol function:

# params: charge: the charge of the overall species? or a list of the charges?
#         atoms: a list of integers that represent the atoms
#         AC: the adjacency matrix
#         mol: a molecular object in RDKit, in their case, they combined it with 
#              the coordinate systm to get a more accurate mol. We are just gonna try
#              to use the default one that is constructed with a list of atoms





# form the RDkit molecular object of formaldehyde
formaldehyde = ac2mol.AC2Mol(formaldehyde_AC, atoms_int)





# ---------------------------- trying a structure that should have stereoisomers ------------------------

# C2F2H2

C2F2H2_AC = np.array([[0, 1, 1, 0, 1, 0],
                      [1, 0, 0, 1, 0, 1],
                      [1, 0, 0, 0, 0, 0], 
                      [0, 1, 0, 0, 0, 0], 
                      [1, 0, 0, 0, 0, 0], 
                      [0, 1, 0, 0, 0, 0]])

C2F2H2_atoms = ['c', 'c', 'f', 'f', 'h', 'h']

C2F2H2_atoms_int = ac2mol.int_atom(C2F2H2_atoms)

C2F2H2 = ac2mol.AC2Mol(C2F2H2_AC, C2F2H2_atoms_int)

best_isomer = C2F2H2[0]
other_isomers = C2F2H2[1]

img = Draw.MolToImage(best_isomer)
img.show()

# for i in other_isomers:
#     img = Draw.MolToImage(i)
#     img.show()

# ---------------- attempting to find all of the intermediates for one cycle of C2F2H2 --------------------------

C2F2H2_intermediates = con.create_intermediates(C2F2H2_AC, C2F2H2_atoms_int)
for i in C2F2H2_intermediates:
    con.draw_molecules(i)