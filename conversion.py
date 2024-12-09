import numpy as np
import ac2mol
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers


def one_cycle_conversion_matrices(molecular_graph, max_bond_changes=2):
    """
    Generate all conversion matrices for a 1 iteration of a given molecular graph.

    Args:
        molecular_graph (numpy matrix): Atom connectivity (AC) matrix of the starting molecular graph.
        max_bond_changes (int): Maximum number of bond changes allowed (formation/dissociation) per step.

    Returns:
        list: A list of valid conversion numpy matrices.
    """

    num_atoms = molecular_graph.shape[0]

    # Iterate over non-diagonal upper triangular matrix indices
    # store those indices in a list
    A_ij_list = []
    for i in range(num_atoms-1):
        # find the unique matrices that only have 1 change
        for j in range(i + 1, num_atoms):  # Upper triangle of the matrix
            A_ij_list.append((i,j))

    
    permutation_size = len(A_ij_list)
    permutation_list = [1]*max_bond_changes + [0]*(permutation_size-max_bond_changes)
    unique_permutations = permute_unique(permutation_list)
    
    # form the permutation matrices
    # 1. for each permutation list
    # 2. match the index of the value in the permutation list with the indices in the adjacency list
    conversion_matrices = []
    for matrix in unique_permutations:
        # create an non-zero upper triangular copy of the reactant graph
        upper = np.triu(np.copy(molecular_graph))
        for x in range(len(matrix)):
            i = A_ij_list[x][0]
            j = A_ij_list[x][1]
            
            # if in the reactant graph it is 1
            # and if in our permutation matrix it is 1
            # switch to -1
            if molecular_graph[i,j] == 1 and matrix[x] == 1:
                upper[i,j] = -1
            else: 
                upper[i,j] = matrix[x]
        
        # add in the lower symmetric component
        lower = upper.T
        conversion = upper + lower
        conversion_matrices.append(conversion)
    
    return conversion_matrices


def permute_unique(nums):
    """ Creates all unique permutations from a list

    Args:
        nums (list): list of values that we want to create permutations for

    Returns:
        list of list: list of the permutation lists
    """
    res = []
    n = len(nums)
    visited = [False] * n
    def dfs(path, i):
        if i == n:
            res.append(path[:])
            return
        for j in range(n):
            if visited[j]: continue
            if j > 0 and nums[j] == nums[j-1] and not visited[j-1]:
                continue
            path.append(nums[j])
            visited[j] = True
            dfs(path, i+1)
            path.pop()
            visited[j] = False
    
    dfs([], 0)
    return res


def create_intermediates(reactant, atoms):
    """ create all of the possible intermediate molecular graphs that stem from one structure

    Args:
        reactant (numpy array): the 2d AC matrix of the reactant
        atoms (list): the interger list of atoms
    """
    # find all of the conversion matrices for that structure
    conversion_matrices = one_cycle_conversion_matrices(reactant)
    
    # find all of the new AC matrices
    intermediate_ac = []
    for i in conversion_matrices:
        intermediate_ac.append(reactant + i)
    
    # create the molecular graphs of them
    intermediate_structures = []
    for i in intermediate_ac:
        intermediates = ac2mol.AC2Mol(i, atoms)
        if intermediates != None: 
            intermediate_structures.append(intermediates[0])
    
    return intermediate_structures
    
    
def print_arrays(array_list):
    """prints out the numpy arrays so we can see them better

    Args:
        array_list (numpy array): list of numpy arrays
    """
    for i in array_list:
        print(i)
        print()
        
def draw_molecules(molecule):
    """print out the molecules

    Args:
        molecule (RDKit molecule object): it is a RDkit molecular object
    """
    img = Draw.MolToImage(molecule)
    img.show()

# ---------------- CONVERSION MATRICES TESTING ----------------------------------------------------
# ------------------------------------------------------------------------------------------
r = np.array([[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,0,1,0,0], [0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]])

test_matrix = np.array([[0,1,0],[0, 0, 1],[0, 0, 0]])

test2_matrix = np.empty((5,5))
result = one_cycle_conversion_matrices(test_matrix)



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

C2F2H2_intermediates = create_intermediates(C2F2H2_AC, C2F2H2_atoms_int)
for i in C2F2H2_intermediates:
    draw_molecules(i)