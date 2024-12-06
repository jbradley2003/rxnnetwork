import numpy as np
import xyz2mol
from rdkit import Chem
from rdkit.Chem import Draw


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


def print_arrays(array_list):
    """prints out the numpy arrays so we can see them better

    Args:
        array_list (numpy array): list of numpy arrays
    """
    for i in array_list:
        print(i)
        print()

# ---------------- CONVERSION MATRICES TESTING ----------------------------------------------------
# ------------------------------------------------------------------------------------------
r = np.array([[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,0,1,0,0], [0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]])

test_matrix = np.array([[0,1,0],[0, 0, 1],[0, 0, 0]])

test2_matrix = np.empty((5,5))
result = one_cycle_conversion_matrices(test_matrix)



# --------------------  For all of our Tests, we will start with FORMALDEHYDE, CH2O, easy structure

# atoms list in the matrix: C, O, H, H
atoms = ['c', 'o', 'h', 'h']

# convert atoms to list of integers
atoms_int = []
for i in atoms:
    atoms_int.append(xyz2mol.int_atom(i))

#                            c  o  h  h
formaldehyde_AC = np.array([[0, 1, 1, 1],
                            [1, 0, 0, 0],
                            [1, 0, 0, 0],
                            [1, 0, 0, 0]])

# ---------------------- TESTING OUT THE BOND ORDER MATRIX CONVERSION ------------------------

# AC2BO function: params: AC: adjacency matrix
#                         atoms: int list of atoms
#                         charge: the charge of the structure
formaldehyde_BO = xyz2mol.AC2BO(formaldehyde_AC, atoms_int, 0)

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


def BO2Mol(BO, atoms):
    """ takes in a bond order matrix and converts it to a RDkit molecular graph

    Args:
        BO (numpy array): the bond order matrix that describes our structure
        atoms (list): the list of atomic numbers in our structure
    """
    rwMol = Chem.RWMol()
    
    BO = BO[0]
    # add the atoms
    atom_indices = []
    for i in atoms:
        atom = Chem.Atom(i)
        idx = rwMol.AddAtom(atom)
        atom_indices.append(idx)
        
    
    # add the bonds
    for i in range(len(BO)):
        for j in range(i+1, len(BO)):  # Only iterate upper triangular matrix
            bond_order = BO[i][j]
            if bond_order > 0:
                # Determine bond type based on bond order
                if bond_order == 1:
                    bond_type = Chem.BondType.SINGLE
                elif bond_order == 2:
                    bond_type = Chem.BondType.DOUBLE
                elif bond_order == 3:
                    bond_type = Chem.BondType.TRIPLE
                else:
                    raise ValueError(f"Unsupported bond order: {bond_order}")
                
                # Add the bond
                rwMol.AddBond(i, j, bond_type)
    
    mol = rwMol.GetMol()
    return mol


# form the RDkit molecular object of formaldehyde
formaldehyde = BO2Mol(formaldehyde_BO, atoms_int)

img = Draw.MolToImage(formaldehyde)
img.show()













