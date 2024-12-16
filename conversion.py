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
        for j in range(i + 1, num_atoms):  # Upper triangle of the matrix
            A_ij_list.append((i,j))

    
    permutation_size = len(A_ij_list)
    # 2 chemical changes
    permutation_list = [1]*max_bond_changes + [0]*(permutation_size-max_bond_changes)
    unique_permutations = permute_unique(permutation_list)
    print("permutation size: ", len(unique_permutations))
    # 1 chemical change
    for i in range(permutation_size):
        temp = permutation_size * [0]
        temp[i] = 1
        unique_permutations.append(temp)
    
    print("permutation list size: ", len(unique_permutations))
    
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

def unique_RDkit_molecules(molecule_list):
    """ takes in a list of RDKit molecules and returns the unique structures and their isomers

    Args:
        molecule_list (list): list of RDKit molecules

    Returns:
        list : list of RDKit molecules that are unique with their isomers
    """
    unique_structures = []
    unique_smiles = set()
    # loop through all of the molecules and get rid of the repeated ones
    for i in molecule_list:
        smiles = Chem.MolToSmiles(i, isomericSmiles = True)
        if smiles not in unique_smiles:
            unique_smiles.add(smiles)
            unique_structures.append(i)
    
    # for i in unique_structures:
    #     isomers = tuple(EnumerateStereoisomers(i))
    #     # loop through all of the isomers and add unique ones to our list
    #     for j in isomers:
    #         print(j)
    #         smiles = Chem.MolToSmiles(j, isomericSmiles = True)
    #         if smiles not in unique_smiles:
    #             unique_smiles.add(smiles)
    #             unique_structures.append(j)

    return unique_structures



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
        intermediate = ac2mol.AC2Mol(i, atoms)
        if intermediate != None: 
            # add the best intermediate
            intermediate_structures.append(intermediate)
    
    intermediate_structures = unique_RDkit_molecules(intermediate_structures)
    print(intermediate_structures)
    
    return intermediate_structures


# ------------------ updating the set and adding it to the reaction Graph --------------------------    

def update_graph(reactant, atoms, smile_list, G, col):
    """ updates our ac set for unique structures

    Args:
        mol_intermediates (list): list of RDKit molecules
        ac_set (set): set of unique AC matrices

    Returns:
        ac_set: the set of unique AC matrices
        G: updated graph
    """
    new_smiles = []
    mol_intermediates = create_intermediates(reactant, atoms)
    reactant_smiles = ac2mol.AC2Smiles(reactant, atoms)
    # update ac_set and nodes
    for i in mol_intermediates:
        mol_smiles = ac2mol.Mol2Smiles(i)
        if mol_smiles not in smile_list:
            new_smiles.append(mol_smiles)
            smile_list.append(mol_smiles)
            G.add_node(mol_smiles, color = col)
        # add edge
        if reactant_smiles != mol_smiles:
            G.add_edge(reactant_smiles, mol_smiles)
    
    return smile_list, G, new_smiles
    



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




