"""This module is used to create all of the conversion matrices and houses the functions that generate
the reaction networks"""

import numpy as np
import ac2mol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
import networkx as nx
from rdkit import DataStructs


# --------------- conversion matrices and creating intermediates ---------------------------
def one_cycle_conversion_matrices(molecular_graph, active_indices, max_bond_changes=2):
    """
    Generate all conversion matrices for a 1 iteration of a given molecular graph.

    Args:
        molecular_graph (numpy matrix): Atom connectivity (AC) matrix of the starting molecular graph.
        max_bond_changes (int): Maximum number of bond changes allowed (formation/dissociation) per step.

    Returns:
        list: A list of valid conversion numpy matrices.
    """


    # Iterate over non-diagonal upper triangular matrix indices
    # store those indices in a list
    
    # TODO: change these to active atom matrices
    active_atom_matrix = create_active_atom_matrix(active_indices, molecular_graph)
    mol_graph_size = molecular_graph.shape[0]
    num_active_atoms = active_atom_matrix.shape[0]
    
    A_ij_list = []
    for i in range(num_active_atoms-1):
        for j in range(i + 1, num_active_atoms):  # Upper triangle of the matrix
            A_ij_list.append((i,j))


    permutation_size = len(A_ij_list)
    # 2 chemical changes
    permutation_list = [1]*max_bond_changes + [0]*(permutation_size-max_bond_changes)
    
    unique_permutations = permute_unique(permutation_list)
    # 1 chemical change
    for i in range(permutation_size):
        temp = permutation_size * [0]
        temp[i] = 1
        unique_permutations.append(temp)
    
    
    # form the permutation matrices
    # 1. for each permutation list
    # 2. match the index of the value in the permutation list with the indices in the adjacency list
    conversion_matrices = []
    for matrix in unique_permutations:
        # create an non-zero upper triangular copy of the active atoms graph
        upper = np.triu(np.copy(active_atom_matrix))

        for x in range(len(matrix)):
            i = A_ij_list[x][0]
            j = A_ij_list[x][1]
            
            # if in the reactant graph it is 1
            # and if in our permutation matrix it is 1
            # switch to -1
            if active_atom_matrix[i,j] == 1 and matrix[x] == 1:
                upper[i,j] = -1
            else: 
                upper[i,j] = matrix[x]
        
        # add in the lower symmetric component
        lower = upper.T
        conversion = upper + lower
        
        
        # add the permuted back into the proper spots in original conversion
        temp_AC = np.zeros((mol_graph_size,mol_graph_size))
        for i in range(num_active_atoms):
            for j in range(num_active_atoms):
                temp_AC[active_indices[i], active_indices[j]] = conversion[i,j]     
        
        conversion_matrices.append(temp_AC)
    
    return conversion_matrices


def permute_unique(nums):
    """ Creates all unique permutations from a list. from leet code :)

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

def create_active_atom_matrix(active_atom_indices, AC):
    """ create the active atom matrix

    Args:
        active_atom_indices (list): positions of active atoms in original adjacency matrix, zero indexing
        AC (np.array): adjacency matrix

    Returns:
        np.array: active atom matrix
    """
    active_atom_matrix = AC[np.ix_(active_atom_indices, active_atom_indices)]
    return active_atom_matrix


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



def create_intermediates(reactant, atoms, active_indices):
    """ create all of the possible intermediate molecular graphs that stem from one structure

    Args:
        reactant (numpy array): the 2d AC matrix of the reactant
        atoms (list): the interger list of atoms
    """
    # find all of the conversion matrices for that structure
    conversion_matrices = one_cycle_conversion_matrices(reactant, active_indices)
    # find all of the new AC matrices
    intermediate_ac = []
    for i in conversion_matrices:
        intermediate_ac.append(reactant + i)
        
    
    # create the molecular graphs of them
    intermediate_structures = []
    for i in intermediate_ac:
        intermediates = ac2mol.AC2Mol(i, atoms)
        if intermediates != None: 
            # add the best intermediate
            intermediate_structures.append(intermediates)
    
    intermediate_structures = unique_RDkit_molecules(intermediate_structures)
    
    return intermediate_structures


def tanimoto_difference(smile1, smile2):
    """ calculates the tanimoto distance between two structures

    Args:
        smile1 (string): first SMILES string
        smile2 (string): second SMILES string

    Returns:
        _type_: tanimoto difference 
    """
    fpgen = AllChem.GetRDKitFPGenerator()
    mol1 = ac2mol.MolFromSmiles(smile1, sanitize = False)
    mol2 = ac2mol.MolFromSmiles(smile2, sanitize = False)
    # default fingerprint, can change later
    fingerprint1 = fpgen.GetFingerprint(mol1)
    fingerprint2 = fpgen.GetFingerprint(mol2)
    
    similarity = DataStructs.TanimotoSimilarity(fingerprint1,fingerprint2)
    distance = 1 - similarity
    return distance

def chemical_distance(smile1, smile2):
    """ calculates the chemical distance between two structures

    Args:
        smile1 (string): first SMILES string
        smile2 (string): second SMILES string

    Returns:
        num: the chemical distance between two structures
    """
    ac1 = np.array(ac2mol.Smiles2AC(smile1)[0])
    ac2 = np.array(ac2mol.Smiles2AC(smile2)[0])


    return np.sum(np.abs(ac1 - ac2))
    


# ------------------ updating the reaction network --------------------------    

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
            G.add_edge(reactant_smiles, mol_smiles, weight = tanimoto_difference(reactant_smiles, mol_smiles))
    
    return smile_list, G, new_smiles



def update_graph_with_distance(reactant, atoms, smile_list, G, col):
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
            G.add_edge(reactant_smiles, mol_smiles, weight = chemical_distance(reactant_smiles, mol_smiles))
    
    return smile_list, G, new_smiles


    
def generateNetwork(ac, ac_list, n):
    """ master equation that generates the graphs

    Args:
        ac (numpy array): the atomic connectivity reactant matrix
        ac_list (list): list of atoms string form
        n (int): the number of intermediate layers

    Returns:
        G: the final reaction network
        inters: all of the structures that are added at each layer
        smiles_list: the list of all of the smiles strucutures
    """
    # color list
    color_list = ["#CC2F00", '#DB6600', '#E39E00', '#76B80D', '#007668', '#006486', '#465AB2', '#6D47B1', '#873B9C', 'brown']
    G = nx.Graph()
    # initiate the graph with reactant
    smiles_list = []
    r_int_list = ac2mol.int_atom(ac_list)
    r_smiles = ac2mol.AC2Smiles(ac, r_int_list)
    smiles_list.append(r_smiles)
    G.add_node(r_smiles, color='purple')
    
    # list that will be used to store new intermediates after each iteration
    inters = [[r_smiles]]
    # n iterations
    for i in range(n):
        new_inter = []
        for j in inters[i]:
            ac_j = ac2mol.Smiles2AC(j)
            _, _, new_smile = update_graph_with_distance(ac_j[0], ac_j[1], smiles_list, G, color_list[i])
            new_inter = new_inter + new_smile
        inters.append(new_inter)
 
        
    return G, inters, smiles_list
    
 
# ------------- display structures and network ---------------------------   
def draw_network(G, node_size = 20):
    """ draws out our graph in NetworkX

    Args:
        G (networkx graph): the networkx graph that we want to display
        node_size (int, optional): the size of our nodes in the graph
    """
    node_colors = [G.nodes[node]['color'] for node in G.nodes()] 
    # pos = nx.spring_layout(G, seed = 42)
    nx.draw(G, node_color = node_colors, node_size = node_size)

    
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





