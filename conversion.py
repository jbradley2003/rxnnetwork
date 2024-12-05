import numpy as np


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

# ----------------------------- TESTING ----------------------------------------------------
# ------------------------------------------------------------------------------------------
r = np.array([[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,0,1,0,0], [0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]])

test_matrix = np.array([[0,1,0],[0, 0, 1],[0, 0, 0]])

test2_matrix = np.empty((5,5))
result = one_cycle_conversion_matrices(test_matrix)

print(test_matrix)
print_arrays(result)










