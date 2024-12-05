import numpy as np

def generate_conversion_matrices(reactant_graph, max_bond_changes=2):
    """
    Generate conversion matrices for a given reactant molecular graph.

    Args:
        reactant_graph (matrix): Atom connectivity (AC) matrix of the reactant.
        max_bond_changes (int): Maximum number of bond changes (formation/dissociation) per step.

    Returns:
        list: A list of valid conversion matrices.
    """

    num_atoms = len(reactant_graph)
    conversion_matrices = []

    # Iterate over all possible changes in the AC matrix
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):  # Upper triangle of the matrix
            current_bond = reactant_graph[i][j]
            # Case 1: Bond formation (0 -> 1)
            if current_bond == 0:
                new_matrix1 = np.copy(reactant_graph)
                new_matrix2 = np.copy(reactant_graph)
                
                new_matrix1[i][j] = 1  # Form bond
                new_matrix1[j][i] = 1  # Symmetric update
                if is_valid_matrix(reactant_graph, new_matrix1, max_bond_changes):
                    conversion_matrices.append(new_matrix1)
                    
                new_matrix2[i][j] = 0  # Form bond
                new_matrix2[j][i] = 0  # Symmetric update
                if is_valid_matrix(reactant_graph, new_matrix2, max_bond_changes):
                    conversion_matrices.append(new_matrix2)
            # Case 2: Bond cleavage (1 -> -1)
            elif current_bond == 1:
                new_matrix1 = np.copy(reactant_graph)
                new_matrix1[i][j] = -1  # Break bond
                new_matrix1[j][i] = -1  # Symmetric update
                if is_valid_matrix(reactant_graph, new_matrix1, max_bond_changes):
                    conversion_matrices.append(new_matrix1)
    return conversion_matrices


# testing another version
def generate_conversion_matrices_v1(reactant_graph, max_bond_changes=2):
    """
    Generate conversion matrices for a given reactant molecular graph.

    Args:
        reactant_graph (numpy matrix): Atom connectivity (AC) matrix of the reactant.
        max_bond_changes (int): Maximum number of bond changes allowed (formation/dissociation) per step.

    Returns:
        list: A list of valid conversion matrices.
    """

    num_atoms = reactant_graph.shape[0]
    print(num_atoms)

    # Iterate over non-diagonal upper triangular matrix indices
    # store those indices in a list
    A_ij_list = []
    for i in range(num_atoms-1):
        # find the unique matrices that only have 1 change
        for j in range(i + 1, num_atoms):  # Upper triangle of the matrix
            A_ij_list.append((i,j))

    print(A_ij_list)
    
    permutation_size = len(A_ij_list)
    permutation_list = [1]*max_bond_changes + [0]*(permutation_size-max_bond_changes)
    unique_permutations = permuteUnique(permutation_list)
    
    # form the permutation matrices
    # 1. for each permutation list
    # 2. match the index of the value in the permutation list with the indices in the adjacency list
    conversion_matrices = []
    for matrix in unique_permutations:
        temp = np.copy(reactant_graph)
        for x in range(len(matrix)):
            i = A_ij_list[x][0]
            j = A_ij_list[x][1]
            # for each element in each permutation matrix
            temp[i,j] = matrix[x]
        conversion_matrices.append(temp)
    
    return conversion_matrices




def permuteUnique(nums):
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





def is_valid_matrix(original_matrix, new_matrix, max_bond_changes):
    """
    Check if a matrix is valid based on maximum bond changes.

    Args:
        original_matrix (matrix): The original AC matrix.
        new_matrix (matrix): The AC matrix after changes.
        max_bond_changes (int): Maximum allowable changes.

    Returns:
        bool: True if valid, False otherwise.
    """
    # Count the number of changes (formation: 0->1, cleavage: 1->-1)
    changes = np.sum(np.abs(original_matrix - new_matrix) > 0)
    return changes <= max_bond_changes



r = [[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,0,1,0,0], [0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]]

test_matrix = np.array([[0,0,0],[0, 0, 0],[0, 0, 0]])
test2_matrix = np.empty((20,20))
print(generate_conversion_matrices_v1(test_matrix))

# test = np.copy(r)


# l = generate_conversion_matrices(r)

# m = np.unique(np.array(l),axis=0)
# i = np.array([]) 

# for val in m:
#     inter = np.add(test, val)
#     print(inter)








