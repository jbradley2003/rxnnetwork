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

test = np.copy(r)

l = generate_conversion_matrices(r)

m = np.unique(np.array(l),axis=0)
i = np.array([]) 

for val in m:
    inter = np.add(test, val)
    print(inter)








