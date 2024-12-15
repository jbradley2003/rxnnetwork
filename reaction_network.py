import conversion as con
import numpy as np
import ac2mol
import networkx as nx
import matplotlib.pyplot as plt

smiles_list = []

C2F2H2_AC = np.array([[0, 1, 1, 0, 1, 0],
                      [1, 0, 0, 1, 0, 1],
                      [1, 0, 0, 0, 0, 0], 
                      [0, 1, 0, 0, 0, 0], 
                      [1, 0, 0, 0, 0, 0], 
                      [0, 1, 0, 0, 0, 0]])

C2F2H2_atoms = ['c', 'c', 'f', 'f', 'h', 'h']

C2F2H2_atoms_int = ac2mol.int_atom(C2F2H2_atoms)


G = nx.Graph()
C2F2H2_smiles = ac2mol.AC2Smiles(C2F2H2_AC, C2F2H2_atoms_int)
G.add_node(C2F2H2_smiles)
smiles_list.append(C2F2H2_smiles)

smile_list, G = con.update_graph(C2F2H2_AC, C2F2H2_atoms_int, smiles_list, G)

nx.draw(G)
plt.savefig('first_iter.png')

# first round of intermediates
inter_1 = con.create_intermediates(C2F2H2_AC, C2F2H2_atoms_int)

for i in inter_1:
    ac = ac2mol.Mol2AC(i)
    smiles_list, G = con.update_graph(ac[0], ac[1], smiles_list, G)

nx.draw(G)
plt.savefig('second_iter.png')