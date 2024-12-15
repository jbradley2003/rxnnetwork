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
G.add_node(C2F2H2_smiles, color = 'purple')
smiles_list.append(C2F2H2_smiles)

con.update_graph(C2F2H2_AC, C2F2H2_atoms_int, smiles_list, G, 'red')


# first round of intermediates
inter_1 = con.create_intermediates(C2F2H2_AC, C2F2H2_atoms_int)

inter_2 = []
for i in inter_1:
    ac = ac2mol.Mol2AC(i)
    smiles_list, G, new_smile = con.update_graph(ac[0], ac[1], smiles_list, G, 'green')
    inter_2 = inter_2 + new_smile
    
# pos = nx.spring_layout(G, seed = 42)
# nx.draw(G, pos, node_size=20, node_color='green')
# plt.savefig('second_iter.png')


# is it failing here
print("Does it get here?")
inter_3 = []
for i in inter_2:
    ac = ac2mol.Smiles2AC(i)
    print(ac[0])
    print(ac[1])
    con.update_graph(ac[0], ac[1], smiles_list, G, 'blue')

node_colors = [G.nodes[node]['color'] for node in G.nodes()] 
pos = nx.spring_layout(G, seed = 42)
nx.draw(G, pos, node_color = node_colors, node_size=20)
plt.savefig('third_iter.png')