import conversion as con
import numpy as np
import ac2mol
import networkx as nx
import matplotlib.pyplot as plt

smiles_list = []
#              h h i cl i cl
r = np.array([[0,0,0,1,0,0],
              [0,0,0,0,0,1],
              [0,0,0,0,1,0], 
              [1,0,0,0,0,0],
              [0,0,1,0,0,0],
              [0,1,0,0,0,0]])

r_atoms = ['h','h','i','cl','i','cl']

r_atoms_int = ac2mol.int_atom(r_atoms)

# print(r_atoms_int)
# r_BO = ac2mol.AC2BO(r, r_atoms_int, 0)
# print(r_BO)


G = nx.Graph()
# r_smiles = ac2mol.AC2Smiles(r, r_atoms_int)
G.add_node("[H][H].[Cl][Cl].[I][I]", color = 'purple')
smiles_list.append("[H][H].[Cl][Cl].[I][I]")

con.update_graph(r, r_atoms_int, smiles_list, G, 'red')

inter_1 = con.create_intermediates(r, r_atoms_int)

inter_2 = []
for i in inter_1:
    ac = ac2mol.Mol2AC(i)
    smiles_list, G, new_smile = con.update_graph(ac[0], ac[1], smiles_list, G, 'green')
    inter_2 = inter_2 + new_smile
    
inter_3 = []
for i in inter_2:
    ac = ac2mol.Smiles2AC(i)
    print(ac[0])
    print(ac[1])
    con.update_graph(ac[0], ac[1], smiles_list, G, 'blue')

node_colors = [G.nodes[node]['color'] for node in G.nodes()] 
pos = nx.spring_layout(G, seed = 42)
nx.draw(G, pos, node_color = node_colors, node_size=20)
plt.savefig('paper_example.png')

print("Number of Cycles: ",len(nx.cycle_basis(G)))
print("Number of Nodes: ",len(G.nodes))