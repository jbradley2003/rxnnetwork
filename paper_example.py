import conversion as con
import numpy as np
import ac2mol
import networkx as nx
import matplotlib.pyplot as plt

# ----------- reactant from the paper ----------------
#              h h i cl i cl
r = np.array([[0,1,0,0,0,0],
              [1,0,0,0,0,0],
              [0,0,0,1,0,0], 
              [0,0,1,0,0,0],
              [0,0,0,0,0,1],
              [0,0,0,0,1,0]])

r_atoms = ['h','h','i','cl','i','cl']

# ------------- product from the paper ------------------
#              h h i cl i cl
p = np.array([[0,0,0,1,0,0],
              [0,0,0,0,0,1],
              [0,0,0,0,1,0], 
              [1,0,0,0,0,0],
              [0,0,1,0,0,0],
              [0,1,0,0,0,0]])

p_atoms = ['h','h','i','cl','i','cl']
p_atoms_int = ac2mol.int_atom(p_atoms)
p_smiles = ac2mol.AC2Smiles(p, p_atoms_int)


G, inters, smiles_list = con.generateNetwork(r, r_atoms, 6)


con.draw_network(G)
plt.savefig('paper_example.png')

# finding the product in our reaction network
for i in range(len(inters)):
    for j in inters[i]:
        if p_smiles == j:
            print("The prouct was found in Intermediate layer ", i+1)
            
            

