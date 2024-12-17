import conversion as con
import numpy as np
import ac2mol
import networkx as nx
import matplotlib.pyplot as plt
import time

# ----------- reactant from the paper ----------------
#              h h i cl i cl
r = np.array([[0,1,0,0,0,0],
              [1,0,0,0,0,0],
              [0,0,0,1,0,0], 
              [0,0,1,0,0,0],
              [0,0,0,0,0,1],
              [0,0,0,0,1,0]])

r_atoms = ['h','h','i','cl','i','cl']
r_atoms_int = ac2mol.int_atom(r_atoms)
r_smiles = ac2mol.AC2Smiles(r, r_atoms_int)

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


start_time = time.time()
G, inters, smiles_list = con.generateNetwork(r, r_atoms, 6)
end_time = time.time()


con.draw_network(G)
plt.savefig('paper_example.png')

# finding the product in our reaction network
for i in range(len(inters)):
    for j in inters[i]:
        if p_smiles == j:
            print("The product was found in Intermediate layer ", i+1)
            
            
            
# Convert edges to a DataFrame
# edges_df = nx.to_pandas_edgelist(G)
# print(edges_df)

# export as GEXF file for gephi
# nx.write_gexf(G, "paper_example.gexf")

# dijkstras algorithm
# distance, path = nx.single_source_dijkstra(G, source=r_smiles, target=p_smiles)
# print(f"Shortest Path from reactant to product: {path} with distance {distance}")
print("It took ", end_time - start_time," seconds to create the example reaction network.")