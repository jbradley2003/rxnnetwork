# simulating the reaction for Kuwata research
import time
import conversion as con
import numpy as np
import ac2mol
import networkx as nx
import matplotlib.pyplot as plt
from rdkit.Chem.rdmolfiles import MolFromSmiles

# ---------------- ethene + ozone -----------------

r = np.array([[0, 1, 1, 1, 0, 0, 0, 0, 0],
              [1, 0, 0, 0, 1, 1, 0, 0, 0],
              [1, 0, 0, 0, 0, 0, 0, 0, 0],
              [1, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 1, 0],
              [0, 0, 0, 0, 0, 0, 1, 0, 1],
              [0, 0, 0, 0, 0, 0, 0, 1, 0]])


r_atoms = ['c','c','h','h','h','h','o','o','o']
r_atoms_int = ac2mol.int_atom(r_atoms)
r_smiles = ac2mol.AC2Smiles(r, r_atoms_int)
reactant = MolFromSmiles(r_smiles, sanitize = False)
#con.draw_molecules(reactant)

# ---------- Criegee Intermediates -------------------
p = np.array([[0, 1, 1, 1, 0, 0, 0, 0, 0],
              [1, 0, 0, 0, 0, 0, 0, 0, 0],
              [1, 0, 0, 0, 0, 0, 0, 0, 0],
              [1, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 1, 0, 1, 1],
              [0, 0, 0, 0, 1, 0, 1, 0, 0],
              [0, 0, 0, 0, 0, 1, 0, 0, 0],
              [0, 0, 0, 0, 1, 0, 0, 0, 0],
              [0, 0, 0, 0, 1, 0, 0, 0, 0]])
p_atoms = ['c','o','h','h','c','o','o','h','h']
p_atoms_int = ac2mol.int_atom(p_atoms)
p_smiles = ac2mol.AC2Smiles(p, p_atoms_int)
criegee = MolFromSmiles(p_smiles, sanitize = False)
#con.draw_molecules(criegee)

start_time = time.time()
G, inters, smiles_list = con.generateNetwork(r, r_atoms, 6)
end_time = time.time()

# finding the product in our reaction network
for i in range(len(inters)):
    for j in inters[i]:
        if p_smiles == j:
            print("The product was found in Intermediate layer ", i+1)


# export as GEXF file for gephi
nx.write_gexf(G, "ozonolysis.gexf")

# dijkstras algorithm
distance, path = nx.single_source_dijkstra(G, source=r_smiles, target=p_smiles)
print(f"Shortest Path from reactant to Criegee: {path} with distance {distance}")

# con.draw_network(G)
#plt.savefig('ozonolysis.png')
print("It took ", end_time - start_time," seconds to create the ozonolysis reaction network.")
