import pandas as pd
import numpy as np
import conversion as con
from rdkit.Chem.rdmolfiles import MolFromSmiles

node_communities = pd.read_csv("ozone_modularity_nodes.csv")
print(node_communities.head())

edges = pd.read_csv("ozonolysis_edges.csv")
print(edges.head())




# -------------- community 0 -------------------------
# it has an epoxide ring

r0 = ['[H][CH3].[H][CH]([H])[H].[O]1[O]=[O]1', '[H][CH3].[H][CH]([H])[H].[O]1[O][O]1', '[H][CH]([H])[H].[H][CH]=[O][O][OH]','[H][CH]([O]1[O][O]1)[C]([H])([H])[H]','[H][C]([H])=[C]([H])[H].[O]1[O][O]1'
      , '[H][CH2][C]([H])([H])[H].[O]1[O][O]1', '[H][CH2][H].[H][CH2][H].[O]1[O]=[O]1', '[HH].[H][CH2][C]([H])([H])[O]1[O][O]1', '[H][CH2][O]([OH])[OH].[H][CH]([H])[H]',
      '[HH].[H][CH]=[C]([H])[H].[O]1[O][O]1', '[H][CH2][H].[H][CH]([H])[O]1[O][O]1', '[H][CH]([H])[C]([H])([H])[O]1[O][O]1',
      '[H][C]([H])([H])[C]1([H])[O][O]1[OH]']





# ------------- community 1 ---------------------------
# they all have ozone core / ozone fragment
# they have C2H2
# most ozone and C2H2 structures are connected
r1 = ['[H][C]([H])=[CH2].[H][O][O]=[O][H]', '[CH3][CH3].[HH].[HH].[HH].[HH].[OH][O][OH]', '[H][CH]([H])[CH3].[H][O][O]([H])[OH]', '[H][C]([H])=[CH2].[H][O]([H])[O]=[O]', '[H][C]([H])([CH3])[O]([OH])[OH].[H][H]',
      '[H][CH]([H])[CH3].[H][H].[OH2].[OH][OH]', '[HH].[H][C]#[CH].[H][H].[O]=[O][OH]', '[H][C]([H])([CH3])[O][O][OH].[H][H]', '[H][C]([H])([H])[CH3].[H][O][O][OH]', '[HH].[HH].[H][C]#[CH].[H][O]([OH])[OH]', '[H][C]([H])=[C]=[O][O][OH].[H][H]',
      '[H][CH]([H])[CH2][O]([OH])[OH].[H][H]', '[HH].[H][CH]([H])[CH3].[H][O][OH].[OH2]', '[HH].[H][O][O]([OH])[C]([H])([H])[CH3]']



# ---------------- community 2 ------------------------
# they all have ozone core / ozone fragment
# they have C2H2
# almost all ozone and C2H2 structures are disconnected
r2 = ['[HH].[HH].[H][C]#[C][H].[OH2].[O]=[O]', '[HH].[HH].[H][C]([H])([CH3])[O]=[O][OH]', '[HH].[H][CH2][C]([H])([H])[O][OH].[OH2]', '[HH].[H][CH]([H])[CH3].[H][O]([OH])[OH]', '[HH].[H][C]([H])=[CH2].[H][O]=[O][OH]', '[H][CH2][CH]([H])[H].[H][O][OH].[OH2]',
      '[H][CH]([H])[CH3].[H][H].[OH][O][OH]', '[H][CH]([H])[O]([O]=[O])[CH]([H])[H]', '[H][C]#[C][H].[H][H].[O]=[O][OH]', '[H][C]1([H])[O]([OH])[C]1([H])[H].[OH2]', '[H][C]1([H])[O]=[O][C]1([H])[H].[OH2]', '[HH].[HH].[HH].[H][C]#[CH].[OH][O][OH]', '[HH].[H][C]#[C][H].[H][O][O][OH]',
      '[HH].[H][C]([H])=[CH2].[H][O][O]=[O]', '[HH].[H][CH2][C]([H])([H])[OH].[OH][OH]', '[HH].[HH].[H][C]([H])([CH3])[O][O]=[O]', '[H][CH2][H].[H][CH]([H])[O][OH].[OH2]', '[HH].[H][C]#[C][H].[H][O]([OH])[OH]']

    
    
# -------------- community 3 --------------------------
# they all have ozone core / ozone fragment
# these ozone core fragments are connected on both ends with something
# have two components 
# that is why it forms 4 membered rings sometimes because each O is connected to the same thing on the ends
r3 = ['[H][CH3].[H][CH]([H])[O]([H])[O][OH]', '[H][CH]([H])[O]1[O]([OH])[C]1([H])[H]', '[HH].[H][CH3].[H][C]1([H])[O][O]=[O]1', '[H][CH]([H])[O]1[O][O][C]1([H])[H]', '[H][CH]([H])[O][O]1[O][C]1([H])[H]', '[H][CH]([H])[O]1[O][C]([H])([H])[O]1', '[HH].[H][CH2][O]([O][OH])[CH]([H])[H]',
      '[HH].[H][CH3].[H][C]1([H])[O][O]1[OH]', '[H][CH2][H].[H][C]([H])([OH])[O][OH]', '[HH].[H][CH2][O]([OH])[O][CH]([H])[H]', '[H][CH2][H].[H][C]1([H])[O][O]1.[OH2]', '[H][CH]([H])[O]([OH])[CH]([H])[H].[OH2]', '[H][C]([H])=[O][CH]([H])[H].[O]=[O]', '[HH].[H][CH]=[O][O][O][CH]([H])[H]',
      '[H][C]([H])=[O][O]=[C]([H])[H].[OH2]', '[HH].[H][CH]=[O][O]([OH])[CH]([H])[H]', '[H][C]([H])=[O].[H][C]([H])=[O][OH]', '[H][CH3].[H][CH]([H])[O]([H])[O]=[O]', '[H][CH2][H].[H][CH2][O]=[O][O][H]', '[H][CH3].[H][O]=[O][O]=[C]([H])[H]', '[H][CH2][H].[H][CH]([H])[O][O][OH]',
      '[H][CH3].[H][CH]([H])[O][O]([H])[OH]']


# --------------- community 4 ---------------------------
# could not really find a pattern...
r4 = ['[H][CH]=[C]([H])[O][O]=[O].[H][H]', '[HH].[HH].[H][C]([H])=[C]=[O].[OH][OH]', '[HH].[H][C]#[C][H].[H][OH].[O]=[O]', '[HH].[HH].[H][CH]1[CH]([H])[O]1[O][OH]', '[HH].[H][C]#[C][H].[H][O]=[O].[OH2]', '[HH].[H][CH]=[C]([H])[O]([H])[O][OH]', '[HH].[H][C]1([H])[O]2[O]([OH])[C]21[H]',
      '[H][C]1([H])[O]([O][OH])[C]1([H])[H]', '[HH].[HH].[H][CH]=[C]([H])[OH].[O]=[O]', '[H][CH]=[C]([H])[H].[H][O]=[O][OH]', '[HH].[H][C]1([H])[O]2[O][O][C]21[H]', '[H][CH]([H])[C]1([H])[O][O][O]1[H]', '[H][C]([H])([OH])[C]([H])([H])[OH].[OH2]', '[HH].[H][CH]=[C]([H])[O][O][O][H]', '[H][CH]([H])[C]1([H])[O]([H])[O]1[OH]',
      '[HH].[HH].[H][C]#[C][H].[OH][O][OH]', '[H][O][O][O]1[CH]([H])[C]1([H])[H]', '[H][O]=[O][C]([H])=[C]([H])[H].[OH2]']

for i in r4:
    struct = MolFromSmiles(i, sanitize = False)
    con.draw_molecules(struct)