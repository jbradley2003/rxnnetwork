import conversion as con
import numpy as np
import ac2mol
import networkx as nx
import matplotlib.pyplot as plt
from rdkit.Chem.rdmolfiles import MolFromSmiles




p1 = MolFromSmiles("[C](=[O])(-[H])-[H].[O](=[C](-[H])-[H])-[OH]", sanitize = False)

# Criegee 1
r = MolFromSmiles("[C](=[C](-[H])-[H])(-[H])-[H].[O]=[O]-[OH]", sanitize = False)
i1 = MolFromSmiles("[CH3]-[H].[CH2](-[H])-[H].[HH].[OH]-[O]-[OH]", sanitize = False)
i2 = MolFromSmiles("[CH3]-[O]-[O]-[OH].[HH].[CH2](-[H])-[H].[HH]", sanitize = False)
i3 = MolFromSmiles("[CH](=[O]-[O]-[O]-[CH](-[H])-[H])-[H].[HH]", sanitize = False)
reactant = MolFromSmiles("[C](=[O])(-[H])-[H].[OH]-[O]=[C](-[H])-[H]", sanitize = False)


# Criegee 2
# [CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]				# 0.7070063948631287
# [CH]1(-[H])-[C](-[H])(-[H])-[O]-[O](-[H])-[O]-1			# 0.7946428656578064
# [CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1-[H]				# 0.8230769038200378
# [CH]1(-[H])-[C](-[H])(-[H])-[O](-[O]-[H])-[O]-1		# 0.8192307949066162
# [CH]1(-[H])-[C](-[H])(-[H])-[O]=[O]-1.[OH2].[HH]			 #0.3564356565475464
# [CH]1(-[H])-[C](-[H])(-[H])-[O](-[OH])-[O]-1-[H]			 #0.8507936596870422
# [CH]1(-[H])-[C](-[H])(-[H])-[O]-1-[OH].[OH2].[HH]				 #0.7428571581840515


# [CH]1(-[H])-[C](-[H])(-[H])-[O](-[OH])-[O]-1.[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	28262			0.450549453496933
# [CH]1(-[H])-[C](-[H])(-[H])-[O]-1-[O]-[OH].[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	28142			0.530386745929718
# [C](-[CH](-[H])-[H])(-[H])(-[O]-[O]=[O])-[H]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	27945			0.6089743375778198
# [C]1(-[CH](-[H])-[H])(-[H])-[O]-[O]-[O]-1.[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	27892			0.5493826866149902
# [C]1(-[CH](-[H])-[H])(-[H])-[O]-[O]-1-[OH].[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	27832			0.6288659572601318
# [CH](=[C](-[H])-[H])-[H].[OH]-[O]=[O].[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	27758			0.9496402740478516
# [CH]1(-[H])-[C](-[H])(-[H])-[O](-[H])-[O]-[O]-1	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	27140			0.46382978558540344
# [CH]1(-[H])-[C](-[H])(-[H])-[O]-[O](-[H])-[O]-1	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	27074			0.3471502661705017
# [CH](-[C](-[H])(-[H])-[OH])(-[H])-[O]=[O].[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	27005			0.7032257914543152
# [CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1-[H]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	26914			0.4497816562652588
# [CH]1(-[H])-[C](-[H])(-[H])-[O]2-[O]-1-[O]-2.[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	26613			0.4298642575740814
# [CH]1(-[H])-[CH](-[H])-[O]-[O]-[O]-1.[HH].[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	25511			0.4126984179019928
# [C]12(-[H])-[C](-[H])(-[H])-[O]-1-[O]-[O]-2.[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	24725			0.5699658989906311
# [C]12(-[H])-[C](-[H])(-[H])-[O]-[O]-1-[O]-2.[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	24208			0.5563380122184753
# [C](-[CH](-[H])-[H])(=[O]-[O]-[OH])-[H].[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	23416			0.8895348906517029
# [CH3]-[C](-[H])(-[H])-[O]-[O]=[O].[HH].[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	22353			0.73758864402771
# [CH](=[C](-[H])-[H])-[O]-[O]-[OH].[HH].[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	22014			0.9060402512550354
# [CH2]1-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH].[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	18598			0.2857142984867096
# [CH3]-[H].[C](-[H])(-[H])=[O]-[O]=[O].[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	7628			0.9452054500579834
# [CH2](-[H])-[O]-[O]=[O].[CH2](-[H])-[H].[HH]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	7489			0.8872180581092834
# [CH](-[H])=[C](-[H])-[H].[HH].[OH]-[O]=[O]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	7192			0.9496402740478516
# [CH](=[C](-[H])-[H])-[H].[HH].[OH]-[O]=[O]	[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]	Undirected	1797			0.9496402740478516



# [CH](=[C](-[H])-[H])-[H].[H]-[O]-[O]=[O]	[CH](=[C](-[H])-[H])-[H].[HH].[OH]-[O]=[O]	Undirected	1180			0.23076923191547394
# [C](=[C](-[H])-[H])(-[H])-[H].[O]=[O]-[OH]	[CH](=[C](-[H])-[H])-[H].[HH].[OH]-[O]=[O]	Undirected	25			0.09090909361839294











# ozonolysis 4 iterations
r = MolFromSmiles("[H][C]([H])=[C]([H])[H].[O]=[O][OH]", sanitize = False)
i1 = MolFromSmiles("[H][CH]([H])[C]([H])([H])[OH].[O]#[O]", sanitize = False)
i2 = MolFromSmiles("[H][C]([H])=[O].[H][C]([H])=[O][OH]", sanitize = False)



r = MolFromSmiles("[C](=[C](-[H])-[H])(-[H])-[H].[O]=[O]-[OH]", sanitize = False)
i1 = MolFromSmiles("[CH](=[C](-[H])-[H])-[H].[HH].[OH]-[O]=[O]", sanitize = False)
i2 = MolFromSmiles("[CH]1(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1.[HH]", sanitize = False)
i3 = MolFromSmiles("[CH](-[C](-[H])(-[H])-[O]=[O])(-[OH])-[H].[HH]", sanitize = False)
p = MolFromSmiles("[C](=[O])(-[H])-[H].[C](-[H])(-[H])=[O]-[OH]", sanitize = False)

# ozonolysis 3 iterations
r = MolFromSmiles("[H][C]([H])=[C]([H])[H].[O]=[O][OH]", sanitize = False)
i1 = MolFromSmiles("[HH].[H][CH]=[C]([H])[H].[O]=[O][OH]", sanitize = False)
i2 = MolFromSmiles("[HH].[H][CH3].[H][C]([H])=[O][O]=[O]", sanitize = False)
i3 = MolFromSmiles("[HH].[H][CH]=[O][OH].[H][C]([H])=[O]", sanitize = False)
p = MolFromSmiles("[H][C]([H])=[O].[H][C]([H])=[O][OH]", sanitize = False)

# paper example dist 7 iterations
r = MolFromSmiles("[Cl][I].[Cl][I].[H][H]", sanitize = False)
i1 = MolFromSmiles("[ClH].[ClH].[H][H].[IH].[IH]", sanitize = False)
i2 = MolFromSmiles("[ClH].[HH].[H][Cl].[IH].[IH]", sanitize = False)
i3 = MolFromSmiles("[H][Cl].[H][Cl].[I][I]", sanitize = False)

# paper example tanimoto 

i2 = MolFromSmiles("[CH]1(-[H])(-[H])-[C](-[H])(-[H])-[O]-[O]-[O]-1", sanitize = False)

con.draw_molecules(i2)



