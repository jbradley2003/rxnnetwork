"""
Module for generating rdkit molobj/smiles/molecular graph from adjacency matrix

Implementation by Jeff Wang & James Bradley based on implementation by Jan H. Jensen and the paper

    Yeonjoon Kim and Woo Youn Kim
    "Universal Structure Conversion Method for Organic Molecules: From Atomic Connectivity
    to Three-Dimensional Geometry"
    Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777
    DOI: 10.1002/bkcs.10334

"""

import copy
import itertools
    
from collections import defaultdict

import numpy as np
import networkx as nx

from rdkit import Chem
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
from rdkit.Chem.rdmolfiles import MolToSmiles
from rdkit.Chem.rdmolfiles import MolFromSmiles
import sys

global __ATOM_LIST__
__ATOM_LIST__ = \
    ['h',  'he',
     'li', 'be', 'b',  'c',  'n',  'o',  'f',  'ne',
     'na', 'mg', 'al', 'si', 'p',  's',  'cl', 'ar',
     'k',  'ca', 'sc', 'ti', 'v ', 'cr', 'mn', 'fe', 'co', 'ni', 'cu',
     'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',
     'rb', 'sr', 'y',  'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag',
     'cd', 'in', 'sn', 'sb', 'te', 'i',  'xe',
     'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy',
     'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w',  're', 'os', 'ir', 'pt',
     'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn',
     'fr', 'ra', 'ac', 'th', 'pa', 'u',  'np', 'pu']


global atomic_valence
global atomic_valence_electrons

atomic_valence = defaultdict(list)
atomic_valence[1] = [1]
atomic_valence[5] = [3,4]
atomic_valence[6] = [4]
atomic_valence[7] = [3,4]
atomic_valence[8] = [2,1,3]
atomic_valence[9] = [1]
atomic_valence[14] = [4]
atomic_valence[15] = [5,3] #[5,4,3]
atomic_valence[16] = [6,3,2] #[6,4,2]
atomic_valence[17] = [1]
atomic_valence[32] = [4]
atomic_valence[35] = [1]
atomic_valence[53] = [1]

atomic_valence_electrons = {}
atomic_valence_electrons[1] = 1
atomic_valence_electrons[5] = 3
atomic_valence_electrons[6] = 4
atomic_valence_electrons[7] = 5
atomic_valence_electrons[8] = 6
atomic_valence_electrons[9] = 7
atomic_valence_electrons[14] = 4
atomic_valence_electrons[15] = 5
atomic_valence_electrons[16] = 6
atomic_valence_electrons[17] = 7
atomic_valence_electrons[32] = 4
atomic_valence_electrons[35] = 7
atomic_valence_electrons[53] = 7


def str_atom(atom):
    """
    convert integer atom to string atom
    """
    global __ATOM_LIST__
    atom = __ATOM_LIST__[atom - 1]
    return atom


def int_atom(atoms):
    """ takes in a list of atoms in string form and converting them all to integer form

    Args:
        atoms (list): list of atoms in string form

    Returns:
        list: same list of atoms but in integer form
    """
    global __ATOM_LIST__

    
    atoms_int = []
    for i in atoms:
        i = i.lower()
        atoms_int.append(__ATOM_LIST__.index(i) + 1)
        
    return atoms_int


def get_UA(maxValence_list, valence_list):
    """
    """
    UA = []
    DU = []
    for i, (maxValence, valence) in enumerate(zip(maxValence_list, valence_list)):
        if not maxValence - valence > 0:
            continue
        UA.append(i)
        DU.append(maxValence - valence)
    return UA, DU


def get_BO(AC, UA, DU, valences, UA_pairs, use_graph=True):
    """
    """
    BO = AC.copy()
    DU_save = []

    while DU_save != DU:
        for i, j in UA_pairs:
            BO[i, j] += 1
            BO[j, i] += 1

        BO_valence = list(BO.sum(axis=1))
        DU_save = copy.copy(DU)
        UA, DU = get_UA(valences, BO_valence)
        UA_pairs = get_UA_pairs(UA, AC, use_graph=use_graph)[0]

    return BO


def valences_not_too_large(BO, valences):
    """
    """
    number_of_bonds_list = BO.sum(axis=1)
    for valence, number_of_bonds in zip(valences, number_of_bonds_list):
        if number_of_bonds > valence:
            return False

    return True

def charge_is_OK(BO, AC, charge, DU, atomic_valence_electrons, atoms, valences,
                 allow_charged_fragments=True):
    # total charge
    Q = 0

    # charge fragment list
    q_list = []

    if allow_charged_fragments:

        BO_valences = list(BO.sum(axis=1))
        for i, atom in enumerate(atoms):
            q = get_atomic_charge(atom, atomic_valence_electrons[atom], BO_valences[i])
            Q += q
            if atom == 6:
                number_of_single_bonds_to_C = list(BO[i, :]).count(1)
                if number_of_single_bonds_to_C == 2 and BO_valences[i] == 2:
                    Q += 1
                    q = 2
                if number_of_single_bonds_to_C == 3 and Q + 1 < charge:
                    Q += 2
                    q = 1

            if q != 0:
                q_list.append(q)

    return (charge == Q)

def BO_is_OK(BO, AC, charge, DU, atomic_valence_electrons, atoms, valences,
    allow_charged_fragments=True):
    """
    Sanity of bond-orders

    args:
        BO -
        AC -
        charge -
        DU - 


    optional
        allow_charges_fragments - 


    returns:
        boolean - true of molecule is OK, false if not
    """

    if not valences_not_too_large(BO, valences):
        return False

    check_sum = (BO - AC).sum() == sum(DU)
    check_charge = charge_is_OK(BO, AC, charge, DU, atomic_valence_electrons, atoms, valences,
                                allow_charged_fragments)

    if check_charge and check_sum: 
        return True

    return False


def get_atomic_charge(atom, atomic_valence_electrons, BO_valence):
    """
    """

    if atom == 1:
        charge = 1 - BO_valence
    elif atom == 5:
        charge = 3 - BO_valence
    elif atom == 15 and BO_valence == 5:
        charge = 0
    elif atom == 16 and BO_valence == 6:
        charge = 0
    else:
        charge = atomic_valence_electrons - 8 + BO_valence

    return charge



def set_atomic_radicals(mol, atoms, atomic_valence_electrons, BO_valences,
                                                use_atom_maps):
    """

    The number of radical electrons = absolute atomic charge

    """
    for i, atom in enumerate(atoms):
        a = mol.GetAtomWithIdx(i)
        if use_atom_maps:
            a.SetAtomMapNum(i+1)
        charge = get_atomic_charge(
            atom,
            atomic_valence_electrons[atom],
            BO_valences[i])

        if (abs(charge) > 0):
            a.SetNumRadicalElectrons(abs(int(charge)))

    return mol


def get_bonds(UA, AC):
    """

    """
    bonds = []

    for k, i in enumerate(UA):
        for j in UA[k + 1:]:
            if AC[i, j] == 1:
                bonds.append(tuple(sorted([i, j])))

    return bonds


def get_UA_pairs(UA, AC, use_graph=True):
    """

    """

    bonds = get_bonds(UA, AC)

    if len(bonds) == 0:
        return [()]

    if use_graph:
        G = nx.Graph()
        G.add_edges_from(bonds)
        UA_pairs = [list(nx.max_weight_matching(G))]
        return UA_pairs

    max_atoms_in_combo = 0
    UA_pairs = [()]
    for combo in list(itertools.combinations(bonds, int(len(UA) / 2))):
        flat_list = [item for sublist in combo for item in sublist]
        atoms_in_combo = len(set(flat_list))
        if atoms_in_combo > max_atoms_in_combo:
            max_atoms_in_combo = atoms_in_combo
            UA_pairs = [combo]

        elif atoms_in_combo == max_atoms_in_combo:
            UA_pairs.append(combo)

    return UA_pairs


def AC2BO(AC, atoms, charge, allow_charged_fragments=True, use_graph=True):
    """
    * Written by Jan H. Jensen *
    
    Todo: 
        still need to figure out a check for invalid BOs
        right now if it is invalid it just crashes the code
    
    implemenation of algorithm shown in Figure 2

    UA: unsaturated atoms

    DU: degree of unsaturation (u matrix in Figure)

    best_BO: Bcurr in Figure

    """

    global atomic_valence
    global atomic_valence_electrons

    # make a list of valences, e.g. for CO: [[4],[2,1]]
    valences_list_of_lists = []
    AC_valence = list(AC.sum(axis=1))
    
    # combines them into pairs and keeps track of their index
    for i,(atomicNum,valence) in enumerate(zip(atoms,AC_valence)):
        # valence can't be smaller than number of neighbourgs
        possible_valence = [x for x in atomic_valence[atomicNum] if x >= valence]
        if not possible_valence:
            #print('Valence of atom',i,'is',valence,'which bigger than allowed max',max(atomic_valence[atomicNum]),'. Stopping')
            # fill in something that will allow for bad structures, just not do anything
            # sys.exit()
            return 
        valences_list_of_lists.append(possible_valence)

    # convert [[4],[2,1]] to [[4,2],[4,1]]
    valences_list = itertools.product(*valences_list_of_lists)

    best_BO = AC.copy()

    for valences in valences_list:

        UA, DU_from_AC = get_UA(valences, AC_valence)

        check_len = (len(UA) == 0)
        if check_len:
            check_bo = BO_is_OK(AC, AC, charge, DU_from_AC,
                atomic_valence_electrons, atoms, valences,
                allow_charged_fragments=allow_charged_fragments)
        else:
            check_bo = None

        if check_len and check_bo:
            return AC, atomic_valence_electrons

        UA_pairs_list = get_UA_pairs(UA, AC, use_graph=use_graph)
        for UA_pairs in UA_pairs_list:
            BO = get_BO(AC, UA, DU_from_AC, valences, UA_pairs, use_graph=use_graph)
            status = BO_is_OK(BO, AC, charge, DU_from_AC,
                        atomic_valence_electrons, atoms, valences,
                        allow_charged_fragments=allow_charged_fragments)
            charge_OK = charge_is_OK(BO, AC, charge, DU_from_AC, atomic_valence_electrons, atoms, valences,
                                     allow_charged_fragments=allow_charged_fragments)

            if status:
                return BO, atomic_valence_electrons
            elif BO.sum() >= best_BO.sum() and valences_not_too_large(BO, valences) and charge_OK:
                best_BO = BO.copy()

    return best_BO, atomic_valence_electrons


def BO2Mol(BO, atoms):
    """ takes in a bond order matrix and converts it to a RDkit molecular graph

    Args:
        BO (numpy array): the bond order matrix that describes our structure
        atoms (list): the list of atomic numbers in our structure
    
    Returns:
        mol (RDKit Molecule Object): returns the best matched molecule for this bond order matrix
    """
    rwMol = Chem.RWMol()
    
    BO = BO[0]
    
    # you cannot have a bond order larger than 3

    if np.any(BO > 3):
        return None
    
    
    # add the atoms
    atom_indices = []
    for i in atoms:
        atom = Chem.Atom(i)
        idx = rwMol.AddAtom(atom)
        atom_indices.append(idx)
        
    
    # add the bonds
    for i in range(len(BO)):
        for j in range(i+1, len(BO)):  # Only iterate upper triangular matrix
            bond_order = BO[i][j]
            if bond_order > 0:
                # Determine bond type based on bond order
                if bond_order == 1:
                    bond_type = Chem.BondType.SINGLE
                elif bond_order == 2:
                    bond_type = Chem.BondType.DOUBLE
                elif bond_order == 3:
                    bond_type = Chem.BondType.TRIPLE
                else:
                    raise ValueError(f"Unsupported bond order: {bond_order}")
                
                # Add the bond
                rwMol.AddBond(i, j, bond_type)
    
    mol = rwMol.GetMol()
    
    # check if this structure can be part of unimolecular reaction
    # if len(Chem.GetMolFrags(mol)) > 1:
    #     return None
    return mol


def AC2Mol(AC, atoms, charge = 0 , allow_charged_fragments=True, use_graph=True):
    """ Create a molecular graph based strictly on an adjacency matrix

    Args:
        AC (numpy array): the adjacency matrix that we are converting to a molecular graph
        atoms (list): list of atomic numbers corresponding to elements in our molecule
        charge (int, optional):the charge of our structure, 0 by default
        allow_charged_fragments (bool, optional): don't really know what this is
        use_graph (bool, optional): uses graph or not
        
    Returns:
        best_mol (RDKit Molecular Object): The ideal structure that came from the bond order matrix
        isomers (list): other stereoisomers in a list form
    """
    
    # create the bond order matrix
    BO = AC2BO(AC, atoms, charge, allow_charged_fragments, use_graph)
    if BO == None:
        return None
    
    # the molecule that matches up with the best bond order matrix
    best_mol = BO2Mol(BO, atoms)
    if best_mol == None:
        return None
    
    return best_mol


def Mol2Smiles(mol):
    return MolToSmiles(mol, isomericSmiles = True, allHsExplicit=True, ignoreAtomMapNumbers=False, canonical=False, allBondsExplicit = True)

def Smiles2AC(smiles):
    mol =  MolFromSmiles(smiles, sanitize = False)
    ac, atoms = Mol2AC(mol)
    return ac, atoms

def AC2Smiles(AC, atoms):
    mol = AC2Mol(AC, atoms)
    smiles = Mol2Smiles(mol)
    return smiles
def Mol2AC(mol):
    """ converts RDKit molecule back into AC

    Args:
        mol (RDKit molecule): molecule in RDKit form

    Returns:
        AC: the AC matrix
        atom_int: the list of atoms in integer form
    """
    # adjacency matrix
    AC = GetAdjacencyMatrix(mol)
    # list of atom symbols
    atom_str = [atom.GetSymbol() for atom in mol.GetAtoms()]
    atom_int = int_atom(atom_str)
    
    return AC, atom_int


def chiral_stereo_check(mol):
    """
    Find and embed chiral information into the model based on the coordinates

    args:
        mol - rdkit molecule, with embeded conformer

    """
    Chem.SanitizeMol(mol)
    Chem.DetectBondStereochemistry(mol, -1)
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True, force=True)
    Chem.AssignAtomChiralTagsFromStructure(mol, -1)

    return