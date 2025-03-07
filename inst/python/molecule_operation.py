# Molecular operations with RDKit
# Barry Song
# 250307

from rdkit import Chem
from rdkit.Chem import rdmolops

# smi SMILES string
# smis: SMILES string list
# SMARTS: SMARTS pattern
# start_atom_idx: start atom index
# end_atom_idx: end atom index
# atom_idx_list: atom index list
# non_traversable_atom_idx: list of atom index that are not allowed to be traversed

# Get SMILES with atom index
def Get_SMILES_with_index_py(smi):
  mol = Chem.MolFromSmiles(smi)
  for atom in mol.GetAtoms():
    atom.SetAtomMapNum(atom.GetIdx())
  return(Chem.MolToSmiles(mol))

# Searching for substructures of molecules using SMARTS
def GetSubstructMatches_py(smis, SMARTS):
  mols = []
  for smi in smis:
    mol = Chem.MolFromSmiles(smi)
    mols.append(mol)
  match = []
  for mol in mols:
    match.append(mol.GetSubstructMatches(Chem.MolFromSmarts(SMARTS)))
  return(match)

# Calculate the shortest distance between two atoms (least number of atoms)
def GetShortestLength_py(smi, start_atom_idx, end_atom_idx):
  mol = Chem.MolFromSmiles(smi)
  start_atom = mol.GetAtomWithIdx(start_atom_idx)
  end_atom = mol.GetAtomWithIdx(end_atom_idx)
  paths = rdmolops.GetShortestPath(mol, start_atom_idx, end_atom_idx)
  return(len(paths) - 2)

# Get symbol of atoms
def GetAtomSymbol_py(smi, atom_idx_list):
  mol = Chem.MolFromSmiles(smi)
  symbolList = []
  for atom_idx in atom_idx_list:
    atom = mol.GetAtomWithIdx(atom_idx)
    symbolList.append(atom.GetSymbol())
  return symbolList

# Break all the rings in the molecule
def BreakRings_py(smi):
  mol = Chem.MolFromSmiles(smi)
  if mol is None:
    raise ValueError("Wrong SMILES")
  
  rw_mol = Chem.RWMol(mol)
  
  ring_info = rw_mol.GetRingInfo()
  rings = ring_info.AtomRings()
  rings = tuple(sorted(rings, key = len))
  
  for ring in rings:
    single_bonds_CC = []
    single_bonds_Other = []
    double_bonds_CC = []
    double_bonds_Other = []
    triple_bonds_CC = []
    triple_bonds_Other = []
    
    for i in range(len(ring)):
      for j in range(i + 1, len(ring)):
        bond = rw_mol.GetBondBetweenAtoms(ring[i],ring[j])
        if bond is not None:
          bond_type = bond.GetBondType()
          atom1 = bond.GetBeginAtom().GetSymbol()
          atom2 = bond.GetEndAtom().GetSymbol()
          
          if bond_type == Chem.BondType.SINGLE:
            if atom1 == "C" and atom2 == "C":
              single_bonds_CC.append(bond)
            else:
              single_bonds_Other.append(bond)
          elif bond_type == Chem.BondType.DOUBLE:
            if atom1 == "C" and atom2 == "C":
              double_bonds_CC.append(bond)
            else:
              double_bonds_Other.append(bond)
          elif bond_type == Chem.BondType.TRIPLE:
            if atom1 == "C" and atom2 == "C":
              triple_bonds_CC.append(bond)
            else:
              triple_bonds_Other.append(bond)
      
    if single_bonds_CC:
      bond_to_break = single_bonds_CC[0]
    elif single_bonds_Other:
      bond_to_break = single_bonds_Other[0]
    elif double_bonds_CC:
      bond_to_break = double_bonds_CC[0]
    elif double_bonds_Other:
      bond_to_break = double_bonds_Other[0]
    elif triple_bonds_CC:
      bond_to_break = triple_bonds_CC[0]
    elif triple_bonds_Other:
      bond_to_break = triple_bonds_Other[0]
    else:
      continue
    
    rw_mol.RemoveBond(bond_to_break.GetBeginAtom().GetIdx(), bond_to_break.GetEndAtom().GetIdx())
    
    mol = rw_mol.GetMol()
    for atom in mol.GetAtoms():
      atom.SetAtomMapNum(atom.GetIdx())
    
    return(Chem.MolToSmiles(mol))

# Traverse the entire molecule starting from the start atom, excluding the specified atoms
def TraverseMolecule_py(smi, start_atom_idx, non_traversable_atom_idx):
  mol = Chem.MolFromSmiles(smi)
  if mol is None:
    raise ValueError("Wrong SMILES")
  
  if start_atom_idx < 0 or start_atom_idx >= mol.GetNumAtoms():
    raise ValueError("Wrong start_atom_idx")
  
  visited_initial = set(non_traversable_atom_ids)
  path = []
  
  def dfs(atom_idx, visited):
    if atom_idx in visited:
      return
    visited.add(atom_idx)
    path.append(atom_idx)
    atom = mol.GetAtomWithIdx(atom_idx)
    for neighbor in atom.GetNeighbors():
      dfs(atom_idx = neighbor.GetIdx(), visited = visited)
  
  dfs(start_atom_idx, visited = visited_initial)
  return(path)
