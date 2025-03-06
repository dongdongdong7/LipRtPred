from rdkit import Chem
from rdkit.Chem import rdmolops

def GetSubstructMatches_py(smis, SMARTS):
  mols = []
  for smi in smis:
    mol = Chem.MolFromSmiles(smi)
    mols.append(mol)
  match = []
  for mol in mols:
    match.append(mol.GetSubstructMatches(Chem.MolFromSmarts(SMARTS)))
  return(match)

def C_C_count_py(smi, start_atom_idx, end_atom_idx):
  mol = Chem.MolFromSmiles(smi)
  C_start_atom = mol.GetAtomWithIdx(start_atom_idx)
  C_end_atom = mol.GetAtomWithIdx(end_atom_idx)
  if C_start_atom.GetSymbol() != "C":
    raise ValueError("Start atom is not C")
  if C_end_atom.GetSymbol() != "C":
    raise ValueError("End atom is not C")
  paths = rdmolops.GetShortestPath(mol, start_atom_idx, end_atom_idx) # a tuple with atom index
  return(len(paths) - 2)

def walk_away_py(smi, start_atom_idx, non_traversable_atom_ids):
  visited_initial = set(non_traversable_atom_ids)
  def dfs(atom_idx, visited):
    visited.add(atom_idx)
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetSymbol() == "C":
      count = 1  # 当前原子是C原子，所以计数1
    else:
      count = 0
    for neighbor in atom.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx not in visited: # and neighbor.GetSymbol() == 'C'
            count += dfs(neighbor_idx, visited)  # 递归遍历邻居
    return count

  mol = Chem.MolFromSmiles(smi)
  return(dfs(atom_idx = start_atom_idx, visited = visited_initial) - 1)


def traverse_molecule_py(smi, start_atom_idx, non_traversable_atom_ids):
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
        
  
def getAtomSymbol_py(smi, atom_idx_list):
  mol = Chem.MolFromSmiles(smi)
  symbolList = []
  for atom_idx in atom_idx_list:
    atom = mol.GetAtomWithIdx(atom_idx)
    symbolList.append(atom.GetSymbol())
  return symbolList

def break_all_rings_bond_py(smi):
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
