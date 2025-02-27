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
  return len(paths) - 2

def walk_away_py(smi, start_atom_idx, main_chains_atom_idx):
  visited_initial = set(main_chains_atom_idx)
  def dfs(atom_idx, visited):
    visited.add(atom_idx)
    count = 1  # 当前原子是C原子，所以计数1
    for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx not in visited and neighbor.GetSymbol() == 'C':
            count += dfs(neighbor_idx, visited)  # 递归遍历邻居
    return count
  
  mol = Chem.MolFromSmiles(smi)
  return(dfs(atom_idx = start_atom_idx, visited = visited_initial) - 1)
  
  
