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
