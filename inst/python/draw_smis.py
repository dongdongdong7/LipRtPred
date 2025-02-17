from rdkit.Chem import Draw
from rdkit import Chem

def flatten_tuple(nested):
  flat_list = []
  for item in nested:
    if isinstance(item, tuple):
      flat_list.extend(flatten_tuple(item))
    else:
      flat_list.append(item)
  return list(flat_list)

def draw_smis(smis, filePath, SMARTS = None, molsPerRow = 1, subImgSize = (400, 400)):
  mols = []
  highlights = []
  for smi in smis:
    mol = Chem.MolFromSmiles(smi)
    if SMARTS is not None:
      substructure = Chem.MolFromSmarts(SMARTS)
      highlights.append(flatten_tuple(mol.GetSubstructMatches(substructure)))
    mols.append(mol)
  img = Draw.MolsToGridImage(mols, highlightAtomLists = highlights, molsPerRow = molsPerRow, subImgSize = subImgSize)
  img.save(filePath)
