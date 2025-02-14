from rdkit.Chem import Draw
from rdkit import Chem

def draw_smis(smis, filePath):
  mols = []
  for smi in smis:
    mol = Chem.MolFromSmiles(smi)
    mols.append(mol)
  img = Draw.MolsToGridImage(mols, returnPNG=False)
  img.save(filePath)
