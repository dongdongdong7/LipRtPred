from mordred import Calculator, descriptors
from rdkit import Chem
import pandas as pd

def getMordredMD_py(smis, ignore_3D = True, thread = 1):
  mols = [Chem.MolFromSmiles(smi) for smi in smis]
  calc = Calculator(descriptors, ignore_3D = ignore_3D)
  res_pd = calc.pandas(mols, nproc = thread)
  return(res_pd.astype(float))
