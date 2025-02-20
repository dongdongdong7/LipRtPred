from rdkit.Chem import Descriptors
from rdkit import Chem
from joblib import Parallel, delayed
from tqdm import tqdm
import pandas as pd

def getRDKitMD_loop_py(smi):
  mol = Chem.MolFromSmiles(smi)
  desc = Descriptors.CalcMolDescriptors(mol)
  return desc

def getRDKitMD_py(smis, thread = 1):
  res = Parallel(n_jobs=thread)(delayed(getRDKitMD_loop_py)(smi) for smi in tqdm(smis))
  return pd.DataFrame(res)
