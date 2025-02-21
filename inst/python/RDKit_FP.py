from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import MACCSkeys
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

# smi is a string
# smis is a list

# Morgan(ECFP)
def getMorganFP_loop_py(smi, radius = 2, fpSize = 2048, count = False):
  mol = Chem.MolFromSmiles(smi)
  mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius = radius, fpSize = fpSize)
  if count:
    fp = mfpgen.GetCountFingerprints(mol)
  else:
    fp = mfpgen.GetFingerprint(mol)
  return(fp.ToList())

def getMorganFP_py(smis, radius = 2, fpSize = 2048, count = False, thread = 1):
  res = Parallel(n_jobs=thread)(delayed(getMorganFP_loop_py)(smi, radius = radius, fpSize = fpSize, count = count) for smi in tqdm(smis))
  return(pd.DataFrame(res))

# Feature Morgan
def getFMorganFP_loop_py(smi, radius = 2, fpSize = 2048, count = False):
  mol = Chem.MolFromSmiles(smi)
  fmgen = rdFingerprintGenerator.GetMorganGenerator(radius = radius, fpSize = fpSize, atomInvariantsGenerator = rdFingerprintGenerator.GetMorganFeatureAtomInvGen())
  if count:
    fp = fmgen.GetCountFingerprint(mol)
  else:
    fp = fmgen.GetFingerprint(mol)
  return(fp.ToList())

def getFMorganFP_py(smis, radius = 2, fpSize = 2048, count = False, thread = 1):
  res = Parallel(n_jobs=thread)(delayed(getFMorganFP_loop_py)(smi, radius = radius, fpSize = fpSize, count = count) for smi in tqdm(smis))
  return(pd.DataFrame(res))

# RDKit
def getRDKitFP_loop_py(smi, minPath = 1, maxPath = 7, useHs = True, fpSize = 2048, count = False):
  mol = Chem.MolFromSmiles(smi)
  rdkgen = rdFingerprintGenerator.GetRDKitFPGenerator(minPath = minPath, maxPath = maxPath, useHs = useHs, fpSize = fpSize)
  if count:
    fp = rdkgen.GetCountFingerprint(mol)
  else:
    fp = rdkgen.GetFingerprint(mol)
  return(fp.ToList())

def getRDKitFP_py(smis, minPath = 1, maxPath = 7, useHs = True, fpSize = 2048, count = False, thread = 1):
  res = Parallel(n_jobs=thread)(delayed(getRDKitFP_loop_py)(smi, minPath=minPath, maxPath=maxPath, useHs=useHs, fpSize=fpSize, count=count) for smi in tqdm(smis))
  return(pd.DataFrame(res))

# Atom Pairs
def getApFP_loop_py(smi, minDistance = 1, maxDistance = 30, fpSize = 2048, count = False):
  mol = Chem.MolFromSmiles(smi)
  apgen = rdFingerprintGenerator.GetAtomPairGenerator(minDistance = minDistance, maxDistance = maxDistance, fpSize = fpSize)
  if count:
    fp = apgen.GetCountFingerprint(mol)
  else:
    fp = apgen.GetFingerprint(mol)
  return(fp.ToList())

def getApFP_py(smis, minDistance = 1, maxDistance = 30, fpSize = 2048, count = False, thread = 1):
  res = Parallel(n_jobs=thread)(delayed(getApFP_loop_py)(smi, minDistance=minDistance, maxDistance=maxDistance,fpSize=fpSize,count=count) for smi in tqdm(smis))
  return(pd.DataFrame(res))

# Topological Torsion
def getTtFP_loop_py(smi, fpSize = 2048, count = False):
  mol = Chem.MolFromSmiles(smi)
  ttgen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize = fpSize)
  if count:
    fp = ttgen.GetCountFingerprint(mol)
  else:
    fp = ttgen.GetFingerprint(mol)
  return(fp.ToList())

def getTtFP_py(smis, fpSize = 2048, count = False, thread = 1):
  res = Parallel(n_jobs=thread)(delayed(getTtFP_loop_py)(smi,fpSize=fpSize,count=count) for smi in tqdm(smis))
  return(pd.DataFrame(res))

# MACCSkeys
def getMaccsFP_loop_py(smi):
  mol = Chem.MolFromSmiles(smi)
  fp = MACCSkeys.GenMACCSKeys(mol)
  return(fp.ToList())

def getMaccsFP_py(smis, thread = 1):
  res = Parallel(n_jobs=thread)(delayed(getMaccsFP_loop_py)(smi) for smi in tqdm(smis))
  return(pd.DataFrame(res))
