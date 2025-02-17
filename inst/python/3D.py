from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def get_atom_distance_py(smi, atom_idx1, atom_idx2):
    # 获取分子的3D坐标
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol) 
    conf = mol.GetConformer()
    
    # 获取两个原子的3D坐标（x, y, z）
    pos1 = conf.GetAtomPosition(atom_idx1)
    pos2 = conf.GetAtomPosition(atom_idx2)
    
    # 计算两点之间的欧几里得距离
    distance = np.linalg.norm(np.array([pos1.x, pos1.y, pos1.z]) - np.array([pos2.x, pos2.y, pos2.z]))
    
    return distance
