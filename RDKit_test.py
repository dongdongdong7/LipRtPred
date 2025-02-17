from rdkit import Chem

smi = "C(O)(=O)CC(O)C/C=C\C/C=C\C/C=C\C/C=C\CCCCC"
mol = Chem.MolFromSmiles(smi)
mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3;$(C=O);$(C-O)](=O)O")) # 检测酯键 0
mol.GetSubstructMatches(Chem.MolFromSmarts("C=C"))
atom1 = mol.GetAtomWithIdx(7)
atom1_neighbors = atom1.GetNeighbors()
atom1_neighbors[0].GetIdx() # 6 -> 0
atom1_neighbors[1].GetIdx()
atom2 = mol.GetAtomWithIdx(6)
atom2_neighbors = atom2.GetNeighbors()
atom2_neighbors[0].GetIdx() # 4 -> 0
atom2_neighbors[0].GetSymbol()
atom2_neighbors[1].GetIdx()
atom3 = mol.GetAtomWithIdx(4)
atom3_neighbors = atom3.GetNeighbors()
atom3_neighbors[0].GetIdx()
atom3_neighbors[0].GetSymbol()
atom3_neighbors[1].GetIdx()
atom3_neighbors[1].GetSymbol()
atom3_neighbors[2].GetIdx()
atom3_neighbors[2].GetSymbol()
atom_O = mol.GetAtomWithIdx(5)
atom_O.GetNeighbors()
atom2.GetSymbol()
