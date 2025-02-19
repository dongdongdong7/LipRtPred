from mordred import Calculator, descriptors
from rdkit import Chem
import pandas as pd

smi = "C(O)(=O)CC(O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC"
mol = Chem.MolFromSmiles(smi)
calc = Calculator(descriptors, ignore_3D = True)
X_mord = pd.DataFrame(calc.pandas([mol]))

