# Search substructure of Glycerophospholipids (GP)
# Barry Song
# 250311

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py

# 1. Phosphate group
# (1) Search phosphate group
# .searchPhosphate(smi = "C(O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](COP(=O)(O)OCCN)O1)[C@]([H])(OC(CCCCCCC/C=C\\CCCCCC)=O)CO/C=C\\CCCCCCCCCCCCCC",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhosphate <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)(O)O",
                       scriptPath = scriptPath)[[1]]
}

# 2. Choline
# (1) Search choline
# .searchCholine(smi = "[C@](COP(=O)([O-])OCC[N+](C)(C)C)([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
#                scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchCholine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[N+](C)(C)(C)CCO",
                       scriptPath = scriptPath)[[1]]
}

# 3. Ethanolamine
# (1) Search ethanolamine
# .searchEthanolamine(smi = "[C@](COP(=O)(O)OCCN)([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
#                     scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchEthanolamine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "OCC[N;NH2]",
                       scriptPath = scriptPath)[[1]]
}

# 4. Serine
# (1) Search serine
# .searchSerine(smi = "C(O)(=O)[C@@]([H])(N)COP(OC[C@]([H])(OC(CCCCCCCCCCCC)=O)COC(CCCCCCCCCCC)=O)(=O)O",
#               scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSerine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "OCC(N)C(=O)O",
                       scriptPath = scriptPath)[[1]]
}

# 5. Inositol
# (1) Search inositol
# .searchInositol(smi = "[C@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)(COP(=O)(O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)COC(CCCCCCCCCCCCCCC)=O",
#                 scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchInositol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)C(O)C(O)C(O)C(O)C1(O)",
                       scriptPath = scriptPath)[[1]]
}

