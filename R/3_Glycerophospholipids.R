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

# (2) Search phosphocholine
# .searchPhosphoCholine(smi = "[C@](COP(=O)([O-])OCC[N+](C)(C)C)([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
#                       scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhosphoCholine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)([O-])OCC[N+](C)(C)C",
                       scriptPath = scriptPath)[[1]]
}

# (3) Search phosphonocholine
# .searchPhosphonoCholine(smi = "[C@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)(COP(CC[N+](C)(C)C)(=O)[O-])COC(CCCCCCCCCCCCCCC)=O",
#                         scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhosphonoCholine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)([O-])CC[N+](C)(C)C",
                       scriptPath = scriptPath)[[1]]
}

# (4) Search phosphoethanolamine
# .searchPhosphoEthanolamine(smi = "[C@](COP(=O)(O)OCCN)([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhosphoEthanolamine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)(O)OCC[N;NH2]",
                       scriptPath = scriptPath)[[1]]
}

# (5) Search phosphonoethanolamine
# .searchPhosphonoEthanolamine(smi = "[C@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)(COP(CCN)(=O)O)COC(CCCCCCCCCCCCCCC)=O",
#                              scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhosphonoEthanolamine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)(O)CC[N;NH2]",
                       scriptPath = scriptPath)[[1]]
}

# (6) Search phosphoserine
# .searchPhosphoSerine(smi = "C(O)(=O)[C@@]([H])(N)COP(OC[C@]([H])(OC(CCCCCCCCCCCC)=O)COC(CCCCCCCCCCC)=O)(=O)O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhosphoSerine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)(O)OCC(N)C(=O)O",
                       scriptPath = scriptPath)[[1]]
}

# (7) Search phosphoinositol
# .searchPhosphoInositol(smi = "[C@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)(COP(=O)(O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)COC(CCCCCCCCCCCCCCC)=O",
#                        scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhosphoInositol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)(O)OC1C(O)C(O)C(O)C(O)C1(O)",
                       scriptPath = scriptPath)[[1]]
}

# (8) Search CDP
# .searchCDP(smi = "[C@]([H])(OC(CCCCCCCCCCC)=O)(COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@@H](O)[C@@H](O)[C@H](N2C(=O)N=C(N)C=C2)O1)COC(CCCCCCCCCCC)=O",
#            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchCDP <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)(O)OP(=O)(O)OCC1C(O)C(O)C(n2ccc(N)nc2(=O))O1",
                       scriptPath = scriptPath)[[1]]
}

# (9) Search phosphoethanolamine glycan
# .searchPEG(smi = "N(CCOP(OC[C@]([H])(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCCCCCCCC)(=O)O)C[C@]1(O)[C@@H](O)[C@H](O)[C@H](O)CO1",
#            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPEG <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)(O)OCCNCC1(O)C(O)C(O)C(O)CO1",
                       scriptPath = scriptPath)[[1]]
}

# (10) Search phosphoethanols
# .searchPhosphoEthanols(smi = "[C@](COP(=O)(O)OCC)([H])(OC(CCCCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
#                        scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhosphoEthanols <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)(O)OC[CH3]",
                       scriptPath = scriptPath)[[1]]
}

# (11) Search phosphothreonines
# .searchPhosphoThreonines(smi = "C(O)(=O)[C@@]([H])(N)[C@@H](C)OP(OC[C@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COC(CCCCCCCCCCCCCCCCC)=O)(=O)O",
#                          scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhosphoThreonines <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)(O)OC(C)C(N)C(=O)O",
                       scriptPath = scriptPath)[[1]]
}

#
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
