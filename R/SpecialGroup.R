# Search special group
# Barry Song
# 250324

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# 1. Betaine
# (1) Search betaine
# .searchBetaine(smi = "[C@](COCCC(C(=O)[O-])[N+](C)(C)C)([H])(O)COC(CCCCCCCCCCCCCCCCC)=O",
#                scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchBetaine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[N+](C)(C)(C)CC(=O)[O-]",
                       scriptPath = scriptPath)[[1]]
}

# 2. Glycolysis
# (1) Pentose
# .searchPentose(smi = "C(O)C(O)COP(=O)(O)OC1C(O)C(O)C(CO)O1",
#                scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPentose <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "CC1C([O,N,P,S])C([O,N,P,S])C([O,N,P,S])O1",
                       scriptPath = scriptPath)[[1]]
}
# (2) Hexose
# .searchHexose(smi = "[C@](COP(=O)(O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1)([H])(OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)COC(CCCCCCCCCCCCCCCCC)=O",
#               scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchHexose <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[CH3,CH2,C&$(C-O)&!$(C=O)]C1C([O,N,P,S])C([O,N,P,S])C([O,N,P,S])C([O,N,P,S])O1",
                       scriptPath = scriptPath)[[1]]
}

# 3. Carnitine
# (1) Search carnitine
# .searchCarnitine(smi = "O=C(O[C@](CC([O-])=O)(C[N+](C)(C)C)[H])CCCCC",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchCarnitine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[N+]([CH3])([CH3])([CH3])[CH2][CX4;CH](O)[CH2][CX3](=O)O",
                       scriptPath = scriptPath)[[1]]
}

# 4. Choline
# (1) Search choline
# .searchCholine(smi = "[C@](COP(=O)([O-])OCC[N+](C)(C)C)([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
#                scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchCholine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[N+]([CH3])([CH3])([CH3])[CH2][CH2][O,P]",
                       scriptPath = scriptPath)[[1]]
}

# 5. Ethanolamine
# (1) Search ethanolamine
# .searchEthanolamine(smi = "[C@](COP(=O)(O)OCCN)([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O",
#                     scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchEthanolamine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[OX2,P][CH2][CH2][NH,NH2]",
                       scriptPath = scriptPath)[[1]]
}

# 6. Serine
# (1) Search serine
# .searchSerine(smi = "C(O)(=O)[C@@]([H])(N)COP(OC[C@]([H])(OC(CCCCCCCCCCCC)=O)COC(CCCCCCCCCCC)=O)(=O)O",
#               scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSerine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "O[CH2][CH]([NH2])C(=O)O",
                       scriptPath = scriptPath)[[1]]
}

# 7. Inositol
# (1) Search inositol
# .searchInositol(smi = "[C@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)(COP(=O)(O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)COC(CCCCCCCCCCCCCCC)=O",
#                 scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchInositol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)C(O)C(O)C(O)C(O)C1(O)",
                       scriptPath = scriptPath)[[1]]
}

# 8. Ethanol
# (1) Search Ethanol
.searchEthanol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[OX2]-[CH2]-[CH3]",
                       scriptPath = scriptPath)[[1]]
}

# 9. Threonine
# (1) Search Threonine
.searchThreonine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(O)(C)C([NH2])C(=O)O",
                       scriptPath = scriptPath)[[1]]
}

# 10. Sulfonyl
# (1) Search Sulfonyl
.searchSulfonyl <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "*-S(=O)(=O)-*",
                       scriptPath = scriptPath)[[1]]
}

# 11. Phosphate group
# (1) Search phosphate group
# .searchPhosphate(smi = "C(O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](COP(=O)(O)OCCN)O1)[C@]([H])(OC(CCCCCCC/C=C\\CCCCCC)=O)CO/C=C\\CCCCCCCCCCCCCC",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchPhosphate(smi = "C(CN)P(OC[C@]([H])(NC(=O)CCCCCCCCCCCCC)[C@]([H])(O)/C=C/CCCCCCCCCCCCC)(=O)O",
#                  scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhosphate <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "P(=O)(O)(O)-*",
                       scriptPath = scriptPath)[[1]]
}

# 12. Glucuronic acid
# (1) Search glucuronic acid
# .searchGluAcid(smi = "C(OC(=O)CC=CCC)C(OCCCCC)C(OC1C(O)C(O)C(O)C(C(=O)O)O1)",
#                scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGluAcid <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(=O)([OH])C1C([O,N,P,S])C([O,N,P,S])C([O,N,P,S])C([O,N,P,S])O1",
                       scriptPath = scriptPath)[[1]]
}

# 13. DEMA
# (1) Search Dimethylethanolamine
# .searchDEMA(smi = "[C@](COP(=O)(O)OCCN(C)C)([H])(OC(CCCCCCCCCCC)=O)COC(CCCCCCCCCCC)=O",
#             scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchDEMA <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(O)C[NX3]([CH3])([CH3])",
                       scriptPath = scriptPath)[[1]]
}

# 14. MMEA
# (1) Search N-Methylethanolamine
# .searchMMEA(smi = "[C@](COP(=O)(O)OCCNC)([H])(OC(CCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCC)=O",
#             scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchMMEA <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(O)C[NH][CH3]",
                       scriptPath = scriptPath)[[1]]
}

# 15. Sialic Acid
# .searchSialicAcid(smi = "CC(=O)N[C@@H]1[C@H](CC(O[C@H]1[C@@H]([C@@H](CO)O)O)(C(=O)O)O)O",
#                   scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSialicAcid <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)C(NC(=O)C)C(C(O)C(O)C(O))OC(O)(C(=O)O)C1",
                       scriptPath = scriptPath)[[1]]
}

# 16. Sulfated Galactosyl
# .searchSulfatedGalactosyl(smi = "O=C(OCC(OC(=O)CCCCC=CCC=CCC=CCC=CCC)COC1OC(CS(=O)(=O)O)C(O)C(O)C1O)CCCCCCCCCCCCCCC",
#                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSulfatedGalactosyl <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "S(=O)(=O)(O)CC1C(O)C(O)C(O)C(O)O1",
                       scriptPath = scriptPath)[[1]]
}

# 17. Trimethylhomoserine
# .searchTrimethylhomoserine(smi = "O=C(OCC(O)COCCC(C(=O)[O-])[N+](C)(C)C)CCCCCCCC=CCC=CCC=CCC",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchTrimethylhomoserine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[N+]([CH3])([CH3])([CH3])C(C(=O)O)CCO",
                       scriptPath = scriptPath)[[1]]
}
