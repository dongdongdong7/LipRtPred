# Calculate CQS Fingerprints
# Barry Song
# 250327

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# Search Phytosphingosine
# .searchPhytosphingosine(smi = "CCCCCCCCCCCCCC[C@H]([C@H]([C@H](CO)N)O)O",
#                         scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPhytosphingosine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(O)C(N)C(O)C(O)CCCCCCCCCCCCC[CH3]",
                       scriptPath = scriptPath)[[1]]
}

# Search Dihydrosphingosine
# .searchDihydrosphingosine(smi = "CCCCCCCCCCCCCCC[C@H]([C@H](CO)N)O",
#                           scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchDihydrosphingosine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(O)C(N)C(O)[CH2]CCCCCCCCCCCCC[CH3]",
                       scriptPath = scriptPath)[[1]]
}

# Search Sphingosine
# .searchSphingosine0(smi = "CCCCCCCCCCCCC/C=C/[C@H]([C@H](CO)N)O",
#                     scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchSphingosine0 <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(O)C(N)C(O)C=CCCCCCCCCCCCC[CH3]",
                       scriptPath = scriptPath)[[1]]
}

# Search Cholesterol
# .searchCholesterol(smi = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
#                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchCholesterol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)CCCC(C)C)CCC4C3CC=C2C1",
                       scriptPath = scriptPath)[[1]]
}

# Search Ergosterol
# .searchErgosterol(smi = "C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
#                   scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchErgosterol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)C=CC(C)C(C)C)CCC4C3=CC=C2C1",
                       scriptPath = scriptPath)[[1]]
}

# Search Brassicasterol
# .searchBrassicasterol(smi = "C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
#                       scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchBrassicasterol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)C=CC(C)C(C)C)CCC4C3CC=C2C1",
                       scriptPath = scriptPath)[[1]]
}

# Search Ethylcholesterol
# alpha
# .searchEthylcholesterol(smi = "CC[C@@H](C(C)C)CC[C@@H](C)[C@@]1([H])CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@@]21C",
#                         scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# beta
# .searchEthylcholesterol(smi = "CC[C@H](C(C)C)CC[C@@H](C)[C@@]1([H])CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@@]21C",
#                         scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchEthylcholesterol <- function(smi, scriptPath){
  posi <- .GetSubstructMatches(smis = smi,
                               SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)CCC(CC)C(C)C)CCC4C3CC=C2C1",
                               scriptPath = scriptPath)[[1]]
  if(length(posi) != 0){
    ch <- .GetAtomCip(smi = smi, atom_idx_vector = unlist(posi)[c(1, 17)], scriptPath = scriptPath)
    if(ch[1] == ch[2]) names(posi) <- "beta"
    else names(posi) <- "alpha"
  }
  return(posi)
}

# Search Stigmasterol
# .searchStigmasterol(smi = "CC[C@H](/C=C/[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C",
#                     scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchStigmasterol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)C=CC(CC)C(C)C)CCC4C3CC=C2C1",
                       scriptPath = scriptPath)[[1]]
}

# Search Cholestanol
# .searchCholestanol(smi = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CC[C@@H](C4)O)C)C",
#                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchCholestanol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C1(O)CCC2(C)C3CCC4(C)C(C(C)CCCC(C)C)CCC4C3CCC2C1",
                       scriptPath = scriptPath)[[1]]
}
