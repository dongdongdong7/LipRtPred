# Search substructure of Glycerolipids (GL)
# Barry Song
# 250307

# system.file("python", "molecule_operation.py", package = "LipRtPred")

# smi SMILES string
# scriptPath: path of molecule_operation.py

# 1. Glycero C(O)C(O)C(O)
# (1) Search glycerol
# .searchGlycerol(smi = "OC[C@]([H])(O)COC(CCCCCCCCCCC)=O",
#                 scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerol <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[CH2;$(C-O)](O)-[CH;$(C-O)](O)-[CH2;$(C-O)](O)",
                       scriptPath = scriptPath)[[1]]
}

# (2) Search glycerol ester
# .searchGlycerolEster(smi = "OC[C@]([H])(O)COC(CCCCCCCCCCC)=O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerolEster <- function(smi, scriptPath){
  glycerol_position <- .searchGlycerol(smi = smi, scriptPath = scriptPath)
  O_sn1 <- sapply(glycerol_position, function(x) {x[2]})
  O_sn2 <- sapply(glycerol_position, function(x) {x[4]})
  O_sn3 <- sapply(glycerol_position, function(x) {x[6]})
  acyloxy_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  sn_info <- sapply(acyloxy_position, function(x) {
    if(x[3] %in% O_sn1) return("sn1")
    else if(x[3] %in% O_sn2) return("sn2")
    else if(x[3] %in% O_sn3) return("sn3")
    else return("none")
  })
  acyloxy_position <- acyloxy_position[sn_info != "none"]
  sn_info <- sn_info[sn_info != "none"]
  names(acyloxy_position) <- sn_info
  return(acyloxy_position)
}

# (3) Search glycerol ester chain
# .searchGlycerolEster_Chain(smi = "OC[C@]([H])(O)COC(CCCCCCCCCCC)=O",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerolEster_Chain <- function(smi, scriptPath){
  glycerolEster_position <- .searchGlycerolEster(smi = smi,
                                                 scriptPath = scriptPath)
  lapply(glycerolEster_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# (4) Search glycerol ether
# .searchGlycerolEther(smi = "OC[C@]([H])(O)CO/C=C\\C#C/C=C\\CCCCCCC(C)C",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
# .searchGlycerolEther(smi = "[C@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)(COP(CC[N+](C)(C)C)(=O)[O-])COC(CCCCCCCCCCCCCCC)=O",
#                      scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerolEther <- function(smi, scriptPath){
  glycerol_position <- .searchGlycerol(smi = smi, scriptPath = scriptPath)
  O_sn1 <- sapply(glycerol_position, function(x) {x[2]})
  O_sn2 <- sapply(glycerol_position, function(x) {x[4]})
  O_sn3 <- sapply(glycerol_position, function(x) {x[6]})
  ethers_position <- .GetSubstructMatches(smis = smi,
                                          SMARTS = "[CH2,CH3]-O-[C;!$(C=O)]",
                                          scriptPath = scriptPath)[[1]]
  sn_info <- sapply(ethers_position, function(x) {
    if(x[2] %in% O_sn1) return("sn1")
    else if(x[2] %in% O_sn2) return("sn2")
    else if(x[2] %in% O_sn3) return("sn3")
    else return("none")
  })
  ethers_position <- ethers_position[sn_info != "none"]
  sn_info <- sn_info[sn_info != "none"]
  names(ethers_position) <- sn_info
  return(ethers_position)
}

# (5) Search glycerol ether chain
# .searchGlycerolEther_Chain(smi = "OC[C@]([H])(O)CO/C=C\\C#C/C=C\\CCCCCCC(C)C",
#                            scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchGlycerolEther_Chain <- function(smi, scriptPath){
  glycerolEther_position <- .searchGlycerolEther(smi = smi,
                                                 scriptPath = scriptPath)
  lapply(glycerolEther_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# 2. Dihydroxyacetone (DHA)
# (1) Search dihydroxyacetone
# .searchDihydroxyacetone(smi = "C(=O)(COP(=O)(O)O)COCCCCCCCCCCCCCCCC",
#                         scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchDihydroxyacetone <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(O)C(=O)C(O)",
                       scriptPath = scriptPath)[[1]]
}

# (2) Search Dihydroxyacetone Ester
# .searchDihydroxyacetoneEster(smi = "C(OC(CCCCCCCCCCCCCCC)=O)C(=O)COP(=O)(O)O",
#                              scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchDihydroxyacetoneEster <- function(smi, scriptPath){
  dihydroxyacetone_position <- .searchDihydroxyacetone(smi = smi,
                                                       scriptPath = scriptPath)
  acyloxy_position <- .searchAcyloxy(smi = smi, scriptPath = scriptPath)
  if(length(acyloxy_position) == 0) return(list())
  dihydroxyacetoneEster_logical <- sapply(acyloxy_position, function(x) {
    if(x[3] %in% unlist(dihydroxyacetone_position)) return(TRUE)
    else return(FALSE)
  })
  dihydroxyacetoneEster_position <- acyloxy_position[dihydroxyacetoneEster_logical]
  return(dihydroxyacetoneEster_position)
}

# (3) Search Dihydroxyacetone Ester Chain
# .searchDihydroxyacetoneEster_Chain(smi = "C(OC(CCCCCCCCCCCCCCC)=O)C(=O)COP(=O)(O)O",
#                                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchDihydroxyacetoneEster_Chain <- function(smi, scriptPath){
  dihydroxyacetoneEster_position <- .searchDihydroxyacetoneEster(smi = smi, scriptPath = scriptPath)
  lapply(dihydroxyacetoneEster_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# (4) Search Dihydroxyacetone Ether
# .searchDihydroxyacetoneEther(smi = "C(=O)(COP(=O)(O)O)COCCCCCCCCCCCCCCCC",
#                              scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchDihydroxyacetoneEther <- function(smi, scriptPath){
  dihydroxyacetone_position <- .searchDihydroxyacetone(smi = smi,
                                                       scriptPath = scriptPath)
  dihydroxyacetoneEther_position <- .GetSubstructMatches(smis = smi,
                                                         SMARTS = "[C;!$(C=O)]-O-[C;$(C-C(=O)-C-O)]",
                                                         scriptPath = scriptPath)[[1]]
  if(length(dihydroxyacetoneEther_position) == 0) return(list())
  dihydroxyacetoneEther_logical <- sapply(dihydroxyacetoneEther_position, function(x) {
    if(all(x[c(2,3)] %in% unlist(dihydroxyacetone_position))) return(TRUE)
    else return(FALSE)
  })
  dihydroxyacetoneEther_position <- dihydroxyacetoneEther_position[dihydroxyacetoneEther_logical]
  return(dihydroxyacetoneEther_position)
}

# (5) Search Dihydroxyacetone Ether Chain
# .searchDihydroxyacetoneEther_Chain(smi = "C(=O)(COP(=O)(O)O)COCCCCCCCCCCCCCCCC",
#                                    scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchDihydroxyacetoneEther_Chain <- function(smi, scriptPath){
  dihydroxyacetoneEther_position <- .searchDihydroxyacetoneEther(smi = smi,
                                                                 scriptPath = scriptPath)
  lapply(dihydroxyacetoneEther_position, function(x) {
    .TraverseMolecule(smi = smi,
                      start_atom_idx = x[1],
                      non_traversable_atom_idx = x[-1],
                      scriptPath = scriptPath)
  })
}

# 3. Betaine
# (1) Search betaine
# .searchBetaine(smi = "[C@](COCCC(C(=O)[O-])[N+](C)(C)C)([H])(O)COC(CCCCCCCCCCCCCCCCC)=O",
#                scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchBetaine <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "[N+](C)(C)(C)CC(=O)[O-]",
                       scriptPath = scriptPath)[[1]]
}

# 4. Glycolysis
# (1) Pentose
# .searchPentose(smi = "C(O)C(O)COP(=O)(O)OC1C(O)C(O)C(CO)O1",
#                scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchPentose <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(O)C1C([O,N])C([O,N])C([O,N])O1",
                       scriptPath = scriptPath)[[1]]
}
# (2) Hexose
# .searchHexose(smi = "[C@](COP(=O)(O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1)([H])(OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)COC(CCCCCCCCCCCCCCCCC)=O",
#               scriptPath = system.file("python", "molecule_operation.py", package = "LipRtPred"))
.searchHexose <- function(smi, scriptPath){
  .GetSubstructMatches(smis = smi,
                       SMARTS = "C(O)C1C([O,N])C([O,N])C([O,N])C([O,N])O1",
                       scriptPath = scriptPath)[[1]]
}
